"""
cyp_suite -- post-prediction pipeline for co-folded CYP3A4 complexes.
(c) Byron DeLaBarre / The Consulting Biochemist 2026
One config, five stages, each a function of (config) -> manifest CSV:

    convert    cif  -> pdb            (gemmi)
    refine     pdb  -> minimized      (OpenMM + OpenFF)
    validate   minimized -> pb.csv    (PoseBusters)
    score      *_full.pdb -> dE       (OrbMol-v2)
    select     dE + ipTM + PB -> best (+ plots)

TWO ENVIRONMENTS ARE REQUIRED. `orbmol` holds orb_models and a torch
2.14.1+cu132 build; `openmm` holds openmmforcefields and openff.toolkit.
They cannot be merged (one torch per env), so `score` runs in `orbmol` and
everything else in `openmm`. Every heavy import is therefore INSIDE its
stage function, so this module imports cleanly in either environment.

    # in the openmm env
    python cyp_suite.py --config run.json --stage convert refine validate select
    # in the orbmol env (or let the driver shell out)
    python cyp_suite.py --config run.json --stage score

Paths derive from ONE root plus the run name. Nothing is read from the
environment except as defaults when a config is created.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass, field, asdict, fields
from pathlib import Path

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# configuration
# ---------------------------------------------------------------------------

#: parameters whose values define the refinement; the tag is their hash, so
#: changing any of them writes to a different directory and nothing is
#: silently reused. This is what makes a parameter sweep safe.
REFINE_PARAMS = (
    "backbone_k", "sidechain_k", "ligand_k", "ligand_tol",
    "wall_k", "wall_rc_default", "wall_rc_fe", "wall_rc_pyrrole",
    "fe_restrain", "fe_restrain_max_d0", "fe_window_half", "fe_k_kcal",
    "max_iter", "min_tol_kj", "forcefields", "small_molecule_ff",
)


#: cache for resolved compound -> arm maps, keyed on the fields that define
#: them, so a mixed run does not re-read the classification csv 100 times.
_ARM_CACHE: dict = {}

ARMS = ("binders", "nonbinders")

#: stages whose output depends on the scoring parameters, so a `variant`
#: gets its own copy; the rest are shared with the base run.
VARIANT_STAGES = ("score", "select", "plot")


@dataclass
class RunConfig:
    # --- identity ----------------------------------------------------------
    run_name: str
    engine: str = "esmfold2"            # "esmfold2" | "openfold3"
    arm: str = "nonbinders"             # "binders" | "nonbinders" | "mixed"

    # --- mixed runs: which compounds get which arm -------------------------
    # Set either of these (or class_csv) and the run becomes mixed: one
    # cif_dir, one refine_dir, one best_dir, but each structure refined with
    # its own compound's arm parameters. Explicit lists win over class_csv.
    binders: tuple = ()                 # ("cmpd_1", "cmpd_4", ...)
    nonbinders: tuple = ()
    class_csv: str = ""                 # csv with a cmpd column + class_col
    class_col: str = "arm"              # values "binders"/"nonbinders"
    default_arm: str = "nonbinders"     # for compounds in neither list
    strict_arms: bool = True            # raise if a compound is unassigned

    # --- the single root everything derives from ---------------------------
    predictions_root: str = ""          # e.g. $CYP_PREDICTIONS
    challenge_csv: str = ""
    cif_subdir: str = "all_predictions"
    refine_subdir: str | None = None    # set to reuse an existing directory

    # --- chains and residue names -----------------------------------------
    prot_chain: str = "A"
    heme_chain: str = "B"
    lig_chain: str = "C"
    lig_resname: str = "LIG"
    heme_resnames: tuple = ("HEM", "HEC", "HEB")
    lig_resnames: tuple = ("LIG", "LIG0", "LIG1", "UNK", "UNL")
    pose_regex: str = (r"^(?P<cmpd>cmpd_\d+)_seed_(?P<seed>\d+)"
                       r"_sample_(?P<sample>\d+)_model$")

    # --- refinement (these define the tag) --------------------------------
    backbone_k: float = 100.0           # kcal/mol/A^2
    sidechain_k: float = 2.0
    ligand_k: float = 15.0
    ligand_tol: float = 0.30            # A, flat-bottom free zone
    wall_k: float = 1.0e5               # kJ/mol/nm^4
    # nm. None -> the arm default. A float pins BOTH arms to that value.
    # A dict pins one arm and leaves the other on its default, which is what
    # a clean single-arm ablation needs:
    #     wall_rc_default = {"nonbinders": 0.29}
    wall_rc_default: float | dict | None = None
    wall_rc_fe: float | dict | None = None
    wall_rc_pyrrole: float | dict | None = None
    fe_restrain: bool | None = None     # None -> (arm == "binders")
    stereo_strict: bool = False         # True restores the old behaviour of
                                        # killing a stereo-mismatched pose
    # Reporting only -- NOT in REFINE_PARAMS, so it never changes cfg.tag.
    # Stiff positional restraints put a systematic floor under the residual
    # force (a constant ~70 kJ/mol/nm says the restrained protein, not the
    # pose), so warning per structure against min_tol_kj fires on everything
    # and hides the outliers that do matter.
    warn_force_kj: float = 250.0
    fe_restrain_max_d0: float = 2.35    # A
    fe_window_half: float = 0.10        # A
    fe_k_kcal: float = 20.0
    max_iter: int = 5000
    min_tol_kj: float = 2.0             # kJ/mol/nm
    forcefields: tuple = ("amber/ff14SB.xml", "implicit/obc2.xml")
    small_molecule_ff: str = "gaff-2.11"
    nagl_model: str = "openff-gnn-am1bcc-1.0.0.pt"

    # --- scoring -----------------------------------------------------------
    shell: float = 3.5                  # A around ligand heavy atoms
    max_atoms: int = 360                # OMol25 training ceiling
    score_device: str = "cuda"

    # --- validation --------------------------------------------------------
    pb_config: str = "dock"
    pb_critical: tuple = (
        "minimum_distance_to_organic_cofactors", "bond_angles",
        "internal_steric_clash", "bond_lengths", "sanitization",
        "all_atoms_connected", "double_bond_stereochemistry",
        "tetrahedral_chirality", "minimum_distance_to_protein",
        "minimum_distance_to_inorganic_cofactors",
        "volume_overlap_with_protein",
        "volume_overlap_with_organic_cofactors")

    # --- selection ---------------------------------------------------------
    de_tol: float = 0.20                # eV; within this of min = undecided
    drop_positive_de: bool = True
    require_pb_pass: bool = True        # walk down the ranking until one passes
    prefer_stereo_ok: bool = True       # rank stereo-correct poses first, but
                                        # fall back rather than submit nothing
    iptm_csv: str = ""                  # defaults to run_dir/harvest_iptm.csv

    # --- scoring variant ----------------------------------------------------
    # A label for a SCORING/SELECTION variant that reuses an existing
    # refine_dir. `shell` is not in REFINE_PARAMS -- changing it must not
    # force a re-refinement -- so it does not change cfg.tag, and without
    # this the new score.csv and best_poses would overwrite the originals.
    # With it set, score/select/plot write to their own directories while
    # convert/refine/validate manifests stay shared, since a PoseBusters
    # verdict depends on the structure, not on the scoring shell.
    variant: str = ""

    # --- environments ------------------------------------------------------
    env_main: str = "openmm"
    env_score: str = "orbmol"

    # -- arm resolution -----------------------------------------------------
    @property
    def is_mixed(self) -> bool:
        return (self.arm == "mixed"
                or bool(self.binders and self.nonbinders)
                or bool(self.class_csv))

    def _id_maps(self):
        """-> (set of valid cmpd_N, {OCNT-N: cmpd_N}) from the challenge csv.

        Empty when challenge_csv is unset, in which case names are normalised
        by shape but not checked against anything.
        """
        if not self.challenge_csv or not Path(self.challenge_csv).is_file():
            return set(), {}
        key = ("_ids", self.challenge_csv)
        if key in _ARM_CACHE:
            return _ARM_CACHE[key]
        df = pd.read_csv(self.challenge_csv)
        id2 = df["id2"].astype(str).str.strip()
        id1 = df["id"].astype(str).str.strip()
        out = (set(id2), {o.upper(): c for o, c in zip(id1, id2)})
        _ARM_CACHE[key] = out
        return out

    def norm_cmpd(self, tok, where="arm list") -> str:
        """Any spelling of a compound id -> the canonical 'cmpd_N'.

        Accepts cmpd_12, cmpd12, CMPD-12, 12, OCNT-2313038, ocnt_2313038.
        Format is forgiving; existence is not. A token that does not resolve
        to a compound in the challenge csv raises here, where the message can
        name it, rather than surfacing later as some *other* compound being
        unassigned.
        """
        s = str(tok).strip()
        if not s:
            raise ValueError(f"{where}: empty compound id")
        known, ocnt2cmpd = self._id_maps()

        m = re.fullmatch(r"(?i)ocnt[-_ ]?(\d+)", s)
        if m:
            cand = f"OCNT-{m.group(1)}"
            if ocnt2cmpd:
                if cand.upper() not in ocnt2cmpd:
                    raise KeyError(
                        f"{where}: {s!r} is not in {Path(self.challenge_csv).name} "
                        f"(looked for id {cand})")
                return ocnt2cmpd[cand.upper()]
            raise KeyError(f"{where}: cannot resolve {s!r} without "
                           f"challenge_csv set")
        m = re.fullmatch(r"(?i)(?:cmpd|compound|cpd)?[-_ ]?(\d+)", s)
        if not m:
            raise ValueError(
                f"{where}: cannot parse {s!r} as a compound id "
                f"(want cmpd_12, 12, or OCNT-2313038)")
        cand = f"cmpd_{int(m.group(1))}"
        if known and cand not in known:
            raise KeyError(
                f"{where}: {s!r} -> {cand}, which is not in "
                f"{Path(self.challenge_csv).name}")
        return cand

    def arm_map(self) -> dict:
        """-> {'cmpd_1': 'binders', ...}. Explicit lists override class_csv.

        Keys are canonical cmpd_N whatever spelling went in, so the two lists
        may be written in OCNT ids, cmpd ids, or a mix of both.
        """
        key = (self.class_csv, self.class_col, self.challenge_csv,
               tuple(self.binders), tuple(self.nonbinders))
        if key in _ARM_CACHE:
            return _ARM_CACHE[key]
        m = {}
        if self.class_csv:
            df = pd.read_csv(self.class_csv)
            ccol = next((c for c in ("cmpd", "id2", "compound", "id")
                         if c in df.columns), None)
            if ccol is None or self.class_col not in df.columns:
                raise KeyError(
                    f"{self.class_csv}: need a cmpd/id2/id column and "
                    f"{self.class_col!r}; found {list(df.columns)}")
            for c, a in zip(df[ccol].astype(str).str.strip(),
                            df[self.class_col].astype(str).str.strip()):
                a = a.strip().lower().replace("_", "-")
                if a in ("binder", "binders", "heme", "heme-binder", "type-ii",
                         "true", "yes", "y", "1"):
                    a = "binders"
                elif a in ("nonbinder", "nonbinders", "non-binder",
                           "non-binders", "false", "no", "n", "0"):
                    a = "nonbinders"
                if a not in ARMS:
                    raise ValueError(
                        f"{self.class_csv}: bad {self.class_col} value {a!r} "
                        f"for {c}; want binders/nonbinders (or true/false)")
                m[self.norm_cmpd(c, where=Path(self.class_csv).name)] = a
        # explicit lists win over class_csv; a compound in BOTH lists is an
        # error, not a precedence question -- normalising means cmpd_12 and
        # OCNT-2313038 in opposite lists is the same collision
        nb = {self.norm_cmpd(c, where="nonbinders") for c in self.nonbinders}
        bd = {self.norm_cmpd(c, where="binders") for c in self.binders}
        both = nb & bd
        if both:
            raise ValueError(
                f"in both binders and nonbinders: {sorted(both)}")
        m.update({c: "nonbinders" for c in nb})
        m.update({c: "binders" for c in bd})
        _ARM_CACHE[key] = m
        return m

    def arm_of(self, cmpd: str) -> str:
        """The arm to refine this compound with. Takes cmpd_N or OCNT-N."""
        if not self.is_mixed:
            return self.arm
        try:
            k = self.norm_cmpd(cmpd, where="arm_of")
        except (KeyError, ValueError):
            k = str(cmpd).strip()
        a = self.arm_map().get(k)
        if a is None:
            if self.strict_arms:
                raise KeyError(
                    f"{cmpd} is in neither binders nor nonbinders. Add it, "
                    f"or set strict_arms=False to fall back to "
                    f"default_arm={self.default_arm!r}.")
            return self.default_arm
        return a

    def arm_table(self) -> pd.DataFrame:
        """The resolved assignment, with both ids, for eyeballing."""
        _, ocnt2cmpd = self._id_maps()
        back = {v: k for k, v in ocnt2cmpd.items()}
        m = self.arm_map()
        return pd.DataFrame(
            [{"cmpd": c, "ocnt": back.get(c, ""), "arm": a}
             for c, a in sorted(m.items(),
                                key=lambda kv: int(kv[0].split("_")[-1]))])

    # -- arm-dependent refinement parameters --------------------------------
    @staticmethod
    def _override(value, arm):
        """float -> both arms; dict -> that arm only; None/absent -> default."""
        if value is None:
            return None
        if isinstance(value, dict):
            return value.get(arm)
        return value

    def wall_rc_default_for(self, arm=None) -> float:
        arm = arm or self.arm
        v = self._override(self.wall_rc_default, arm)
        if v is not None:
            return float(v)
        return 0.34 if arm == "nonbinders" else 0.26

    def wall_rc_fe_for(self, arm=None) -> float:
        arm = arm or self.arm
        v = self._override(self.wall_rc_fe, arm)
        if v is not None:
            return float(v)
        return 0.19 if arm == "binders" else self.wall_rc_default_for(arm)

    def wall_rc_pyrrole_for(self, arm=None) -> float:
        arm = arm or self.arm
        v = self._override(self.wall_rc_pyrrole, arm)
        if v is not None:
            return float(v)
        return 0.24 if arm == "binders" else self.wall_rc_default_for(arm)

    def fe_restrain_for(self, arm=None) -> bool:
        arm = arm or self.arm
        v = self._override(self.fe_restrain, arm)
        if v is not None:
            return bool(v)
        return arm == "binders"

    # back-compat: the single-arm properties the earlier cells use
    @property
    def arm_wall_rc_default(self) -> float:
        return self.wall_rc_default_for()

    @property
    def arm_wall_rc_fe(self) -> float:
        return self.wall_rc_fe_for()

    @property
    def arm_wall_rc_pyrrole(self) -> float:
        return self.wall_rc_pyrrole_for()

    @property
    def arm_fe_restrain(self) -> bool:
        return self.fe_restrain_for()

    @property
    def tag(self) -> str:
        """Short hash of the refinement parameters.

        In a mixed run the compound->arm assignment is part of the hash:
        moving a compound between arms changes what refinement it gets, so it
        must not silently reuse the old directory. Single-arm tags are
        unchanged by this, so existing refine_dir names stay valid.
        """
        blob = {k: _jsonable(getattr(self, k)) for k in REFINE_PARAMS}
        if self.is_mixed:
            blob["_arms"] = sorted(self.arm_map().items())
        return hashlib.sha1(
            json.dumps(blob, sort_keys=True).encode()).hexdigest()[:8]

    @property
    def run_dir(self) -> Path:
        return Path(self.predictions_root) / self.run_name

    @property
    def cif_dir(self) -> Path:
        return self.run_dir / self.cif_subdir

    @property
    def pdb_dir(self) -> Path:
        return self.cif_dir / "pdb_format"

    @property
    def work_dir(self) -> Path:
        return self.run_dir / "min_work"

    @property
    def refine_dir(self) -> Path:
        if self.refine_subdir:
            return self.run_dir / self.refine_subdir
        return self.run_dir / f"refined_{self.arm_label}_{self.tag}"

    @property
    def arm_label(self) -> str:
        return "mixed" if self.is_mixed else self.arm

    @property
    def base_stage_dir(self) -> Path:
        """Where the shared (refinement-level) manifests live."""
        return self.refine_dir / "_suite"

    @property
    def stage_dir(self) -> Path:
        return self.refine_dir / (f"_suite_{self.variant}" if self.variant
                                  else "_suite")

    @property
    def best_dir(self) -> Path:
        return self.refine_dir / (f"best_poses_{self.variant}" if self.variant
                                  else "best_poses")

    def ensure_dirs(self) -> "RunConfig":
        """Create every directory a stage may write to.

        base_stage_dir matters even when a variant is set: the shared
        manifests (convert/refine/validate) live there, so creating only
        stage_dir leaves refine.csv with nowhere to go on a fresh refine_dir.
        """
        for d in (self.refine_dir, self.base_stage_dir, self.stage_dir,
                  self.work_dir):
            d.mkdir(parents=True, exist_ok=True)
        return self

    # manifests -- one per stage, all keyed on (cmpd, seed, sample)
    def manifest(self, name: str) -> Path:
        """Variant owns score/select/plot; everything else is shared.

        convert/refine/validate describe the structures in refine_dir and do
        not depend on the scoring parameters, so a variant reads and appends
        to the originals instead of recomputing them.
        """
        if self.variant and name in VARIANT_STAGES:
            return self.stage_dir / f"{name}.csv"
        return self.base_stage_dir / f"{name}.csv"

    # -- io -----------------------------------------------------------------
    def save(self, path: str | Path) -> Path:
        path = Path(path)
        path.parent.mkdir(parents=True, exist_ok=True)
        d = {k: _jsonable(v) for k, v in asdict(self).items()}
        d["_derived"] = {"tag": self.tag, "run_dir": str(self.run_dir),
                         "refine_dir": str(self.refine_dir)}
        path.write_text(json.dumps(d, indent=1))
        return path

    @classmethod
    def load(cls, path: str | Path) -> "RunConfig":
        d = json.loads(Path(path).read_text())
        d.pop("_derived", None)
        names = {f.name for f in fields(cls)}
        tuple_fields = {f.name for f in fields(cls)
                        if getattr(f.type, "__name__", str(f.type)) == "tuple"
                        or f.type is tuple}
        kw = {}
        for k, v in d.items():
            if k not in names:
                continue
            kw[k] = tuple(v) if (k in tuple_fields and isinstance(v, list)) else v
        return cls(**kw)

    @classmethod
    def from_env(cls, run_name: str, **over) -> "RunConfig":
        """Defaults from the environment, for creating a config the first time."""
        try:
            from dotenv import load_dotenv, find_dotenv
            load_dotenv(find_dotenv(), override=True)
        except Exception:
            pass
        kw = dict(
            run_name=run_name,
            predictions_root=os.getenv("CYP_PREDICTIONS", ""),
            challenge_csv=os.getenv("CHALLENGE_CSV", ""),
        )
        kw.update(over)
        return cls(**kw)

    def describe(self) -> str:
        L = [f"run        {self.run_name}   engine={self.engine}  "
             f"arm={self.arm_label}",
             f"tag        {self.tag}  (refinement parameter hash"
             + (" + arm assignment)" if self.is_mixed else ")"),
             f"cif        {self.cif_dir}",
             f"pdb        {self.pdb_dir}",
             f"refined    {self.refine_dir}",
             f"best       {self.best_dir}",
             f"stage      {self.stage_dir}"
             + (f"   (variant {self.variant!r}; refine/validate shared with "
                f"{self.base_stage_dir.name})" if self.variant else ""),
             f"ligand_k   {self.ligand_k}   tol {self.ligand_tol} A"]
        arms = ARMS if self.is_mixed else (self.arm,)
        if self.is_mixed:
            m = self.arm_map()
            n_b = sum(1 for v in m.values() if v == "binders")
            L.append(f"arms       {n_b} binders / {len(m)-n_b} nonbinders"
                     f"   (strict={self.strict_arms})")
        for a in arms:
            L.append(f"  {a:<10} wall k={self.wall_k:g} "
                     f"rc={self.wall_rc_default_for(a)} "
                     f"fe={self.wall_rc_fe_for(a)} "
                     f"pyr={self.wall_rc_pyrrole_for(a)} nm   "
                     f"fe_restr={self.fe_restrain_for(a)}")
        L.append(f"shell      {self.shell} A   max_atoms {self.max_atoms}")
        return "\n".join(L)


def _jsonable(v):
    if isinstance(v, Path):
        return str(v)
    if isinstance(v, tuple):
        return list(v)
    if isinstance(v, (np.floating, np.integer)):
        return v.item()
    return v


# ---------------------------------------------------------------------------
# shared helpers -- gemmi only, so usable in both environments
# ---------------------------------------------------------------------------

AA3 = set("ALA ARG ASN ASP CYS CYX GLN GLU GLY HIS HID HIE HIP ILE LEU LYS "
          "MET PHE PRO SER THR TRP TYR VAL".split())


def _log(msg, flush=True):
    print(msg, flush=flush)


def parse_pose(path, cfg: RunConfig):
    """-> {'cmpd','seed','sample'} or None. Tolerates the _full suffix."""
    stem = Path(path).stem
    if stem.endswith("_full"):
        stem = stem[:-5]
    m = re.match(cfg.pose_regex, stem)
    return m.groupdict() if m else None


def challenge_maps(cfg: RunConfig):
    """-> (smiles_by_id2, ocnt_by_id2). id2 holds cmpd_<N>."""
    df = pd.read_csv(cfg.challenge_csv)
    for col in ("id2", "id", "SMILES"):
        if col not in df.columns:
            raise KeyError(f"no {col!r} in {cfg.challenge_csv}; "
                           f"found {list(df.columns)}")
    key = df["id2"].astype(str).str.strip()
    return (dict(zip(key, df["SMILES"])),
            dict(zip(key, df["id"].astype(str).str.strip())))


def _chain_xyz(pdb_path, chain_name, heavy_only=True):
    import gemmi
    st = gemmi.read_structure(str(pdb_path)); st.setup_entities()
    return np.array([[a.pos.x, a.pos.y, a.pos.z]
                     for ch in st[0] if ch.name == chain_name
                     for r in ch for a in r
                     if not (heavy_only and a.element.name == "H")]
                    ).reshape(-1, 3)


def _min_dist(a, b):
    if a.size == 0 or b.size == 0:
        return np.nan
    return float(np.sqrt(((a[:, None, :] - b[None, :, :]) ** 2).sum(-1)).min())


def min_lig_contact(p, cfg, heavy_only=False):
    """Minimum protein-ligand distance. heavy_only=False includes hydrogens,
    which is the check that caught Modeller placing protein H on the ligand."""
    return _min_dist(_chain_xyz(p, cfg.prot_chain, heavy_only),
                     _chain_xyz(p, cfg.lig_chain, heavy_only))


def min_lig_heme(p, cfg):
    return _min_dist(_chain_xyz(p, cfg.heme_chain),
                     _chain_xyz(p, cfg.lig_chain))


def lig_rmsd(a_pdb, b_pdb, cfg):
    xa = _chain_xyz(a_pdb, cfg.lig_chain)
    xb = _chain_xyz(b_pdb, cfg.lig_chain)
    if xa.shape != xb.shape or xa.size == 0:
        return np.nan
    return float(np.sqrt(((xa - xb) ** 2).sum(1).mean()))


def heme_penetration(pdb_path, cfg):
    """Is a ligand heavy atom inside the porphyrin macrocycle?

    The four pyrrole N define the ring plane; the Fe->atom vector splits into
    axial and in-plane parts. Near-zero axial with small lateral means the
    atom is THROUGH the ring.
    """
    import gemmi
    st = gemmi.read_structure(str(pdb_path)); st.setup_entities()
    N4, fe, lig = [], None, []
    for ch in st[0]:
        for r in ch:
            rn = r.name.strip().upper()
            if rn in cfg.heme_resnames:
                for a in r:
                    if a.element.name.upper() == "FE":
                        fe = np.array(a.pos.tolist())
                    elif a.name.strip() in ("NA", "NB", "NC", "ND"):
                        N4.append(a.pos.tolist())
            elif rn in cfg.lig_resnames:
                lig += [a.pos.tolist() for a in r if a.element.name != "H"]
    if fe is None or len(N4) != 4 or not lig:
        return {"min_fe_dist": np.nan, "through_ring": np.nan}
    N4, lig = np.array(N4), np.array(lig)
    normal = np.linalg.svd(N4 - N4.mean(0))[2][-1]
    v = lig - fe
    axial = v @ normal
    lateral = np.linalg.norm(v - np.outer(axial, normal), axis=1)
    return {"min_fe_dist": float(np.linalg.norm(v, axis=1).min()),
            "through_ring": int(((np.abs(axial) < 1.2) & (lateral < 3.2)).sum())}


def read_heme(pdb_path, cfg):
    """-> (list[(ELEMENT, atom_name)], positions in Angstrom)."""
    import gemmi
    st = gemmi.read_structure(str(pdb_path)); st.setup_entities()
    atoms, pos = [], []
    for ch in st[0]:
        for r in ch:
            if r.name.strip().upper() not in cfg.heme_resnames:
                continue
            for a in r:
                atoms.append((a.element.name.upper(), a.name.strip()))
                pos.append(a.pos.tolist())
    return atoms, np.array(pos).reshape(-1, 3)


# ---------------------------------------------------------------------------
# stage 1: convert
# ---------------------------------------------------------------------------

def _insert_remarks(pdb_path, lines, remark_no=999):
    txt = open(pdb_path).readlines()
    cut = next((i for i, l in enumerate(txt)
                if l.startswith(("MODEL", "ATOM", "HETATM"))), len(txt))
    with open(pdb_path, "w") as fh:
        fh.writelines(txt[:cut] + [l.rstrip() + "\n" for l in lines] + txt[cut:])


def _cif_to_pdb(cif_path, out_pdb, cfg, report=False):
    """mmCIF -> PDB with chains prot/heme/lig and ONE ligand residue.

    A co-folding engine names a SMILES ligand itself and may split it across
    residues, and cif residue names can exceed the 3 characters PDB allows
    (which the competition validator rejects). Everything that is neither
    protein, heme nor water is merged into one residue.
    """
    import gemmi
    st = gemmi.read_structure(str(cif_path)); st.setup_entities()
    prot = gemmi.Chain(cfg.prot_chain)
    heme = gemmi.Residue(); heme.name = "HEM"
    heme.seqid = gemmi.SeqId("1"); heme.het_flag = "H"
    lig = gemmi.Residue(); lig.name = cfg.lig_resname
    lig.seqid = gemmi.SeqId("1"); lig.het_flag = "H"
    seen = []

    for ch in st[0]:
        for r in ch:
            nm = r.name.strip().upper()
            seen.append((ch.name, nm, len(r)))
            if nm in AA3:
                prot.add_residue(r)
            elif nm in cfg.heme_resnames:
                for a in r:
                    heme.add_atom(a)
            elif nm in ("HOH", "WAT"):
                continue
            else:
                for a in r:
                    lig.add_atom(a)

    if report:
        _log(f"    non-protein source residues: "
             f"{[(c, n, k) for c, n, k in seen if n not in AA3]}")
        _log(f"    -> protein {sum(len(r) for r in prot)} atoms, "
             f"HEM {len(heme)}, {cfg.lig_resname} {len(lig)}")
    if not len(lig):
        raise ValueError("no ligand residue found")

    new = gemmi.Structure(); new.name = Path(out_pdb).stem
    model = gemmi.Model("1")
    if len(prot):
        model.add_chain(prot)
    for cname, res in ((cfg.heme_chain, heme), (cfg.lig_chain, lig)):
        if len(res):
            c = gemmi.Chain(cname); c.add_residue(res); model.add_chain(c)
    new.add_model(model); new.setup_entities()
    new.write_pdb(str(out_pdb))
    return sum(len(r) for r in prot), len(heme), len(lig)


def stage_convert(cfg: RunConfig, force=False) -> pd.DataFrame:
    import gemmi
    cfg.pdb_dir.mkdir(parents=True, exist_ok=True)
    cfg.ensure_dirs()
    files = sorted(cfg.cif_dir.glob("*.cif"))
    _log(f"[convert] {len(files)} cif in {cfg.cif_dir}")
    if not files:
        raise FileNotFoundError(f"no .cif files in {cfg.cif_dir}")

    rows, failed = [], []
    for i, cif in enumerate(files, 1):
        out = cfg.pdb_dir / f"{cif.stem}.pdb"
        info = parse_pose(cif, cfg) or {}
        try:
            if out.is_file() and out.stat().st_size > 0 and not force:
                n_p = n_h = n_l = np.nan
            elif cfg.engine in ("esmfold2", "openfold3"):
                # both engines now go through the same gemmi path; the old
                # OF3 route only renamed LIG0 -> LIG, which this subsumes
                n_p, n_h, n_l = _cif_to_pdb(cif, out, cfg, report=(i == 1))
            else:
                raise ValueError(f"unknown engine {cfg.engine!r}")
            if not out.read_text().startswith("REMARK 999"):
                _insert_remarks(out, [
                    f"REMARK 999 RUN_NAME {cfg.run_name}",
                    f"REMARK 999 ENGINE   {cfg.engine}",
                    f"REMARK 999 SOURCE   {cif.name}"])
            rows.append({**info, "cif": cif.name, "pdb": out.name,
                         "path": str(out), "n_protein": n_p,
                         "n_heme": n_h, "n_lig": n_l})
        except Exception as ex:
            failed.append((cif.name, str(ex)))
            _log(f"  [!] {cif.name}: {ex}")

    df = pd.DataFrame(rows)
    df.to_csv(cfg.manifest("convert"), index=False)

    bad = []
    for p in sorted(cfg.pdb_dir.glob("*.pdb")):
        st = gemmi.read_structure(str(p)); st.setup_entities()
        ids = {ch.name for ch in st[0]}
        names = {r.name.strip() for ch in st[0] for r in ch}
        long_ = {n for n in names if len(n) > 3}
        need = {cfg.prot_chain, cfg.heme_chain, cfg.lig_chain}
        if not ids >= need or long_ or "HEM" not in names:
            bad.append((p.name, sorted(ids), sorted(long_), "HEM" in names))
    _log(f"[convert] {len(df)} written/present, {len(failed)} failed, "
         f"{len(bad)} with layout problems")
    for n, ids, long_, hem in bad:
        _log(f"  [!] {n}: chains={ids} over-3-char={long_ or '-'} heme={hem}")
    return df


# ---------------------------------------------------------------------------
# stage 2: refine  (OpenMM + OpenFF -- env_main)
# ---------------------------------------------------------------------------

def _refine_helpers(cfg):
    """Build the OpenMM-dependent closures. Imported here, not at module top,
    so this module loads in the orbmol env too."""
    import gemmi
    import openmm as mm
    from openmm import app, unit
    from openmmforcefields.generators import SystemGenerator
    from openff.toolkit.topology import Molecule
    from openff.toolkit.utils.nagl_wrapper import NAGLToolkitWrapper
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from pdbfixer import PDBFixer

    NAGL = NAGLToolkitWrapper()
    WALL_GROUP = 31
    PYRROLE = ("NA", "NB", "NC", "ND")
    DONOR_ELEMS = ("N", "O", "S")

    def _write_subset(src, dst, keep_chains):
        st = gemmi.read_structure(str(src)); st.setup_entities()
        for model in st:
            for cname in [c.name for c in model if c.name not in keep_chains]:
                model.remove_chain(cname)
        st.write_pdb(str(dst))

    def _fix_protein(prot_pdb, fixed_pdb):
        """PDBFixer: complete sidechains, add OXT, protonate at pH 7.4.

        The heme-ligating Cys is deprotonated in the written output, not here:
        renaming it CYM breaks Topology.createStandardBonds(), which is keyed
        on residue name and would sever the backbone at that position.
        """
        fx = PDBFixer(filename=str(prot_pdb))
        fx.findMissingResidues()
        fx.missingResidues = {}        # keep the prediction; build no loops
        fx.findMissingAtoms()
        fx.addMissingAtoms()
        fx.addMissingHydrogens(7.4)
        with open(str(fixed_pdb), "w") as fh:
            app.PDBFile.writeFile(fx.topology, fx.positions, fh)
        return fixed_pdb

    def _axial_cys_from_fe(st, fe_pos, max_d=3.2):
        best, bd = None, np.inf
        for ch in st[0]:
            for r in ch:
                if r.name.strip().upper() not in ("CYS", "CYX", "CYM"):
                    continue
                sg = next((a for a in r if a.name.strip() == "SG"), None)
                if sg is None:
                    continue
                d = float(np.linalg.norm(np.array(sg.pos.tolist()) - fe_pos))
                if d < bd:
                    best, bd = (ch.name, r.seqid.num), d
        return (best, bd) if best and bd <= max_d else (None, bd)

    def strip_cys_hg_pdb(pdb_path, fe_pos, max_d=3.2, verbose=False):
        """Drop the HG from the heme-ligating cysteine, by text edit.

        The residue keeps the name CYS on purpose -- see _fix_protein.
        """
        if fe_pos is None:
            return 0
        st = gemmi.read_structure(str(pdb_path)); st.setup_entities()
        key, d = _axial_cys_from_fe(st, np.asarray(fe_pos), max_d)
        if key is None:
            return 0
        chain_id, seq = key
        kept, dropped = [], 0
        for line in open(pdb_path):
            if line.startswith(("ATOM", "HETATM")):
                if (line[17:20].strip().upper() in ("CYS", "CYX", "CYM")
                        and line[21] == chain_id
                        and line[22:26].strip() == str(seq)
                        and line[12:16].strip() in ("HG", "HG1")):
                    dropped += 1
                    continue
            kept.append(line)
        with open(pdb_path, "w") as fh:
            fh.writelines(kept)
        if verbose:
            _log(f"      Cys{seq} -> thiolate: {dropped} HG removed "
                 f"(SG-Fe {d:.2f} A)")
        return dropped

    def _gemmi_atom(name, sym, pos_ang):
        at = gemmi.Atom()
        at.name = name[:4]
        at.element = gemmi.Element(sym or "X")
        at.pos = gemmi.Position(*(float(c) for c in pos_ang))
        at.occ, at.b_iso = 1.0, 0.0
        return at

    def _write_full_complex(topology, positions, out_path,
                            heme_atoms=(), heme_xyz_ang=None):
        """Explicit chains, ONE ligand residue with its hydrogens.

        PDBFile.writeFile split the ligand into a heavy-atom residue and a
        separate all-hydrogen UNK residue in another chain. The heme is an
        external field, not part of the system, so its (frozen, hence
        unchanged) input coordinates are passed in.
        """
        xyz = np.asarray(positions.value_in_unit(unit.angstrom))
        prot = gemmi.Chain(cfg.prot_chain)
        lig_res = gemmi.Residue()
        lig_res.name = cfg.lig_resname
        lig_res.seqid = gemmi.SeqId("1"); lig_res.het_flag = "H"

        for res in topology.residues():
            rn = res.name.strip().upper()
            if rn in AA3:
                gr = gemmi.Residue()
                gr.name = rn
                gr.seqid = gemmi.SeqId(str(res.id))
                gr.het_flag = "A"
                for a in res.atoms():
                    gr.add_atom(_gemmi_atom(
                        a.name, a.element.symbol if a.element else "X",
                        xyz[a.index]))
                prot.add_residue(gr)
            elif rn not in cfg.heme_resnames:
                for a in res.atoms():
                    lig_res.add_atom(_gemmi_atom(
                        a.name, a.element.symbol if a.element else "X",
                        xyz[a.index]))

        heme_res = gemmi.Residue()
        heme_res.name = "HEM"
        heme_res.seqid = gemmi.SeqId("1"); heme_res.het_flag = "H"
        if heme_xyz_ang is not None and len(heme_atoms):
            for (sym, name), p in zip(heme_atoms, np.asarray(heme_xyz_ang)):
                heme_res.add_atom(_gemmi_atom(name, sym.capitalize(), p))

        st = gemmi.Structure(); st.name = Path(out_path).stem
        model = gemmi.Model("1")
        if len(prot):
            model.add_chain(prot)
        for cname, res in ((cfg.heme_chain, heme_res), (cfg.lig_chain, lig_res)):
            if len(res):
                c = gemmi.Chain(cname); c.add_residue(res); model.add_chain(c)
        st.add_model(model); st.setup_entities()
        st.write_pdb(str(out_path))
        return len(lig_res), len(heme_res)

    def _strip_ligand_h(full_pdb, out_pdb):
        """Submission format: same chains, ligand heavy atoms only."""
        st = gemmi.read_structure(str(full_pdb)); st.setup_entities()
        for ch in st[0]:
            if ch.name != cfg.lig_chain:
                continue
            for r in ch:
                for nm in [a.name for a in r if a.element.name == "H"]:
                    while True:
                        i = next((k for k, a in enumerate(r)
                                  if a.name == nm), None)
                        if i is None:
                            break
                        del r[i]
        st.setup_entities()
        st.write_pdb(str(out_pdb))

    def _lig_indices(topology):
        return [a.index for a in topology.atoms()
                if a.residue.name.strip().upper() in cfg.lig_resnames]

    def _nearest_donor_to_fe(topology, positions, heme_atoms, heme_xyz):
        fe = next((np.array(p) for (s, _n), p in zip(heme_atoms, heme_xyz)
                   if s == "FE"), None)
        if fe is None:
            return None, np.nan, None
        xyz = np.asarray(positions.value_in_unit(unit.angstrom))
        cand = [a for a in topology.atoms()
                if a.residue.name.strip().upper() in cfg.lig_resnames
                and a.element is not None
                and a.element.symbol.upper() in DONOR_ELEMS]
        if not cand:
            return None, np.nan, None
        d = [float(np.linalg.norm(xyz[a.index] - fe)) for a in cand]
        j = int(np.argmin(d))
        return cand[j].index, d[j], cand[j].element.symbol

    def _wall_radii(heme_atoms, arm):
        rc = []
        for sym, name in heme_atoms:
            if sym == "FE":
                rc.append(cfg.wall_rc_fe_for(arm))
            elif name in PYRROLE:
                rc.append(cfg.wall_rc_pyrrole_for(arm))
            else:
                rc.append(cfg.wall_rc_default_for(arm))
        return np.array(rc)

    def _add_heme_wall(system, topology, heme_atoms, heme_xyz, arm):
        """Soft repulsion from each heme atom, on ligand atoms only.

        E = k*(rc^2 - r^2)^2 for r < rc. A polynomial, bounded by k*rc^4 even
        at exact overlap, with no 1/r anywhere, so it cannot diverge. rc is
        per-centre so the iron can be permissive (binders) or fully steric
        (non-binders) without touching anything else.
        """
        lig_idx = _lig_indices(topology)
        if not lig_idx:
            raise RuntimeError("no ligand atoms in the topology")
        heme_nm = np.asarray(heme_xyz) / 10.0
        rc_i = _wall_radii(heme_atoms, arm)
        wall = mm.CustomExternalForce(
            "step(rci - rh)*k*(rci^2 - rh^2)^2;"
            "rh = sqrt((x-hx)^2 + (y-hy)^2 + (z-hz)^2)")
        wall.addGlobalParameter("k", cfg.wall_k)
        for p in ("hx", "hy", "hz", "rci"):
            wall.addPerParticleParameter(p)
        for i in lig_idx:
            for (hx, hy, hz), rc in zip(heme_nm, rc_i):
                wall.addParticle(i, [float(hx), float(hy), float(hz),
                                     float(rc)])
        wall.setForceGroup(WALL_GROUP)
        system.addForce(wall)
        return lig_idx, heme_nm, rc_i

    def _wall_energy_numpy(lig_xyz_nm, heme_nm, rc_i):
        d = np.linalg.norm(np.asarray(lig_xyz_nm)[:, None, :]
                           - np.asarray(heme_nm)[None, :, :], axis=-1)
        lim = np.asarray(rc_i)[None, :]
        inside = d < lim
        return float((cfg.wall_k * (lim ** 2 - d ** 2)[inside] ** 2).sum())

    def _add_fe_restraint(system, topology, positions, heme_atoms, heme_xyz):
        """Flat-bottom restraint holding an EXISTING Fe-donor contact.

        ff14SB/GAFF have no iron, so a coordinated pose has nothing keeping it
        coordinated and the wall alone pushes the donor off the metal. The
        window is centred on the distance the pose arrives with, so 2.1 A
        stays 2.1 A; a fixed window instead pulls borderline poses in and
        drags the whole ligand toward the heme.
        """
        fe = next((np.array(p) for (s, _n), p in zip(heme_atoms, heme_xyz)
                   if s == "FE"), None)
        donor, d0, _el = _nearest_donor_to_fe(topology, positions,
                                              heme_atoms, heme_xyz)
        if fe is None or donor is None or d0 > cfg.fe_restrain_max_d0:
            return None, d0
        rmin = max(0.15, (d0 - cfg.fe_window_half) / 10.0)
        rmax = (d0 + cfg.fe_window_half) / 10.0
        k = cfg.fe_k_kcal * unit.kilocalories_per_mole / unit.angstrom**2
        f = mm.CustomExternalForce(
            "kfe*(max(0, rfe - rmax)^2 + max(0, rmin - rfe)^2);"
            "rfe = sqrt((x-fx)^2 + (y-fy)^2 + (z-fz)^2)")
        f.addGlobalParameter("kfe", k)
        f.addGlobalParameter("rmin", rmin)
        f.addGlobalParameter("rmax", rmax)
        for p in ("fx", "fy", "fz"):
            f.addPerParticleParameter(p)
        f.addParticle(donor, [float(c) for c in fe / 10.0])
        system.addForce(f)
        return donor, d0

    def minimize_complex(pdb_path, smiles, out_path, arm=None, verbose=False):
        pdb_path, out_path = Path(pdb_path), Path(out_path)
        arm = arm or cfg.arm
        if arm not in ARMS:
            raise ValueError(f"arm must be one of {ARMS}, got {arm!r}; "
                             f"in a mixed run pass the compound's arm")
        info = {"arm": arm, "fe_donor_raw": np.nan, "fe_donor_elem": None,
                "fe_restrained": False, "n_heme_centres": 0,
                "rms_force": np.nan, "max_force": np.nan, "converged": None,
                "lig_stereo": None, "stereo_ok": True,
                "stereo_mismatch": "", "stereo_undefined": ""}

        heme_atoms, heme_ang = read_heme(pdb_path, cfg)
        fe_pos = next((np.array(p) for (s, _n), p in zip(heme_atoms, heme_ang)
                       if s == "FE"), None)

        cfg.work_dir.mkdir(parents=True, exist_ok=True)
        prot_pdb = cfg.work_dir / f"{pdb_path.stem}_prot.pdb"
        lig_pdb = cfg.work_dir / f"{pdb_path.stem}_lig.pdb"
        _write_subset(pdb_path, prot_pdb, {cfg.prot_chain})
        prot_pdb = _fix_protein(
            prot_pdb, cfg.work_dir / f"{pdb_path.stem}_prot_fixed.pdb")
        _write_subset(pdb_path, lig_pdb, {cfg.lig_chain})

        # --- ligand: bond orders from the reference SMILES -----------------
        ref = Chem.MolFromSmiles(smiles)
        raw = Chem.MolFromPDBFile(str(lig_pdb), removeHs=True, sanitize=False)
        lig = AllChem.AssignBondOrdersFromTemplate(ref, raw)
        lig = Chem.AddHs(lig, addCoords=True)

        # Constitution must match. Stereochemistry is EXCLUDED: the template
        # may leave centres unspecified while the 3D pose defines them, which
        # is what failed all four cmpd_19 poses.
        if (Chem.MolToSmiles(Chem.RemoveHs(lig), isomericSmiles=False)
                != Chem.MolToSmiles(ref, isomericSmiles=False)):
            raise RuntimeError(
                f"{pdb_path.name}: bond-order transfer changed the molecule\n"
                f"  template:  {Chem.MolToSmiles(ref, isomericSmiles=False)}\n"
                f"  perceived: "
                f"{Chem.MolToSmiles(Chem.RemoveHs(lig), isomericSmiles=False)}")
        # Stereochemistry: RECORD it, do not kill the pose. Three outcomes,
        # and the old code conflated the last two:
        #   match      the pose has the configuration the template specifies
        #   mismatch   the generator produced the wrong enantiomer/diastereomer
        #   undefined  RDKit cannot determine it from THESE coordinates ("?"),
        #              which is not the same as wrong. A pseudo-asymmetric
        #              centre (lowercase r/s -- e.g. a 1,4-disubstituted
        #              cyclohexane, which is achiral and whose labels encode
        #              only cis/trans) goes undefined when the ring geometry
        #              is degenerate. RDKit returns the opposite letter, not
        #              "?", when the pose is genuinely the other diastereomer.
        # stage_select ranks stereo_ok poses first and falls back only when a
        # compound has none, so a bad pose is demoted rather than deleted and
        # there is always something to submit.
        ref_cc = dict(Chem.FindMolChiralCenters(
            ref, includeUnassigned=True, useLegacyImplementation=False))
        lig_cc = dict(Chem.FindMolChiralCenters(
            Chem.RemoveHs(lig), includeUnassigned=True,
            useLegacyImplementation=False))
        mism, undef = [], []
        for i, t in ref_cc.items():
            if t == "?":
                continue                       # template itself unspecified
            g = lig_cc.get(i)
            if g is None:
                continue                       # centre not perceived at all
            if g == "?":
                undef.append(f"{i}:{t}->?")
            elif g != t:
                mism.append(f"{i}:{t}->{g}")
        info["lig_stereo"] = ",".join(f"{i}{s}" for i, s in sorted(lig_cc.items()))
        info["stereo_mismatch"] = ";".join(mism)
        info["stereo_undefined"] = ";".join(undef)
        info["stereo_ok"] = not mism
        if mism or undef:
            _log(f"      [stereo] {pdb_path.name}: "
                 + (f"MISMATCH {mism}" if mism else "")
                 + ("  " if mism and undef else "")
                 + (f"undetermined {undef}" if undef else ""))
        if mism and cfg.stereo_strict:
            raise RuntimeError(
                f"{pdb_path.name}: stereocentre(s) {mism} (template->pose)")

        off_mol = Molecule.from_rdkit(lig, allow_undefined_stereo=True)
        off_mol.assign_partial_charges(cfg.nagl_model, toolkit_registry=NAGL)

        sysgen = SystemGenerator(
            forcefields=list(cfg.forcefields),
            small_molecule_forcefield=cfg.small_molecule_ff,
            molecules=[off_mol],
            cache=str(cfg.work_dir / "ligand_params.json"),
            forcefield_kwargs={"constraints": app.HBonds, "rigidWater": True},
            nonperiodic_forcefield_kwargs={"nonbondedMethod": app.NoCutoff},
        )

        prot = app.PDBFile(str(prot_pdb))
        modeller = app.Modeller(prot.topology, prot.positions)
        modeller.add(off_mol.to_topology().to_openmm(),
                     off_mol.conformers[0].to_openmm())
        modeller.addHydrogens(sysgen.forcefield,
                              platform=mm.Platform.getPlatformByName("CPU"))
        modeller.topology.setPeriodicBoxVectors(None)

        system = sysgen.create_system(modeller.topology, molecules=[off_mol])
        for f in system.getForces():
            if isinstance(f, mm.NonbondedForce):
                f.setNonbondedMethod(mm.NonbondedForce.NoCutoff)

        # --- positional restraints, flat-bottomed --------------------------
        # A plain harmonic on every ligand heavy atom pins the ligand's
        # INTERNAL geometry too, holding generator bond-angle distortions and
        # intramolecular clashes against GAFF -- which failed PoseBusters'
        # bond_angles and internal_steric_clash. A per-atom tolerance allows
        # local relaxation while still stopping drift. Protein tol is 0.
        force = mm.CustomExternalForce(
            "k0*max(0, sqrt((x-x0)^2 + (y-y0)^2 + (z-z0)^2) - tol0)^2")
        for p in ("k0", "tol0", "x0", "y0", "z0"):
            force.addPerParticleParameter(p)
        pos = modeller.positions
        bb = {"N", "CA", "C", "O"}
        for atom in modeller.topology.atoms():
            if atom.element == app.element.hydrogen:
                continue
            resn = atom.residue.name.strip().upper()
            if resn in cfg.lig_resnames:
                kk, tol = cfg.ligand_k, cfg.ligand_tol
            elif atom.name in bb:
                kk, tol = cfg.backbone_k, 0.0
            else:
                kk, tol = cfg.sidechain_k, 0.0
            x, y, z = pos[atom.index].value_in_unit(unit.nanometer)
            force.addParticle(atom.index, [
                kk * unit.kilocalories_per_mole / unit.angstrom**2,
                tol * unit.angstrom,
                x * unit.nanometer, y * unit.nanometer, z * unit.nanometer])
        system.addForce(force)

        # --- heme wall, then the Fe-donor restraint ------------------------
        lig_idx = heme_nm = rc_i = None
        if heme_atoms:
            info["n_heme_centres"] = len(heme_atoms)
            lig_idx, heme_nm, rc_i = _add_heme_wall(
                system, modeller.topology, heme_atoms, heme_ang, arm)
            di, d0, el = _nearest_donor_to_fe(
                modeller.topology, modeller.positions, heme_atoms, heme_ang)
            info["fe_donor_raw"], info["fe_donor_elem"] = d0, el
            if di is None:
                _log(f"      [!] no N/O/S donor in {pdb_path.name}")
            elif cfg.fe_restrain_for(arm):
                donor, _ = _add_fe_restraint(
                    system, modeller.topology, modeller.positions,
                    heme_atoms, heme_ang)
                info["fe_restrained"] = donor is not None
            if verbose:
                _log(f"      arm={arm} rc_fe={cfg.wall_rc_fe_for(arm)} "
                     f"rc_def={cfg.wall_rc_default_for(arm)} nm")
                _log(f"      wall {len(heme_atoms)}x{len(lig_idx)}  "
                     f"Fe-donor {el}{di} {d0:.2f} A  "
                     f"restrained={info['fe_restrained']}")
        else:
            _log(f"      [!] no heme in {pdb_path.name}")

        integrator = mm.LangevinMiddleIntegrator(
            300 * unit.kelvin, 1 / unit.picosecond, 0.002 * unit.picoseconds)
        sim = app.Simulation(modeller.topology, system, integrator)
        sim.context.setPositions(modeller.positions)

        if verbose and heme_atoms:
            xyz = np.asarray(modeller.positions.value_in_unit(unit.nanometer))
            e_mm = sim.context.getState(
                getEnergy=True, groups={WALL_GROUP}).getPotentialEnergy(
                ).value_in_unit(unit.kilojoules_per_mole)
            e_np = _wall_energy_numpy(xyz[lig_idx], heme_nm, rc_i)
            _log(f"      wall energy: OpenMM {e_mm:.3f} vs numpy {e_np:.3f}")
            if (not np.isfinite(e_mm)
                    or abs(e_mm - e_np) > max(1.0, 1e-4 * abs(e_np))):
                raise RuntimeError("wall energy mismatch")

        e0 = sim.context.getState(getEnergy=True).getPotentialEnergy()
        if not np.isfinite(e0.value_in_unit(unit.kilojoules_per_mole)):
            raise RuntimeError(f"starting energy is {e0}")
        sim.minimizeEnergy(
            tolerance=cfg.min_tol_kj * unit.kilojoules_per_mole / unit.nanometer,
            maxIterations=cfg.max_iter)
        state = sim.context.getState(getPositions=True, getEnergy=True,
                                     getForces=True)
        e1 = state.getPotentialEnergy()

        # minimizeEnergy() returns nothing and never says WHICH stopping
        # condition fired -- tolerance reached, or maxIterations hit. OpenMM's
        # tolerance is on the RMS over all force COMPONENTS, so recompute that
        # and compare. rms_force > min_tol_kj means it ran out of iterations
        # and the structure is NOT at a minimum.
        _f = np.asarray(state.getForces().value_in_unit(
            unit.kilojoules_per_mole / unit.nanometer))
        info["rms_force"] = float(np.sqrt((_f ** 2).sum() / _f.size))
        info["max_force"] = float(np.linalg.norm(_f, axis=1).max())
        info["converged"] = bool(info["rms_force"] <= cfg.min_tol_kj)
        if info["rms_force"] > cfg.warn_force_kj:
            _log(f"      [!] {pdb_path.name}: large residual force "
                 f"{info['rms_force']:.1f} kJ/mol/nm (max "
                 f"{info['max_force']:.0f}) after {cfg.max_iter} iterations "
                 f"-- well above this run's usual floor")

        full_path = out_path.with_name(out_path.stem + "_full.pdb")
        n_lig, n_heme = _write_full_complex(
            sim.topology, state.getPositions(), full_path,
            heme_atoms, heme_ang)
        strip_cys_hg_pdb(full_path, fe_pos, verbose=verbose)
        _strip_ligand_h(full_path, out_path)
        info["e_raw"] = e0.value_in_unit(unit.kilocalories_per_mole)
        info["e_min"] = e1.value_in_unit(unit.kilocalories_per_mole)
        return info

    return minimize_complex


def stage_refine(cfg: RunConfig, force=False, limit=None, only=None,
                 verbose_first=True) -> pd.DataFrame:
    """only: restrict to these compounds, in either notation
    (("cmpd_12",) or ("OCNT-2313038",)). None = every structure in pdb_dir."""
    cfg.ensure_dirs()
    smiles_of, _ = challenge_maps(cfg)
    minimize_complex = _refine_helpers(cfg)

    pdbs = sorted(cfg.pdb_dir.glob("*.pdb"))
    if only:
        want = {cfg.norm_cmpd(c, where="only") for c in
                ([only] if isinstance(only, str) else only)}
        n0 = len(pdbs)
        pdbs = [p for p in pdbs
                if (m := parse_pose(p, cfg)) and m["cmpd"] in want]
        _log(f"[refine] only={sorted(want)}: {len(pdbs)} of {n0} structures")
        if not pdbs:
            raise ValueError(f"no structures in {cfg.pdb_dir.name} match "
                             f"{sorted(want)}")
    if limit:
        pdbs = pdbs[:limit]
    _log(f"[refine] {len(pdbs)} structures -> {cfg.refine_dir}")

    # resolve every arm up front: an unassigned compound is a configuration
    # error, and finding it 40 minutes into a refinement is no use.
    if cfg.is_mixed:
        seen = {}
        for p in pdbs:
            m = parse_pose(p, cfg)
            if m:
                seen[m["cmpd"]] = seen.get(m["cmpd"], 0) + 1
        unassigned = [c for c in seen if c not in cfg.arm_map()]
        if unassigned and cfg.strict_arms:
            raise KeyError(
                f"[refine] unassigned compounds in {cfg.pdb_dir.name}: "
                f"{sorted(unassigned)}. Add them to cfg.binders / "
                f"cfg.nonbinders, or set strict_arms=False.")
        if unassigned:
            _log(f"  [!] {len(unassigned)} compounds -> default_arm="
                 f"{cfg.default_arm}: {sorted(unassigned)}")
        split = {}
        for c, n in seen.items():
            split[cfg.arm_of(c)] = split.get(cfg.arm_of(c), 0) + n
        _log("  arms: " + "  ".join(f"{a}={split.get(a, 0)} poses"
                                    for a in ARMS))

    rows, t0 = [], time.time()
    first_of_arm = set()
    for i, p in enumerate(pdbs, 1):
        info = parse_pose(p, cfg)
        if info is None:
            _log(f"  [!] {p.stem}: filename does not parse"); continue
        smi = smiles_of.get(info["cmpd"])
        if smi is None:
            _log(f"  [!] {p.stem}: no SMILES for {info['cmpd']}"); continue
        arm = cfg.arm_of(info["cmpd"])
        out = cfg.refine_dir / p.name
        full = out.with_name(out.stem + "_full.pdb")
        pen0 = heme_penetration(p, cfg)
        rec = {"cmpd": info["cmpd"], "seed": int(info["seed"]),
               "sample": int(info["sample"]), "structure": p.stem,
               "arm": arm,
               "clash_raw": min_lig_contact(p, cfg, heavy_only=True),
               "clash_all_raw": min_lig_contact(p, cfg),
               "heme_dist_raw": min_lig_heme(p, cfg),
               "through_ring_raw": pen0["through_ring"],
               "fe_donor_raw_in": pen0["min_fe_dist"]}
        try:
            if out.is_file() and out.stat().st_size > 0 and not force:
                rec["skipped"] = True
            else:
                # verbose on the first structure of each arm, so a mixed run
                # prints the wall radii actually used for both
                vb = verbose_first and arm not in first_of_arm
                first_of_arm.add(arm)
                rec.update(minimize_complex(p, smi, out, arm=arm, verbose=vb))
                rec["skipped"] = False
            pen = heme_penetration(out, cfg)
            rec.update(clash_min=min_lig_contact(out, cfg, heavy_only=True),
                       clash_all_min=min_lig_contact(out, cfg),
                       heme_dist_min=min_lig_heme(out, cfg),
                       through_ring=pen["through_ring"],
                       min_fe_dist=pen["min_fe_dist"],
                       lig_rmsd=lig_rmsd(p, out, cfg),
                       submission=str(out), full=str(full))
            el = time.time() - t0
            _log(f"  [{i}/{len(pdbs)}] {p.stem} [{arm[:2]}]: clash "
                 f"{rec['clash_raw']:.2f}->{rec['clash_min']:.2f}  heme "
                 f"{rec['heme_dist_raw']:.2f}->{rec['heme_dist_min']:.2f}  "
                 f"ring {rec['through_ring_raw']}->{rec['through_ring']}  "
                 f"RMSD {rec['lig_rmsd']:.2f}"
                 f"   (~{el/i*(len(pdbs)-i)/60:.0f} min left)")
        except Exception as ex:
            rec["error"] = str(ex)
            _log(f"  [{i}/{len(pdbs)}] [!] {p.stem}: FAILED {ex}")
        rows.append(rec)

    df = pd.DataFrame(rows)
    df.to_csv(cfg.manifest("refine"), index=False)
    _summarise_refine(df)
    return df


def _summarise_refine(df, label="refine"):
    if df.empty:
        return
    # a mixed run gets one block per arm: "heme <3.0 A" means coordination for
    # the binders and a failure for the non-binders, so a pooled count of it
    # cannot be read either way
    if "arm" in df and df["arm"].nunique() > 1:
        for a in ARMS:
            sub = df[df["arm"] == a]
            if len(sub):
                _summarise_refine(sub, label=f"refine:{a}")
        _log(f"[refine] {len(df)} structures total across "
             f"{df['arm'].nunique()} arms")
        return
    ok = df[df.get("error", pd.Series(index=df.index, dtype=object)).isna()]
    _log(f"[{label}] {len(ok)}/{len(df)} refined")
    if "lig_rmsd" in ok and ok["lig_rmsd"].notna().any():
        _log(f"  ligand RMSD : median {ok['lig_rmsd'].median():.2f}  "
             f"max {ok['lig_rmsd'].max():.2f} A   (want <~0.3)")
        _log(f"  clash heavy <2.2 A: {int((ok['clash_min'] < 2.2).sum())} "
             f"(was {int((ok['clash_raw'] < 2.2).sum())})")
        _log(f"  clash all   <1.9 A: {int((ok['clash_all_min'] < 1.9).sum())} "
             f"(was {int((ok['clash_all_raw'] < 1.9).sum())})")
        _log(f"  heme <3.0 A : {int((ok['heme_dist_min'] < 3.0).sum())} "
             f"(was {int((ok['heme_dist_raw'] < 3.0).sum())})")
        _log(f"  through ring: {int((ok['through_ring'] > 0).sum())} "
             f"(was {int((ok['through_ring_raw'] > 0).sum())})   <-- must be 0")
    if "rms_force" in ok and ok["rms_force"].notna().any():
        f = ok["rms_force"].dropna()
        med = f.median()
        # MAD, not sd: outliers inflate sd and then hide inside their own
        # threshold. 1.4826*MAD is the sd-equivalent for a normal core.
        mad = float((f - med).abs().median()) * 1.4826
        _log(f"  rms force   : median {med:.1f}  robust sd {mad:.2f}  "
             f"range {f.min():.1f}-{f.max():.1f} kJ/mol/nm")
        if med > 0 and mad / med < 0.10:
            _log(f"                spread is {mad/med:.1%} of the median -- a "
                 f"systematic floor from the positional restraints, not pose "
                 f"strain; iteration count will not move it")
        thr = max(2.0 * med, med + 8 * (mad if mad > 0 else 1.0))
        hi = f[f > thr]
        if len(hi):
            _log(f"  [!] {len(hi)} structure(s) above {thr:.0f} -- far enough "
                 f"off the floor to be real strain:")
            for i in hi.sort_values(ascending=False).index[:5]:
                _log(f"        {ok.loc[i, 'structure']}  "
                     f"rms {ok.loc[i, 'rms_force']:.0f}  "
                     f"heme {ok.loc[i, 'heme_dist_min']:.2f}  "
                     f"rmsd {ok.loc[i, 'lig_rmsd']:.2f}")
    if "converged" in ok and ok["converged"].notna().any():
        cv = ok["converged"].fillna(False).astype(bool)
        _log(f"  reached tol : {int(cv.sum())}/{len(cv)}"
             + ("" if cv.any() else "   (see the note above)"))
    if "error" in df and df["error"].notna().any():
        _log(f"  failures: {int(df['error'].notna().sum())}")
        for r in df[df["error"].notna()].itertuples():
            _log(f"    {r.structure}: {r.error}")


# ---------------------------------------------------------------------------
# stage 3: validate  (PoseBusters -- env_main)
# ---------------------------------------------------------------------------

def stage_validate(cfg: RunConfig, paths=None, force=False) -> pd.DataFrame:
    """PoseBusters on the submission-format PDBs.

    Ligand as an SDF with bond orders from the reference SMILES, and the heme
    KEPT in the conditioning file -- without it
    minimum_distance_to_organic_cofactors has no cofactor to measure, which
    is the check the validator failed three compounds on.
    """
    import gemmi
    from posebusters import PoseBusters
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from rdkit import RDLogger
    RDLogger.DisableLog("rdApp.*")

    cfg.ensure_dirs()
    pb_work = cfg.stage_dir / "pb_work"; pb_work.mkdir(exist_ok=True)
    smiles_of, _ = challenge_maps(cfg)

    if paths is None:
        paths = sorted(p for p in cfg.refine_dir.glob("*.pdb")
                       if not p.stem.endswith("_full"))
        # keep poses already validated unless force; PoseBusters is cheap but
        # there is no reason to redo 200 of them to pick up 11
        pbf = cfg.manifest("validate")
        if pbf.is_file() and not force:
            PV = pd.read_csv(pbf)
            if {"cmpd", "seed", "sample"} <= set(PV.columns):
                seen = set(zip(PV.cmpd, PV.seed.astype(int),
                               PV["sample"].astype(int)))
                n0 = len(paths)
                paths = [p for p in paths
                         if (m := parse_pose(p, cfg)) and
                         (m["cmpd"], int(m["seed"]),
                          int(m["sample"])) not in seen]
                _log(f"[validate] {n0 - len(paths)} already validated, "
                     f"{len(paths)} new")
                if not paths:
                    _log("[validate] nothing new; manifest unchanged")
                    return PV
    paths = [Path(p) for p in paths]
    _log(f"[validate] {len(paths)} structures")

    pb = PoseBusters(config=cfg.pb_config)
    rows, bool_checks = [], set()
    for p in paths:
        info = parse_pose(p, cfg)
        if info is None:
            _log(f"  [!] {p.stem}: does not parse"); continue
        smi = smiles_of.get(info["cmpd"])
        if smi is None:
            _log(f"  [!] {p.stem}: no SMILES"); continue
        rec = {"cmpd": info["cmpd"], "seed": int(info["seed"]),
               "sample": int(info["sample"]), "structure": p.stem,
               "pb_error": None}
        try:
            cond = pb_work / f"{p.stem}_cond.pdb"
            st = gemmi.read_structure(str(p)); st.setup_entities()
            for model in st:
                drop = [c.name for c in model
                        if any(r.name.strip() == cfg.lig_resname for r in c)]
                for cname in drop:
                    model.remove_chain(cname)
            st.write_pdb(str(cond))

            cx = Chem.MolFromPDBFile(str(p), sanitize=False, removeHs=False)
            frags = Chem.SplitMolByPDBResidues(cx)
            raw = frags.get(cfg.lig_resname)
            if raw is None or raw.GetNumAtoms() == 0:
                raise ValueError(f"no {cfg.lig_resname} residue "
                                 f"(present: {sorted(frags)})")
            # template without explicit H: the submission ligand has none
            mol = AllChem.AssignBondOrdersFromTemplate(
                Chem.MolFromSmiles(smi), raw)
            Chem.SanitizeMol(mol)
            sdf = pb_work / f"{p.stem}_lig.sdf"
            w = Chem.SDWriter(str(sdf)); w.write(mol); w.close()

            res = pb.bust(mol_pred=mol, mol_true=None,
                          mol_cond=str(cond), full_report=True)
            ignore = {"mol_pred", "mol_true", "mol_cond", "mol_true_loaded"}
            n_pass = n_total = 0
            for col in [c for c in res.columns if c not in ignore]:
                v = res[col].iloc[0]
                key = f"pb_{col.replace('pb_', '')}"
                if pd.isna(v):                 # bool(pd.NA) raises
                    rec[key] = None
                elif isinstance(v, (bool, np.bool_)):
                    bool_checks.add(key)
                    rec[key] = bool(v)
                    n_pass += int(v); n_total += 1
                else:
                    rec[key] = v               # measured value
            rec.update(pb_n_pass=n_pass, pb_n_total=n_total,
                       pb_all_pass=(n_total > 0 and n_pass == n_total))
            rec["lig_heme_min_dist"] = min_lig_heme(p, cfg)
        except Exception as ex:
            rec["pb_error"] = str(ex)
            _log(f"  [!] {p.stem}: {ex}")
        rows.append(rec)

    df = pd.DataFrame(rows)
    crit = [f"pb_{c}" for c in cfg.pb_critical if f"pb_{c}" in df.columns]
    if crit:
        df["pb_critical_pass"] = df[crit].apply(
            lambda r: all(v is not False for v in r), axis=1)
        df.loc[df.pb_error.notna(), "pb_critical_pass"] = False
    else:
        df["pb_critical_pass"] = df.get("pb_all_pass", False)
    df.attrs["bool_checks"] = sorted(bool_checks)

    # MERGE, never overwrite. stage_select calls this with paths=[one_pose]
    # for on-demand validation, and a plain write would replace the whole
    # manifest with that single row -- silently discarding every other
    # pose's verdict and leaving the next select with nothing to gate on.
    mf = cfg.manifest("validate")
    if mf.is_file():
        try:
            OLD = pd.read_csv(mf)
            if {"cmpd", "seed", "sample"} <= set(OLD.columns):
                n_new = len(df)
                df = (pd.concat([OLD, df], ignore_index=True)
                        .drop_duplicates(["cmpd", "seed", "sample"],
                                         keep="last"))
                if len(df) > n_new:
                    _log(f"[validate] merged {n_new} new into "
                         f"{len(df)} total")
        except Exception as ex:
            _log(f"[validate] could not merge existing manifest ({ex}); "
                 f"writing the new rows only")
    df.to_csv(mf, index=False)

    _log(f"[validate] critical pass {int(df.pb_critical_pass.sum())}/{len(df)}"
         f"   ({len(crit)} of {len(cfg.pb_critical)} critical checks present)")
    bc = sorted(bool_checks)
    if bc:
        fails = df[bc].apply(lambda c: (c == False).sum())
        fails = fails[fails > 0].sort_values(ascending=False)
        _log("  failures by check:\n" + (fails.to_string() if len(fails)
                                         else "    none"))
    for r in df[~df.pb_critical_pass].itertuples():
        bad = [c.replace("pb_", "") for c in crit if getattr(r, c, None) is False]
        _log(f"    {r.structure}: {bad or r.pb_error}")
    return df


# ---------------------------------------------------------------------------
# stage 4: score  (OrbMol-v2 -- env_score)
# ---------------------------------------------------------------------------

RES_CHARGE = {"ASP": -1, "GLU": -1, "LYS": +1, "ARG": +1}   # neutral HIS
CAP_NH, CAP_CH = 1.01, 1.09


def _score_helpers(cfg):
    import gemmi
    from ase import Atoms
    from scipy.spatial.distance import cdist

    def _atom(res, name):
        if res is None:
            return None
        for a in res:
            if a.name == name:
                return a
        return None

    def _cap_position(anchor, ref, fallback_ref, bond):
        if anchor is None:
            return None
        a = np.array(anchor.pos.tolist())
        if ref is not None:
            v = np.array(ref.pos.tolist()) - a
        elif fallback_ref is not None:
            v = a - np.array(fallback_ref.pos.tolist())
        else:
            return None
        n = np.linalg.norm(v)
        return None if n < 1e-6 else (a + bond * v / n).tolist()

    def residues_near(model, lig_xyz, shell, exclude_chains=()):
        keep = set()
        for ch in model:
            if ch.name in exclude_chains:
                continue
            for res in ch:
                if res.name not in AA3:
                    continue
                xyz = np.array([a.pos.tolist() for a in res
                                if a.element.name != "H"])
                if xyz.size == 0:
                    continue
                d = np.linalg.norm(xyz[:, None, :] - lig_xyz[None, :, :],
                                   axis=-1)
                if d.min() <= shell:
                    keep.add((ch.name, res.seqid.num))
        return keep

    def build_pocket(chain_residues, keep_keys):
        """-> (Atoms, formal_charge, n_caps). Caps each cut backbone bond."""
        syms, pos, charge, ncap = [], [], 0, 0
        for (cname, seqid) in sorted(keep_keys):
            res = chain_residues[cname][seqid]
            charge += RES_CHARGE.get(res.name, 0)
            for a in res:
                syms.append(a.element.name); pos.append(a.pos.tolist())
            if (cname, seqid - 1) not in keep_keys:
                h = _cap_position(
                    _atom(res, "N"),
                    _atom(chain_residues[cname].get(seqid - 1), "C"),
                    _atom(res, "CA"), CAP_NH)
                if h:
                    syms.append("H"); pos.append(h); ncap += 1
            if (cname, seqid + 1) not in keep_keys:
                h = _cap_position(
                    _atom(res, "C"),
                    _atom(chain_residues[cname].get(seqid + 1), "N"),
                    _atom(res, "CA"), CAP_CH)
                if h:
                    syms.append("H"); pos.append(h); ncap += 1
        at = Atoms(symbols=syms, positions=np.array(pos))
        at.info.update(charge=charge, spin=1)
        return at, charge, ncap

    def prep_pose(path, lig_charge):
        """Refined *_full.pdb -> (model, chain_residues, ligand Atoms, heavy xyz)."""
        st = gemmi.read_structure(str(path)); st.setup_entities()
        model = st[0]
        syms, pos = [], []
        for ch in model:
            for r in ch:
                if r.het_flag != "H":
                    continue
                nm = r.name.strip().upper()
                if nm in cfg.heme_resnames or nm in ("HOH", "WAT"):
                    continue
                if nm not in cfg.lig_resnames:
                    continue
                for a in r:
                    syms.append(a.element.name); pos.append(a.pos.tolist())
        if not syms:
            return None, "no ligand residue"
        pos = np.array(pos)
        heavy = pos[[s != "H" for s in syms]]
        nH = sum(1 for s in syms if s == "H")
        if len(heavy) < 6:
            return None, f"only {len(heavy)} heavy atoms"
        if nH:
            dmin = cdist(pos[[s == "H" for s in syms]], heavy).min(axis=1)
            orphans = int((dmin > 1.3).sum())
            if orphans:
                return None, (f"{orphans}/{nH} hydrogens >1.3 A from any "
                              f"ligand heavy atom (max {dmin.max():.2f})")
        lig = Atoms(symbols=syms, positions=pos)
        lig.info.update(charge=lig_charge, spin=1)
        return (model, {ch.name: {r.seqid.num: r for r in ch} for ch in model},
                lig, heavy), None

    def combine(a, b):
        at = a + b
        at.info.update(charge=a.info["charge"] + b.info["charge"], spin=1)
        return at

    return residues_near, build_pocket, prep_pose, combine


def stage_score(cfg: RunConfig, force=False, only=None) -> pd.DataFrame:
    """dE = E(complex) - E(pocket) - E(ligand) at the pose geometry.

    The pocket residue set is the UNION over a compound's poses, so the atom
    set is identical within a compound and the energies are comparable.
    Interaction energy, not binding free energy: no desolvation, no entropy,
    and not comparable BETWEEN compounds.
    """
    os.environ.setdefault("TORCHDYNAMO_DISABLE", "1")
    from rdkit import Chem
    from rdkit import RDLogger
    RDLogger.DisableLog("rdApp.*")
    import torch._dynamo
    torch._dynamo.config.disable = True
    from orb_models.forcefield.pretrained import orbmol_v2

    cfg.ensure_dirs()
    out_csv = cfg.manifest("score")

    # Per-pose resume, not all-or-nothing. Previously this returned the whole
    # existing manifest untouched, so a pose refined AFTER the first scoring
    # pass was never scored -- and force=True was the only alternative, which
    # rescored every pose on the GPU. Now already-scored poses are kept and
    # only the new ones run.
    # A stamp of the parameters that produced this manifest. Without it, a
    # re-score with a different `shell` into the same stage_dir would find
    # every (cmpd, seed, sample) already present and skip the lot, silently
    # returning the OLD numbers -- worse than overwriting, because nothing
    # looks wrong.
    stamp_f = cfg.stage_dir / "score_params.json"
    stamp = {"shell": cfg.shell, "max_atoms": cfg.max_atoms,
             "lig_chain": cfg.lig_chain, "heme_chain": cfg.heme_chain}

    PREV = pd.DataFrame()
    done = set()
    if out_csv.is_file() and not force:
        if stamp_f.is_file():
            old_stamp = json.loads(stamp_f.read_text())
            if old_stamp != stamp:
                diff = {k: (old_stamp.get(k), stamp[k]) for k in stamp
                        if old_stamp.get(k) != stamp[k]}
                raise RuntimeError(
                    f"{out_csv} was produced with different scoring "
                    f"parameters {diff} (was -> now). Re-scoring into it "
                    f"would skip the existing rows and return the OLD "
                    f"numbers. Give this config a distinct `variant` (or "
                    f"pass force=True to discard and rescore).")
        else:
            _log(f"[score] {out_csv.name} has no parameter stamp; assuming it "
                 f"matches shell={cfg.shell}")
        PREV = pd.read_csv(out_csv)
        if {"cmpd", "seed", "sample"} <= set(PREV.columns):
            done = set(zip(PREV.cmpd, PREV.seed.astype(int),
                           PREV["sample"].astype(int)))
            _log(f"[score] {len(done)} poses already scored in {out_csv.name}")
        else:
            PREV = pd.DataFrame()

    smiles_of, ocnt_of = challenge_maps(cfg)
    lig_charge = {}
    for c, s in smiles_of.items():
        m = Chem.MolFromSmiles(s)
        lig_charge[c] = Chem.GetFormalCharge(m) if m else None

    residues_near, build_pocket, prep_pose, combine = _score_helpers(cfg)
    exclude = (cfg.lig_chain, cfg.heme_chain)

    want = None
    if only:
        want = {cfg.norm_cmpd(c, where="only") for c in
                ([only] if isinstance(only, str) else only)}

    rows, n_seen, n_skip = [], 0, 0
    for f in sorted(cfg.refine_dir.glob("*_full.pdb")):
        info = parse_pose(f, cfg)
        if not info or info["cmpd"] not in smiles_of:
            continue
        n_seen += 1
        key = (info["cmpd"], int(info["seed"]), int(info["sample"]))
        if want is not None and info["cmpd"] not in want:
            continue
        if key in done:
            n_skip += 1
            continue
        rows.append({"cmpd": info["cmpd"], "id": ocnt_of[info["cmpd"]],
                     "seed": int(info["seed"]),
                     "sample": int(info["sample"]), "path": str(f)})
    POSES = pd.DataFrame(rows)
    _log(f"[score] {n_seen} poses on disk, {n_skip} already scored, "
         f"{len(POSES)} to score"
         + (f"   only={sorted(want)}" if want else ""))
    if POSES.empty:
        if len(PREV):
            _log("[score] nothing new to score; manifest unchanged")
            return PREV
        raise SystemExit(f"[score] no parseable *_full.pdb in {cfg.refine_dir}")

    MODEL, ADAPTER = orbmol_v2(device=cfg.score_device)

    def energy(at):
        g = ADAPTER.from_ase_atoms(at, device=cfg.score_device)
        return float(MODEL.predict(g, fp64_energy=True)["energy"].item())

    out = []
    for cmpd, grp in POSES.groupby("cmpd"):
        _log(f"  === {cmpd} ({len(grp)} poses) ===")
        prepped = {}
        for r in grp.itertuples():
            p, why = prep_pose(r.path, lig_charge[cmpd])
            if p is None:
                _log(f"    prep failed s{r.seed}_p{r.sample}: {why}"); continue
            prepped[(r.seed, r.sample)] = (p, r.path)
        if not prepped:
            continue
        keep = set()
        for (pm, _cr, _lg, lx), _p in prepped.values():
            keep |= residues_near(pm, lx, cfg.shell, exclude_chains=exclude)
        _log(f"    union pocket: {len(keep)} residues")

        for (seed, sample), ((pm, chres, lig, _lx), path) in sorted(prepped.items()):
            if any(k[1] not in chres.get(k[0], {}) for k in keep):
                _log(f"    s{seed}_p{sample}: pocket residue absent, skipped")
                continue
            pocket, q, ncap = build_pocket(chres, keep)
            cplx = combine(pocket, lig)
            over = len(cplx) > cfg.max_atoms
            e_c, e_p, e_l = energy(cplx), energy(pocket), energy(lig)
            out.append({"cmpd": cmpd, "seed": seed, "sample": sample,
                        "path": path, "n_pocket_res": len(keep),
                        "n_atoms_complex": len(cplx), "n_caps": ncap,
                        "over_max_atoms": over, "charge_pocket": q,
                        "charge_lig": lig.info["charge"],
                        "E_complex": e_c, "E_pocket": e_p, "E_ligand": e_l,
                        "dE": e_c - e_p - e_l})
            _log(f"    s{seed}_p{sample}: {len(cplx):4d} atoms  "
                 f"dE = {out[-1]['dE']:+9.3f} eV"
                 f"{'  [over max_atoms]' if over else ''}")

    df = pd.DataFrame(out)
    if len(PREV):
        df = (pd.concat([PREV, df], ignore_index=True)
                .drop_duplicates(["cmpd", "seed", "sample"], keep="last"))
        _log(f"[score] {len(PREV)} kept + {len(out)} new = {len(df)}")
    df.to_csv(out_csv, index=False)
    stamp_f.write_text(json.dumps(stamp, indent=1))
    _log(f"[score] {len(df)} poses -> {out_csv}")
    if len(df):
        _log(f"  dE range {df.dE.min():+.2f} to {df.dE.max():+.2f} eV")
        if df.dE.min() > 0:
            _log("  ALL dE POSITIVE -- the systems are clashing or "
                 "mis-assembled; check the refinement before trusting this")
    return df


# ---------------------------------------------------------------------------
# stage 5: select (+ plots)
# ---------------------------------------------------------------------------

def _load_iptm(cfg: RunConfig):
    path = Path(cfg.iptm_csv) if cfg.iptm_csv else cfg.run_dir / "harvest_iptm.csv"
    if not path.is_file():
        _log(f"[select] no ipTM file at {path}; ties will fall back to dE")
        return None
    ip = pd.read_csv(path)
    col = next((c for c in ("iptm_AC", "iptm", "pair_chains_iptm")
                if c in ip.columns), None)
    if col is None:
        _log(f"[select] no ipTM column in {path}; found {list(ip.columns)}")
        return None
    ip[col] = pd.to_numeric(ip[col], errors="coerce")

    # Two layouts in the wild:
    #   ESMFold2  cmpd,seed,sample,iptm_AC,...      -- already keyed
    #   OpenFold3 Target_ID,Model,seed,iptm_AC,...  -- `seed` is a STRING
    #             ("seed_3163119785"), there is no `sample`, and the Model
    #             stem carries all three. Derive the key from the stem.
    need = {"cmpd", "seed", "sample"}
    keyed = need <= set(ip.columns) and pd.to_numeric(
        ip["seed"], errors="coerce").notna().all()
    if not keyed:
        stem_col = None
        for c in ip.columns:
            vals = ip[c].astype(str).head(20)
            if vals.map(lambda v: parse_pose(v, cfg) is not None).all():
                stem_col = c; break
        if stem_col is None:
            _log(f"[select] {path.name}: no cmpd/seed/sample and no parseable "
                 f"model-stem column; found {list(ip.columns)}")
            return None
        k = ip[stem_col].astype(str).map(lambda v: parse_pose(v, cfg))
        ip = ip.assign(cmpd=[d["cmpd"] for d in k],
                       seed=[int(d["seed"]) for d in k],
                       sample=[int(d["sample"]) for d in k])
        _log(f"[select] {path.name}: keyed from the {stem_col!r} column "
             f"(no usable cmpd/seed/sample)")

    ip["seed"] = pd.to_numeric(ip["seed"], errors="coerce").astype("Int64")
    ip["sample"] = pd.to_numeric(ip["sample"], errors="coerce").astype("Int64")
    ip = ip.dropna(subset=["cmpd", "seed", "sample"])
    ip["seed"] = ip["seed"].astype(int)
    ip["sample"] = ip["sample"].astype(int)
    _log(f"[select] ipTM from {path.name} column {col!r} "
         f"({len(ip)} rows, {ip.cmpd.nunique()} compounds)")
    return ip[["cmpd", "seed", "sample", col]].rename(columns={col: "iptm"})


def stage_select(cfg: RunConfig, validate_candidates=True) -> pd.DataFrame:
    """Rank by dE, break ties within de_tol by ipTM, then walk down the
    ranking until a candidate passes the critical PoseBusters checks."""
    cfg.ensure_dirs()
    sc = pd.read_csv(cfg.manifest("score"))
    sc = sc.drop_duplicates(["cmpd", "seed", "sample"])
    # carry the arm each pose was refined with through to the selection table
    rf = cfg.manifest("refine")
    if rf.is_file():
        R = pd.read_csv(rf)
        carry = [c for c in ("arm", "stereo_ok", "stereo_mismatch",
                             "stereo_undefined")
                 if c in R.columns and c not in sc.columns]
        if carry:
            sc = sc.merge(
                R[["cmpd", "seed", "sample"] + carry].drop_duplicates(
                    ["cmpd", "seed", "sample"]),
                on=["cmpd", "seed", "sample"], how="left")
    if "arm" not in sc.columns:
        sc["arm"] = cfg.arm_label
    if "stereo_ok" not in sc.columns:
        sc["stereo_ok"] = True
    # a pose with no recorded verdict is treated as acceptable, not excluded
    sc["stereo_ok"] = sc["stereo_ok"].fillna(True).astype(bool)
    ip = _load_iptm(cfg)
    M = sc.merge(ip, on=["cmpd", "seed", "sample"], how="left") if ip is not None \
        else sc.assign(iptm=np.nan)

    pbf = cfg.manifest("validate")
    PB = pd.read_csv(pbf) if pbf.is_file() else None
    if PB is not None:
        M = M.merge(PB[["cmpd", "seed", "sample", "pb_critical_pass"]],
                    on=["cmpd", "seed", "sample"], how="left")
    else:
        M["pb_critical_pass"] = np.nan
        if cfg.require_pb_pass:
            _log("[select] no validate manifest; PoseBusters gating disabled")

    rows = []
    for cmpd, g in M.groupby("cmpd"):
        g = g.copy()
        n_all = len(g)
        if cfg.drop_positive_de and (g.dE < 0).any():
            g = g[g.dE < 0]

        # rank: dE ascending, ipTM descending within de_tol of the running best
        order = _rank_candidates(g, cfg)
        pick, by, reason, n_tried = None, None, "", 0
        for cand, how in order:
            n_tried += 1
            ok = cand.get("pb_critical_pass")
            if (not cfg.require_pb_pass) or PB is None or ok is True:
                pick, by = cand, how
                reason = "first candidate" if n_tried == 1 else \
                         f"candidate {n_tried} (earlier ones failed PoseBusters)"
                break
            if pd.isna(ok) and validate_candidates:
                sub = cfg.refine_dir / (f"{cmpd}_seed_{int(cand['seed'])}"
                                        f"_sample_{int(cand['sample'])}_model.pdb")
                if sub.is_file():
                    v = stage_validate(cfg, paths=[sub])
                    if len(v) and bool(v.pb_critical_pass.iloc[0]):
                        pick, by = cand, how
                        reason = f"candidate {n_tried}, validated on demand"
                        break
                    continue
            if ok is True:
                pick, by = cand, how
                break
        if pick is None and len(order):
            pick, by = order[0]
            reason = "no candidate passed PoseBusters; took rank 1"

        if pick is None:
            _log(f"  [!] {cmpd}: no candidates"); continue
        de_best = g.loc[g.dE.idxmin()]
        rows.append({
            "cmpd": cmpd, "arm": pick.get("arm"),
            "stereo_ok": pick.get("stereo_ok"),
            "stereo_mismatch": pick.get("stereo_mismatch", ""),
            "stereo_undefined": pick.get("stereo_undefined", ""),
            "n_stereo_ok": int(g.stereo_ok.fillna(True).astype(bool).sum())
                           if "stereo_ok" in g.columns else None,
            "n_poses": n_all, "n_scored": len(g),
            "seed": int(pick["seed"]), "sample": int(pick["sample"]),
            "decided_by": by, "reason": reason, "n_tried": n_tried,
            "dE": pick["dE"], "iptm": pick.get("iptm"),
            "dE_min": g.dE.min(), "dE_spread": g.dE.max() - g.dE.min(),
            "pb_critical_pass": pick.get("pb_critical_pass"),
            "same_as_dE_only": int(pick["seed"]) == int(de_best.seed),
            "path": pick.get("path")})

    SEL = pd.DataFrame(rows)
    if len(SEL):
        SEL = SEL.sort_values("cmpd",
                              key=lambda s: s.str.split("_").str[1].astype(int))
    cfg.best_dir.mkdir(parents=True, exist_ok=True)
    _, ocnt_of = challenge_maps(cfg)
    SEL["file"] = ""
    for i, r in SEL.iterrows():
        src = cfg.refine_dir / (f"{r.cmpd}_seed_{int(r.seed)}_sample_"
                                f"{int(r['sample'])}_model.pdb")
        if not src.is_file():
            _log(f"  [!] missing {src.name}"); continue
        dst = cfg.best_dir / f"{ocnt_of.get(r.cmpd, r.cmpd)}.pdb"
        shutil.copy2(str(src), str(dst))
        SEL.at[i, "file"] = dst.name
    SEL.to_csv(cfg.manifest("select"), index=False)

    _log(f"[select] {(SEL.file != '').sum()}/{len(SEL)} -> {cfg.best_dir}")
    if len(SEL):
        _log(SEL[["cmpd", "arm", "seed", "decided_by", "n_tried", "dE",
                  "iptm", "pb_critical_pass"]].round(3).to_string(index=False))
        n_ip = int((SEL.decided_by == "iptm").sum())
        _log(f"  decided by dE {len(SEL)-n_ip}, by ipTM tie-break {n_ip}")
        _log(f"  differs from dE-only pick: "
             f"{int((~SEL.same_as_dE_only).sum())}/{len(SEL)}")
        if SEL["arm"].notna().any():
            _log("  by arm: " + "  ".join(
                f"{a}={int((SEL.arm == a).sum())}" for a in ARMS))
        if "stereo_ok" in SEL.columns:
            n_bad = int((~SEL.stereo_ok.fillna(True).astype(bool)).sum())
            _log(f"  stereo: {len(SEL) - n_bad}/{len(SEL)} picks match the "
                 f"specified configuration")
            for _i, r in SEL[~SEL.stereo_ok.fillna(True).astype(bool)].iterrows():
                _log(f"    [!] {r.cmpd}: no stereo-correct pose existed; "
                     f"took best of {r.n_scored} ({r.get('stereo_mismatch','')})")
    audit_best(cfg)
    return SEL


# ---------------------------------------------------------------------------
# completeness audit and cross-run merge
# ---------------------------------------------------------------------------

def audit_best(cfg: RunConfig, expect=None, quiet=False) -> dict:
    """Is there exactly one best pose per challenge compound?

    A compound can go missing anywhere -- the refinement raised, the score
    stage dropped it, the filename did not parse -- and a short directory
    listing does not say which. This names the missing ones and the stage
    that lost them.
    """
    _, ocnt_of = challenge_maps(cfg)
    want = {c: ocnt_of.get(c, c) for c in ocnt_of}
    if expect is None:
        expect = len(want)
    have = {p.stem for p in cfg.best_dir.glob("*.pdb")} \
        if cfg.best_dir.is_dir() else set()
    missing = {c: o for c, o in want.items() if o not in have}
    extra = sorted(have - set(want.values()))

    # where each missing compound was last seen
    trail = {}
    for stage in ("convert", "refine", "validate", "score", "select"):
        f = cfg.manifest(stage)
        if not f.is_file():
            continue
        try:
            d = pd.read_csv(f)
        except Exception:
            continue
        if "cmpd" not in d.columns:
            continue
        for c in missing:
            sub = d[d.cmpd == c]
            if len(sub):
                err = ""
                if "error" in sub.columns and sub["error"].notna().any():
                    err = str(sub["error"].dropna().iloc[0])[:110]
                trail[c] = (stage, len(sub), err)

    out = {"n_have": len(have & set(want.values())), "n_expect": expect,
           "missing": missing, "extra": extra, "trail": trail}
    if not quiet:
        _log(f"[audit] {out['n_have']}/{expect} compounds in {cfg.best_dir}")
        for c, o in sorted(missing.items(),
                           key=lambda kv: int(kv[0].split('_')[-1])
                           if kv[0].split('_')[-1].isdigit() else 0):
            stage, n, err = trail.get(c, ("none", 0, ""))
            arm = cfg.arm_of(c) if cfg.is_mixed else cfg.arm
            _log(f"  [!] missing {c} ({o})  arm={arm}  "
                 f"last seen in {stage!r} manifest ({n} rows)"
                 + (f"  error: {err}" if err else ""))
        for e in extra:
            _log(f"  [?] extra file not in the challenge csv: {e}.pdb")
    return out


def merge_best(dest, *cfgs, expect=None, strict=True):
    """Union the best_poses of several runs into one submission directory.

    For arms already refined as separate runs. A compound selected by more
    than one run is an error, not a preference: it means the arms were not
    disjoint, and silently taking one would hide that.
    """
    dest = Path(dest)
    dest.mkdir(parents=True, exist_ok=True)
    seen, rows = {}, []
    for cfg in cfgs:
        for p in sorted(cfg.best_dir.glob("*.pdb")):
            if p.stem in seen:
                msg = (f"{p.stem} selected by both {seen[p.stem]} and "
                       f"{cfg.refine_dir.name}: the arms overlap")
                if strict:
                    raise RuntimeError(msg)
                _log(f"  [!] {msg}; keeping the first")
                continue
            seen[p.stem] = cfg.refine_dir.name
            shutil.copy2(str(p), str(dest / p.name))
            rows.append({"file": p.name, "arm": cfg.arm_label,
                         "run": cfg.run_name, "refine_dir": cfg.refine_dir.name,
                         "tag": cfg.tag})
    M = pd.DataFrame(rows)
    M.to_csv(dest / "merge_manifest.csv", index=False)
    _log(f"[merge] {len(M)} files -> {dest}")
    if len(M):
        _log("  " + "  ".join(f"{a}={int((M.arm == a).sum())}"
                              for a in sorted(M.arm.unique())))
    if expect is not None:
        _, ocnt_of = challenge_maps(cfgs[0])
        missing = sorted(set(ocnt_of.values()) - set(M.file.str[:-4]))
        if missing:
            _log(f"  [!] missing {len(missing)}: {missing}")
        if strict and len(M) != expect:
            raise RuntimeError(f"merged {len(M)} files, expected {expect}")
    return M


def _rank_candidates(g, cfg):
    """Ordered [(row_dict, 'dE'|'iptm'), ...] best first.

    Stereochemistry is a HARD PARTITION above the energy ranking: every
    stereo-correct pose outranks every mismatched one, and within each block
    dE and ipTM decide as usual. The competition scores against the specified
    configuration, so a wrong enantiomer is worth little however well it
    scores -- but a compound whose poses are ALL wrong still yields its best
    pose rather than nothing, which is the whole point of partitioning rather
    than filtering.
    """
    if cfg.prefer_stereo_ok and "stereo_ok" in g.columns:
        ok = g[g.stereo_ok.fillna(True).astype(bool)]
        bad = g[~g.stereo_ok.fillna(True).astype(bool)]
        if len(ok) and len(bad):
            return (_rank_block(ok, cfg) + _rank_block(bad, cfg))
    return _rank_block(g, cfg)


def _rank_block(g, cfg):
    rem = g.sort_values("dE").to_dict("records")
    order = []
    while rem:
        lo = rem[0]["dE"]
        tie = [r for r in rem if r["dE"] <= lo + cfg.de_tol]
        rest = [r for r in rem if r["dE"] > lo + cfg.de_tol]
        if len(tie) > 1 and any(pd.notna(r.get("iptm")) for r in tie):
            by_de = list(tie)                      # already dE-ascending
            by_ip = sorted(tie, key=lambda r: (-(r.get("iptm")
                           if pd.notna(r.get("iptm")) else -np.inf)))
            # Label by whether ipTM actually MOVED a pose, not by whether a
            # tie group existed. A pose that ipTM left where dE had put it
            # was decided by dE; only one promoted past a lower-dE pose was
            # decided by ipTM.
            de_rank = {(r["cmpd"], r["seed"], r["sample"]): i
                       for i, r in enumerate(by_de)}
            for i, r in enumerate(by_ip):
                k = (r["cmpd"], r["seed"], r["sample"])
                order.append((r, "iptm" if i < de_rank[k] else "dE"))
        else:
            order += [(r, "dE") for r in tie]
        rem = rest
    return order


def stage_plot(cfg: RunConfig):
    """Small multiples of dE vs rank per compound, plus the rank1->rank2 gap."""
    # Do NOT call matplotlib.use() unconditionally: it is global and
    # permanent for the kernel, so a notebook that runs this stage once finds
    # every later plt.show() silently drawing nothing. Modern matplotlib
    # already selects Agg when there is no display, and savefig works under
    # any backend, so only fall back if importing pyplot actually fails.
    try:
        import matplotlib
        import matplotlib as mpl
        try:
            import matplotlib.pyplot as plt
        except Exception:
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
    except Exception as ex:
        _log(f"[plot] matplotlib unavailable ({ex}); skipping")
        return None

    DF = pd.read_csv(cfg.manifest("score"))
    if DF.empty:
        _log("[plot] nothing to plot"); return None
    ip = _load_iptm(cfg)
    if ip is not None:
        DF = DF.merge(ip, on=["cmpd", "seed", "sample"], how="left")
    pbf = cfg.manifest("validate")
    if pbf.is_file():
        PB = pd.read_csv(pbf)
        DF = DF.merge(PB[["cmpd", "seed", "sample", "pb_critical_pass"]],
                      on=["cmpd", "seed", "sample"], how="left")

    ACCENT, EMPH, FAIL = "#3D63C9", "#C2410C", "#9a9a95"
    mpl.rcParams.update({"font.size": 9, "axes.linewidth": 0.6,
                         "axes.edgecolor": "#9a9a95", "figure.dpi": 110})
    DF["rank"] = DF.groupby("cmpd")["dE"].rank(method="first").astype(int)
    DF["ddE"] = DF["dE"] - DF.groupby("cmpd")["dE"].transform("min")

    cmpds = sorted(DF.cmpd.unique(), key=lambda c: int(c.split("_")[1]))
    ncol = 3
    nrow = -(-len(cmpds) // ncol)
    fig, axes = plt.subplots(nrow, ncol, figsize=(3.0 * ncol, 2.3 * nrow),
                             sharex=True, sharey=True)
    axes = np.atleast_1d(axes).ravel()
    for ax, c in zip(axes, cmpds):
        g = DF[DF.cmpd == c].sort_values("rank")
        ax.plot(g["rank"], g.ddE, "-", lw=1.4, color=ACCENT, zorder=2)
        for r in g.itertuples():
            failed = getattr(r, "pb_critical_pass", True) is False
            ax.plot(r.rank, r.ddE, "o", ms=4.5,
                    color=FAIL if failed else ACCENT,
                    mec="white", mew=0.8, zorder=3)
        b = g.iloc[0]
        ax.plot(b["rank"], b.ddE, "o", ms=8, color=EMPH, mec="white",
                mew=1.2, zorder=4)
        ax.set_title(f"{c}   best {b.dE:+.2f} eV", fontsize=9, pad=4,
                     loc="left")
        ax.axhline(cfg.de_tol, lw=1.0, ls="--", color=EMPH, zorder=1)
        ax.grid(axis="y", lw=0.4, color="#e3e3df", zorder=0)
        ax.set_axisbelow(True)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    for ax in axes[len(cmpds):]:
        ax.axis("off")
    for ax in axes[-ncol:]:
        ax.set_xlabel("rank (1 = lowest dE)")
    for ax in axes[::ncol]:
        ax.set_ylabel("dE - best  (eV)")
    fig.suptitle(f"dE profiles  {cfg.run_name} {cfg.arm_label} tag={cfg.tag}   "
                 f"(grey = fails critical PoseBusters; dashed = de_tol)",
                 fontsize=9, y=0.999, x=0.01, ha="left")
    fig.tight_layout()
    p1 = cfg.stage_dir / f"dE_profiles_{cfg.arm_label}_{cfg.tag}.png"
    fig.savefig(p1, dpi=200, bbox_inches="tight")
    plt.close(fig)

    two = DF[DF["rank"] <= 2].pivot(index="cmpd", columns="rank", values="dE")
    if two.shape[1] >= 2:
        gaps = two.assign(gap=lambda d: d[2] - d[1]).sort_values("gap")
        fig2, ax = plt.subplots(figsize=(5.6, 0.42 * len(gaps) + 1.1))
        ax.barh(gaps.index, gaps.gap, height=0.62, color=ACCENT, zorder=2)
        for y, v in zip(gaps.index, gaps.gap):
            ax.text(v + gaps.gap.max() * 0.02, y, f"{v:.2f}", va="center",
                    fontsize=8, color="#44443f")
        ax.axvline(cfg.de_tol, lw=1.2, ls="--", color=EMPH, zorder=3)
        ax.set_xlabel("dE(rank 2) - dE(rank 1)   (eV)")
        ax.set_title("Separation of the top pose from the runner-up",
                     fontsize=10, loc="left", pad=8)
        ax.grid(axis="x", lw=0.4, color="#e3e3df", zorder=0)
        ax.set_axisbelow(True)
        for s in ("top", "right", "left"):
            ax.spines[s].set_visible(False)
        ax.tick_params(left=False)
        fig2.tight_layout()
        p2 = cfg.stage_dir / f"dE_gap_{cfg.arm_label}_{cfg.tag}.png"
        fig2.savefig(p2, dpi=200, bbox_inches="tight")
        plt.close(fig2)
        _log(f"[plot] {p1}\n       {p2}")
        return p1, p2
    _log(f"[plot] {p1}")
    return (p1,)


# ---------------------------------------------------------------------------
# driver
# ---------------------------------------------------------------------------

STAGES = ("convert", "refine", "validate", "score", "select", "plot")


def describe_runs(predictions_root, run_name=None, params=None):
    """Scan refine directories and say which tag used which parameters.

    Every stage writes config.json beside its outputs, so a directory name
    like refined_mixed_4612ed35 can always be resolved back to the settings
    that produced it -- weeks later, with no kernel state.
    """
    root = Path(predictions_root)
    roots = [root / run_name] if run_name else sorted(
        d for d in root.iterdir() if d.is_dir())
    params = params or ("shell", "wall_rc_default", "wall_rc_fe", "ligand_k",
                        "ligand_tol", "max_iter", "min_tol_kj", "arm")
    rows = []
    for r in roots:
        for cj in sorted(r.glob("refined_*/**/config.json")):
            try:
                d = json.loads(cj.read_text())
            except Exception:
                continue
            rec = {"run": r.name, "refine_dir": cj.parent.parent.name,
                   "stage_dir": cj.parent.name,
                   "tag": (d.get("_derived") or {}).get("tag", "?")}
            for k in params:
                v = d.get(k)
                rec[k] = json.dumps(v) if isinstance(v, (dict, list)) else v
            rows.append(rec)
    T = pd.DataFrame(rows).drop_duplicates(
        subset=["refine_dir", "stage_dir"]).reset_index(drop=True)
    if T.empty:
        _log(f"[runs] no config.json under {root}; they are written from now "
             f"on, so older directories may have none")
    return T


def run_stage(cfg: RunConfig, name: str, cross_env=True, **kw):
    """Run one stage, shelling out to env_score for `score` when needed."""
    # Provenance: a directory name must always be resolvable back to the
    # settings that produced it, without the kernel that made it.
    try:
        cfg.ensure_dirs()
        cfg.save(cfg.stage_dir / "config.json")
    except Exception as ex:
        _log(f"[{name}] could not write config.json ({ex})")

    if name == "score" and cross_env and not _can_import("orb_models"):
        cfg_path = cfg.stage_dir / "config.json"
        cfg.save(cfg_path)
        cmd = ["conda", "run", "--no-capture-output", "-n", cfg.env_score,
               "python", str(Path(__file__).resolve()),
               "--config", str(cfg_path), "--stage", "score"]
        _log(f"[score] -> {cfg.env_score}: {' '.join(cmd)}")
        env = dict(os.environ, PYTHONPATH=str(Path(__file__).resolve().parent))
        subprocess.run(cmd, check=True, env=env)
        return pd.read_csv(cfg.manifest("score"))
    return {"convert": stage_convert, "refine": stage_refine,
            "validate": stage_validate, "score": stage_score,
            "select": stage_select, "plot": stage_plot}[name](cfg, **kw)


def _can_import(mod):
    import importlib.util
    try:
        return importlib.util.find_spec(mod) is not None
    except Exception:
        return False


def run_all(cfg: RunConfig, stages=STAGES, cross_env=True, **kw):
    cfg.ensure_dirs()
    cfg.save(cfg.stage_dir / "config.json")
    _log(cfg.describe()); _log("")
    out = {}
    for s in stages:
        t0 = time.time()
        out[s] = run_stage(cfg, s, cross_env=cross_env,
                           **(kw.get(s) or {}))
        _log(f"[{s}] done in {time.time()-t0:.0f}s\n")
    return out


def _slug(v):
    """Short filesystem-safe label for a parameter value."""
    if isinstance(v, bool):
        return "T" if v else "F"
    if isinstance(v, float):
        return f"{v:g}".replace(".", "p").replace("-", "m")
    if isinstance(v, int):
        return str(v)
    if isinstance(v, dict):
        return "".join(f"{k[:2]}{_slug(x)}" for k, x in sorted(v.items()))
    return re.sub(r"[^A-Za-z0-9]+", "", str(v))[:10]


def sweep(base: RunConfig, grid: dict, prefix: str = ""):
    """Yield configs over a parameter grid, each writing somewhere distinct.

    Two kinds of key, and they need different handling:

      in REFINE_PARAMS (ligand_k, wall_rc_*, ...) -- changes cfg.tag, so each
        combination already gets its own refine_dir. Nothing more to do.

      NOT in REFINE_PARAMS (shell, de_tol, max_atoms, require_pb_pass, ...) --
        deliberately excluded so that re-scoring does not force a
        re-refinement. That means the tag does NOT change, so every
        combination would write over the same score.csv and best_poses. These
        are folded into `variant` instead, which redirects score/select/plot
        while the refined structures stay shared.
    """
    import itertools
    keys = list(grid)
    for combo in itertools.product(*(grid[k] for k in keys)):
        over = dict(zip(keys, combo))
        d = asdict(base)
        d.update(over)
        d = {k: (tuple(v) if isinstance(v, list) else v) for k, v in d.items()}
        soft = {k: v for k, v in over.items() if k not in REFINE_PARAMS}
        if soft:
            lab = "_".join(f"{k}{_slug(v)}" for k, v in sorted(soft.items()))
            d["variant"] = f"{prefix}{lab}" if not base.variant else \
                           f"{base.variant}_{prefix}{lab}"
        yield RunConfig(**d)


def sweep_plan(cfgs, strict=True):
    """Print where every config in a sweep will write, and refuse collisions.

    Run this BEFORE a sweep. Two configs collide when they differ in some
    parameter but resolve to the same output directory -- the second then
    overwrites the first, silently. That happens when a parameter outside
    REFINE_PARAMS is changed without a `variant`, since those deliberately do
    not affect cfg.tag. cs.sweep() assigns the variant for you; a hand-rolled
    dataclasses.replace does not.
    """
    rows, by_dir = [], {}
    for c in cfgs:
        key = str(c.stage_dir)
        sf = c.manifest("score")
        stamp_f = c.stage_dir / "score_params.json"
        conflict = ""
        if sf.is_file() and stamp_f.is_file():
            try:
                old_stamp = json.loads(stamp_f.read_text())
                if old_stamp.get("shell") != c.shell:
                    conflict = f"shell {old_stamp.get('shell')}!={c.shell}"
            except Exception:
                pass
        elif sf.is_file():
            conflict = "unstamped"
        rows.append({"tag": c.tag, "variant": c.variant or "-",
                     "refine_dir": c.refine_dir.name,
                     "stage_dir": c.stage_dir.name,
                     "refine_exists": c.refine_dir.is_dir(),
                     "score_exists": sf.is_file(),
                     "conflict": conflict or "-"})
        by_dir.setdefault(key, []).append(c)

    T = pd.DataFrame(rows)
    _log(T.to_string(index=False))

    clashes = []
    for key, group in by_dir.items():
        if len(group) < 2:
            continue
        seen = {}
        for c in group:
            sig = json.dumps({k: _jsonable(v) for k, v in asdict(c).items()},
                             sort_keys=True)
            seen.setdefault(sig, []).append(c)
        if len(seen) > 1:
            diff = set()
            sigs = [json.loads(k) for k in seen]
            for k in sigs[0]:
                if len({json.dumps(s.get(k), sort_keys=True) for s in sigs}) > 1:
                    diff.add(k)
            clashes.append((Path(key).name, sorted(diff)))

    soft = [f.name for f in fields(RunConfig)
            if f.name not in REFINE_PARAMS
            and f.name in ("shell", "max_atoms", "de_tol", "drop_positive_de",
                           "require_pb_pass", "prefer_stereo_ok", "pb_config",
                           "score_device", "stereo_strict")]
    if clashes:
        for d, diff in clashes:
            _log(f"  [!] COLLISION in {d}: configs differ in {diff} but write "
                 f"to the same place; the later one overwrites the earlier")
        _log(f"      Those parameters do not change cfg.tag. Pass a distinct "
             f"`variant`, or build the sweep with cs.sweep(), which does it.")
        if strict:
            raise RuntimeError(f"{len(clashes)} output collision(s); nothing run")
    else:
        _log(f"  OK: {len(cfgs)} configs -> {len(by_dir)} distinct output "
             f"directories, no collisions")
    stale = T[T.conflict.str.startswith("shell")]
    if len(stale):
        for _i, r in stale.iterrows():
            _log(f"  [!] {r.stage_dir} already holds results from a different "
                 f"{r.conflict}; stage_score will refuse. Give this config a "
                 f"distinct `variant`.")
        if strict:
            raise RuntimeError(f"{len(stale)} config(s) target a stage_dir "
                               f"holding other parameters' results")
    reuse = T[T.refine_exists]
    if len(reuse):
        _log(f"  {len(reuse)} config(s) reuse an existing refine_dir "
             f"(no re-refinement): {list(reuse.tag)}")
    _log(f"  parameters that need a variant if varied by hand: {soft}")
    return T


def compare_runs(cfgs):
    """One row per config: refinement metrics plus PoseBusters pass rate."""
    rows = []
    for c in cfgs:
        # resolved per-arm radii, which are what actually act -- more use
        # than the raw wall_rc_default field, which may be None or a dict
        rec = {"tag": c.tag, "variant": c.variant or "-", "shell": c.shell,
               "rc_nonbind": c.wall_rc_default_for("nonbinders"),
               "rc_bind": c.wall_rc_default_for("binders"),
               "fe_nonbind": c.wall_rc_fe_for("nonbinders"),
               "fe_bind": c.wall_rc_fe_for("binders"),
               "ligand_k": c.ligand_k, "ligand_tol": c.ligand_tol,
               "wall_k": c.wall_k, "arm": c.arm_label}
        rf = c.manifest("refine")
        if rf.is_file():
            R = pd.read_csv(rf)
            rec.update(n=len(R),
                       rmsd_med=R.get("lig_rmsd", pd.Series(dtype=float)).median(),
                       rmsd_max=R.get("lig_rmsd", pd.Series(dtype=float)).max(),
                       clash_all_lt19=int((R.get("clash_all_min",
                                                 pd.Series(dtype=float)) < 1.9).sum()),
                       heme_lt30=int((R.get("heme_dist_min",
                                            pd.Series(dtype=float)) < 3.0).sum()),
                       # the wall's ACTUAL job -- must stay 0 at every radius
                       through_ring=int((R.get("through_ring",
                                               pd.Series(dtype=float)) > 0).sum()),
                       heme_min=R.get("heme_dist_min",
                                      pd.Series(dtype=float)).min(),
                       converged=(int(R["converged"].fillna(False).astype(bool).sum())
                                  if "converged" in R.columns else None))
        vf = c.manifest("validate")
        if vf.is_file():
            V = pd.read_csv(vf)
            rec["pb_critical_pass"] = int(V.pb_critical_pass.sum())
            rec["pb_n"] = len(V)
        sf = c.manifest("score")
        if sf.is_file():
            S = pd.read_csv(sf)
            rec["dE_min"] = S.dE.min()
            rec["dE_med"] = S.dE.median()
            rec["n_scored"] = len(S)
            if "n_atoms_complex" in S.columns:
                rec["atoms_med"] = S.n_atoms_complex.median()
                rec["atoms_max"] = S.n_atoms_complex.max()
            if "over_max_atoms" in S.columns:
                rec["over_max"] = int(S.over_max_atoms.astype(bool).sum())
        lf = c.manifest("select")
        if lf.is_file():
            L = pd.read_csv(lf)
            rec["n_picked"] = len(L)
            if "pb_critical_pass" in L.columns:
                rec["picks_pb_pass"] = int((L.pb_critical_pass == True).sum())
        rows.append(rec)
    return pd.DataFrame(rows)


def preview_de_tol(cfg, values=(0.02, 0.05, 0.10, 0.20, 0.40)):
    """How would the picks change at other de_tol? READ-ONLY, no re-scoring.

    de_tol only affects ranking, so sweeping it through cs.sweep would force
    a needless re-score (score is variant-owned). This re-ranks the existing
    score manifest in memory and writes nothing.

    Reports, per de_tol: how many compounds are decided by dE vs ipTM, and
    how many picks differ from the current setting.
    """
    sc = pd.read_csv(cfg.manifest("score")).drop_duplicates(
        ["cmpd", "seed", "sample"])
    rf = cfg.manifest("refine")
    if rf.is_file():
        R = pd.read_csv(rf)
        carry = [c for c in ("stereo_ok",) if c in R.columns and c not in sc]
        if carry:
            sc = sc.merge(R[["cmpd", "seed", "sample"] + carry].drop_duplicates(
                ["cmpd", "seed", "sample"]), on=["cmpd", "seed", "sample"],
                how="left")
    if "stereo_ok" not in sc:
        sc["stereo_ok"] = True
    sc["stereo_ok"] = sc["stereo_ok"].fillna(True).astype(bool)

    ip = _load_iptm(cfg)
    M = sc.merge(ip, on=["cmpd", "seed", "sample"], how="left") \
        if ip is not None else sc.assign(iptm=np.nan)
    pbf = cfg.manifest("validate")
    if pbf.is_file():
        PB = pd.read_csv(pbf)
        if "pb_critical_pass" in PB.columns:
            M = M.merge(PB[["cmpd", "seed", "sample", "pb_critical_pass"]]
                        .drop_duplicates(["cmpd", "seed", "sample"]),
                        on=["cmpd", "seed", "sample"], how="left")
    if "pb_critical_pass" not in M:
        M["pb_critical_pass"] = np.nan

    import dataclasses as _dc
    out, picks = [], {}
    for t in values:
        c2 = _dc.replace(cfg, de_tol=float(t))
        by, chosen = {"dE": 0, "iptm": 0}, {}
        for cmpd, g in M.groupby("cmpd"):
            if c2.drop_positive_de and (g.dE < 0).any():
                g = g[g.dE < 0]
            ranked = _rank_candidates(g, c2)
            pick = None
            for cand, how in ranked:
                ok = cand.get("pb_critical_pass")
                if (not c2.require_pb_pass) or ok is True or pd.isna(ok):
                    pick = (cand, how); break
            if pick is None and ranked:
                pick = ranked[0]
            if pick is None:
                continue
            cand, how = pick
            by[how] = by.get(how, 0) + 1
            chosen[cmpd] = (int(cand["seed"]), int(cand["sample"]), how)
        picks[t] = chosen
        out.append({"de_tol": t, "by_dE": by.get("dE", 0),
                    "by_iptm": by.get("iptm", 0), "n": len(chosen)})

    cur = cfg.de_tol
    ref = picks.get(cur) or picks[values[0]]
    for row in out:
        pk = picks[row["de_tol"]]
        row["differs_from_current"] = sum(
            1 for c in ref if c in pk and pk[c][:2] != ref[c][:2])
    T = pd.DataFrame(out)
    _log(f"[de_tol] current = {cur}; "
         f"{len(ref)} compounds, read-only preview")
    _log(T.to_string(index=False))
    _log("  by_dE rises as de_tol falls: ties shrink and the energy decides.")
    _log("  If by_iptm stays high even at 0.02, the dE values really are "
         "degenerate and ipTM is the only thing separating the poses.")
    return T, picks


def compare_best_structures(cfgs, labels=None):
    """Compare the actual submitted COORDINATES across a sweep.

    For a SCORING sweep (shell, de_tol), identical picks mean an identical
    submission -- the structures were never re-refined. For a REFINEMENT
    sweep (wall_rc_*, ligand_k) that is false: every variant re-minimised,
    so the same pick is a DIFFERENT set of coordinates. compare_selections
    would report "no change" while the files on disk differ.

    -> one row per compound: whether the pick moved, how far the ligand moved
    from the baseline's submitted pose, and its closest heme approach.
    """
    labels = labels or [f"{c.variant or ''}:{c.tag}".lstrip(":") for c in cfgs]
    _, ocnt_of = challenge_maps(cfgs[0])

    picks = {}
    for c, lab in zip(cfgs, labels):
        f = c.manifest("select")
        picks[lab] = ({r["cmpd"]: (int(r["seed"]), int(r["sample"]))
                       for r in pd.read_csv(f).to_dict("records")}
                      if f.is_file() else {})

    base_lab, base_cfg = labels[0], cfgs[0]
    rows = []
    for cmpd, ocnt in sorted(ocnt_of.items(),
                             key=lambda kv: int(kv[0].split("_")[-1])
                             if kv[0].split("_")[-1].isdigit() else 0):
        ref = base_cfg.best_dir / f"{ocnt}.pdb"
        if not ref.is_file():
            continue
        rec = {"cmpd": cmpd, "ocnt": ocnt}
        for c, lab in zip(cfgs, labels):
            f = c.best_dir / f"{ocnt}.pdb"
            if not f.is_file():
                rec[f"rmsd_{lab}"] = np.nan; continue
            rec[f"pick_{lab}"] = "s{}p{}".format(*picks[lab][cmpd]) \
                if cmpd in picks[lab] else "?"
            rec[f"rmsd_{lab}"] = (0.0 if lab == base_lab
                                  else lig_rmsd(ref, f, c))
            rec[f"heme_{lab}"] = min_lig_heme(f, c)
        pk = {v for k, v in rec.items() if k.startswith("pick_")}
        rec["pick_changed"] = len(pk) > 1
        rows.append(rec)

    T = pd.DataFrame(rows)
    if T.empty:
        _log("[compare] no best_poses found"); return T

    rcols = [c for c in T.columns if c.startswith("rmsd_")]
    hcols = [c for c in T.columns if c.startswith("heme_")]
    _log(f"[compare] {len(T)} compounds, baseline = {base_lab}")
    _log(f"  picks changed      : {int(T.pick_changed.sum())}/{len(T)}")
    moved = (T[rcols].drop(columns=[f"rmsd_{base_lab}"], errors="ignore")
             > 0.25).any(axis=1)
    _log(f"  coordinates moved  : {int(moved.sum())}/{len(T)} "
         f"(ligand RMSD > 0.25 A from the baseline submission)")
    for cname in rcols:
        if cname == f"rmsd_{base_lab}":
            continue
        v = T[cname].dropna()
        if len(v):
            _log(f"    {cname[5:]:14s} median {v.median():.3f}  "
                 f"max {v.max():.3f} A")
    _log("  closest heme approach of the SUBMITTED poses:")
    for cname in hcols:
        v = T[cname].dropna()
        if len(v):
            _log(f"    {cname[5:]:14s} median {v.median():.2f}  "
                 f"min {v.min():.2f} A   "
                 f"({int((v < 3.0).sum())}/{len(v)} under 3.0)")
    return T


def compare_selections(cfgs, labels=None):
    """Which compound picks change across a sweep? -> (table, summary).

    For a scoring sweep the absolute energies do not matter; only whether the
    ranking puts a different pose first. This lines the picks up per compound
    so a shell that changes nothing is obvious at a glance.
    """
    labels = labels or [c.variant or c.tag for c in cfgs]
    cols = {}
    for c, lab in zip(cfgs, labels):
        f = c.manifest("select")
        if not f.is_file():
            _log(f"[compare] no select manifest for {lab}"); continue
        L = pd.read_csv(f)
        cols[lab] = {r["cmpd"]: f"s{int(r['seed'])}p{int(r['sample'])}"
                     for r in L.to_dict("records")}
    if not cols:
        return pd.DataFrame(), pd.DataFrame()
    T = pd.DataFrame(cols)
    T = T.loc[sorted(T.index, key=lambda c: int(str(c).split("_")[-1])
                     if str(c).split("_")[-1].isdigit() else 0)]
    first = T.columns[0]
    T["same_as_" + str(first)] = T.apply(
        lambda r: all(r[c] == r[first] for c in T.columns if c != first), axis=1)
    summary = pd.DataFrame([
        {"variant": c, "n": int(T[c].notna().sum()),
         "differs_from_" + str(first): int((T[c] != T[first]).sum())}
        for c in cols])
    return T, summary


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("--config", required=True)
    ap.add_argument("--stage", nargs="+", default=list(STAGES),
                    choices=list(STAGES))
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--limit", type=int, default=None)
    ap.add_argument("--no-cross-env", action="store_true")
    a = ap.parse_args(argv)

    cfg = RunConfig.load(a.config)
    kw = {}
    if a.force:
        for s in ("convert", "refine", "validate", "score"):
            kw[s] = {"force": True}
    if a.limit:
        kw.setdefault("refine", {})["limit"] = a.limit
    run_all(cfg, stages=tuple(a.stage), cross_env=not a.no_cross_env, **kw)
    return 0


if __name__ == "__main__":
    sys.exit(main())
