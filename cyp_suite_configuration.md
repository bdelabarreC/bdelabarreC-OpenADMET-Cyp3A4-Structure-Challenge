# `cyp_suite` configuration reference

Every stage in `cyp_suite` takes a single `RunConfig`. There are no other
arguments that change what gets computed — the stage functions accept only
execution controls (resume, restrict, limit). So this document is mostly a
reference for `RunConfig` fields, followed by the entry points.

```python
import cyp_suite as cs

cfg = cs.RunConfig.from_env(
    run_name   = "esmfold2_cyp3a4_run9",
    engine     = "esmfold2",
    arm        = "mixed",
    binders    = ("cmpd_0", "cmpd_1", ...),
    nonbinders = ("cmpd_2", "cmpd_3", ...),
    shell      = 3.5,
    variant    = "shell3p5",
)
results = cs.run_all(cfg)
```

---

## 1. Three mechanisms worth understanding first

Most confusion with this configuration comes from these three, so they are
worth reading before the field tables.

### `cfg.tag` — what identifies a refinement

`tag` is a SHA1 (first 8 hex characters) of the fields listed in
`REFINE_PARAMS`, plus the resolved compound→arm map when the run is mixed:

```python
REFINE_PARAMS = (
    "backbone_k", "sidechain_k", "ligand_k", "ligand_tol",
    "wall_k", "wall_rc_default", "wall_rc_fe", "wall_rc_pyrrole",
    "fe_restrain", "fe_restrain_max_d0", "fe_window_half", "fe_k_kcal",
    "max_iter", "min_tol_kj", "forcefields", "small_molecule_ff",
)
```

`refine_dir` is `run_dir / f"refined_{arm_label}_{tag}"`, so changing any
refinement parameter automatically produces a separate output tree and cannot
overwrite earlier work. Changing the arm assignment does the same.

A consequence worth knowing: `wall_rc_default = 0.26` and
`wall_rc_default = {"nonbinders": 0.26}` hash to **different tags** even when
they produce identical radii for the run in question, because the raw field
value is what is hashed.

### `variant` — isolating a scoring change

`shell`, `de_tol`, `max_atoms` and `require_pb_pass` are deliberately **not**
in `REFINE_PARAMS`: re-scoring must not force a re-refinement. That means they
do not change `tag`, and without further action two different `shell` values
would write over the same `score.csv` and `best_poses`.

`variant` solves that. With it set:

| stage | output location |
|---|---|
| `convert`, `refine`, `validate` | `refine_dir/_suite/` (shared) |
| `score`, `select`, `plot` | `refine_dir/_suite_<variant>/` |
| best poses | `refine_dir/best_poses_<variant>/` |

(`VARIANT_STAGES = ("score", "select", "plot")`.)

Sharing `validate` across variants is intentional — a PoseBusters verdict
depends on the structure, not on the scoring shell.

`stage_score` also writes a `score_params.json` stamp recording `shell`,
`max_atoms`, `lig_chain` and `heme_chain`. Re-scoring into a manifest produced
with different parameters **raises** rather than silently returning the old
numbers. Give the new settings their own `variant`, or pass `force=True`.

### Per-arm overrides — float versus dict

Four refinement fields accept a float, a dict, or `None`:

```python
wall_rc_default = 0.26                    # pins BOTH arms to 0.26
wall_rc_default = {"nonbinders": 0.26}    # pins non-binders; binders keep default
wall_rc_default = None                    # both arms keep their defaults
```

The dict form is what a clean single-arm ablation needs. The resolution rule
is in `RunConfig._override`.

---

## 2. Identity and arms

| field | type | default | meaning |
|---|---|---|---|
| `run_name` | `str` | *required* | directory under `predictions_root` |
| `engine` | `str` | `"esmfold2"` | `"esmfold2"` or `"openfold3"`; affects pose-name parsing |
| `arm` | `str` | `"nonbinders"` | `"binders"`, `"nonbinders"` or `"mixed"` |

Setting `binders`, `nonbinders` or `class_csv` makes the run mixed: one
`cif_dir`, one `refine_dir`, one `best_dir`, but each structure refined with
its own compound's arm parameters.

| field | type | default | meaning |
|---|---|---|---|
| `binders` | `tuple` | `()` | compounds in the binder arm |
| `nonbinders` | `tuple` | `()` | compounds in the non-binder arm |
| `class_csv` | `str` | `""` | csv with a `cmpd` column plus `class_col` |
| `class_col` | `str` | `"arm"` | column holding `binders`/`nonbinders` |
| `default_arm` | `str` | `"nonbinders"` | arm for compounds in neither list |
| `strict_arms` | `bool` | `True` | raise if a compound is unassigned |

Explicit lists win over `class_csv`. Compound naming is forgiving —
`cmpd_12`, `cmpd12`, `CMPD-12`, `12` and `OCNT-2313038` all resolve to the
same compound via `cfg.norm_cmpd()`, and the lists may mix notations without
changing `tag`, since everything is normalised to `cmpd_N` first.

Existence, by contrast, is strict: a compound in neither list, in both, or
absent from the challenge csv raises immediately and names itself.

Useful methods: `cfg.is_mixed`, `cfg.arm_map()`, `cfg.arm_of(cmpd)`,
`cfg.arm_table()` (a DataFrame of the resolved assignment with both ids).

---

## 3. Paths

| field | type | default | meaning |
|---|---|---|---|
| `predictions_root` | `str` | `""` | root holding run directories; `from_env` reads `$CYP_PREDICTIONS` |
| `challenge_csv` | `str` | `""` | compound table (id, id2, SMILES); `from_env` reads `$CHALLENGE_CSV` |
| `cif_subdir` | `str` | `"all_predictions"` | co-folding output directory inside the run |
| `refine_subdir` | `str \| None` | `None` | set to reuse an existing refinement directory instead of the tag-derived name |

Derived, all read-only properties:

```
run_dir         = predictions_root / run_name
cif_dir         = run_dir / cif_subdir
pdb_dir         = cif_dir / "pdb_format"
work_dir        = run_dir / "min_work"
refine_dir      = run_dir / f"refined_{arm_label}_{tag}"   (or refine_subdir)
base_stage_dir  = refine_dir / "_suite"
stage_dir       = refine_dir / "_suite_<variant>"          (or _suite)
best_dir        = refine_dir / "best_poses_<variant>"      (or best_poses)
```

`cfg.ensure_dirs()` creates all of them. `cfg.manifest(name)` returns the path
of a stage manifest, routing variant-owned stages to `stage_dir` and the rest
to `base_stage_dir`.

---

## 4. Chains and residue names

| field | type | default |
|---|---|---|
| `prot_chain` | `str` | `"A"` |
| `heme_chain` | `str` | `"B"` |
| `lig_chain` | `str` | `"C"` |
| `lig_resname` | `str` | `"LIG"` |
| `heme_resnames` | `tuple` | `("HEM", "HEC", "HEB")` |
| `lig_resnames` | `tuple` | `("LIG", "LIG0", "LIG1", "UNK", "UNL")` |
| `pose_regex` | `str` | `^(?P<cmpd>cmpd_\d+)_seed_(?P<seed>\d+)_sample_(?P<sample>\d+)_model$` |

`pose_regex` must expose named groups `cmpd`, `seed` and `sample`. OpenFold3
filenames carry a `cmpd_N__` prefix and arbitrary large seeds; stripping the
prefix lets them parse under the default regex.

---

## 5. Refinement (these define `tag`)

| field | type | default | units | meaning |
|---|---|---|---|---|
| `backbone_k` | `float` | `100.0` | kcal/mol/Å² | backbone positional restraint |
| `sidechain_k` | `float` | `2.0` | kcal/mol/Å² | sidechain positional restraint |
| `ligand_k` | `float` | `15.0` | kcal/mol/Å² | ligand positional restraint |
| `ligand_tol` | `float` | `0.30` | Å | flat-bottom free zone before the ligand restraint bites |
| `wall_k` | `float` | `1.0e5` | kJ/mol/nm⁴ | heme wall force constant |
| `wall_rc_default` | `float \| dict \| None` | `None` | nm | general ligand–heme exclusion radius |
| `wall_rc_fe` | `float \| dict \| None` | `None` | nm | exclusion radius at the iron |
| `wall_rc_pyrrole` | `float \| dict \| None` | `None` | nm | exclusion radius at the pyrrole nitrogens |
| `fe_restrain` | `bool \| dict \| None` | `None` | — | apply the Fe coordination restraint |
| `fe_restrain_max_d0` | `float` | `2.35` | Å | upper bound on the restrained Fe–donor distance |
| `fe_window_half` | `float` | `0.10` | Å | half-width of the flat-bottom Fe window |
| `fe_k_kcal` | `float` | `20.0` | kcal/mol/Å² | Fe restraint force constant |
| `max_iter` | `int` | `5000` | — | minimiser iteration cap |
| `min_tol_kj` | `float` | `2.0` | kJ/mol/nm | minimiser convergence tolerance |
| `forcefields` | `tuple` | `("amber/ff14SB.xml", "implicit/obc2.xml")` | — | protein force field |
| `small_molecule_ff` | `str` | `"gaff-2.11"` | — | ligand force field |
| `nagl_model` | `str` | `"openff-gnn-am1bcc-1.0.0.pt"` | — | partial-charge model |

### Heme wall

The wall is a polynomial penalty, so it cannot diverge:

```
E = step(rc − r) · wall_k · (rc² − r²)²
```

### Per-arm defaults

When the corresponding field is `None`, these apply:

| | binders | non-binders |
|---|---|---|
| `wall_rc_default` | 0.26 nm | 0.34 nm |
| `wall_rc_fe` | 0.19 nm | falls back to `wall_rc_default` |
| `wall_rc_pyrrole` | 0.24 nm | falls back to `wall_rc_default` |
| `fe_restrain` | `True` | `False` |

Resolve them with `cfg.wall_rc_default_for(arm)`, `cfg.wall_rc_fe_for(arm)`,
`cfg.wall_rc_pyrrole_for(arm)` and `cfg.fe_restrain_for(arm)`. The
`cfg.arm_*` properties are the single-arm equivalents.

### Stereochemistry and convergence reporting

| field | type | default | meaning |
|---|---|---|---|
| `stereo_strict` | `bool` | `False` | `True` discards a stereo-mismatched pose instead of recording it |
| `warn_force_kj` | `float` | `250.0` | residual-force warning threshold — **reporting only**, not in `REFINE_PARAMS` |

With `stereo_strict = False` the refinement records `stereo_ok`,
`stereo_mismatch` and `stereo_undefined` rather than gating on them. Note that
RDKit's lowercase `r`/`s` are pseudo-asymmetric centres and `?` means
undeterminable, which is not the same as wrong — a 1,4-disubstituted
cyclohexane is achiral and its labels encode only cis/trans.

`warn_force_kj` is set well above `min_tol_kj` on purpose. Stiff backbone
restraints put a systematic floor of roughly 70 kJ/mol/nm under the residual
force, which reflects the restrained protein rather than the pose, so warning
against `min_tol_kj` fires on nearly every structure and buries the outliers
that matter.

---

## 6. Scoring

| field | type | default | meaning |
|---|---|---|---|
| `shell` | `float` | `3.5` | Å around ligand heavy atoms defining the pocket |
| `max_atoms` | `int` | `360` | OMol25 training ceiling; exceeding it sets `over_max_atoms` |
| `score_device` | `str` | `"cuda"` | torch device for OrbMol |

`dE = E(complex) − E(pocket) − E(ligand)` at the pose geometry. This is an
interaction energy, not a binding free energy: no desolvation, no entropy, and
**not comparable between compounds**.

The pocket residue set is the **union over a compound's poses**, so the atom
set is identical within a compound and the energies are comparable across its
poses. Coordinates still come from each pose individually, so sidechains move.

One consequence, learned the hard way: because the pocket is a property of the
pose *set*, pooling poses from two different sources changes the pocket and
therefore changes `dE` for every pose. `dE` is comparable only within a pose
set whose union pocket is coherent.

---

## 7. Validation

| field | type | default | meaning |
|---|---|---|---|
| `pb_config` | `str` | `"dock"` | PoseBusters configuration |
| `pb_critical` | `tuple` | 12 check names | which failures block a pose |

Default `pb_critical`:

```
minimum_distance_to_organic_cofactors, bond_angles, internal_steric_clash,
bond_lengths, sanitization, all_atoms_connected,
double_bond_stereochemistry, tetrahedral_chirality,
minimum_distance_to_protein, minimum_distance_to_inorganic_cofactors,
volume_overlap_with_protein, volume_overlap_with_organic_cofactors
```

`pb_all_pass` reports all 31 checks; `pb_critical_pass` reports only these.
Selection gates on `pb_critical_pass`. The heme is kept in the conditioning
file — without it `minimum_distance_to_organic_cofactors` has no cofactor to
measure.

---

## 8. Selection

| field | type | default | meaning |
|---|---|---|---|
| `de_tol` | `float` | `0.20` | eV; poses within this of the best are treated as tied, and ipTM breaks the tie |
| `drop_positive_de` | `bool` | `True` | discard positive-dE poses when any negative one exists |
| `require_pb_pass` | `bool` | `True` | walk down the ranking until a pose passes the critical checks |
| `prefer_stereo_ok` | `bool` | `True` | rank stereo-correct poses first, falling back rather than submitting nothing |
| `iptm_csv` | `str` | `""` | defaults to `run_dir/harvest_iptm.csv` |

The ipTM loader accepts two layouts: `cmpd,seed,sample,iptm_AC,...`
(ESMFold2) and `Target_ID,Model,seed,iptm_AC,...` (OpenFold3, where `seed` is
a string and the key is derived from the model stem). It takes the first of
`iptm_AC`, `iptm` or `pair_chains_iptm` that is present, so a custom file with
a plain `iptm` column works.

`de_tol` is an absolute tolerance, which is worth keeping in mind: per-compound
dE spreads vary by orders of magnitude, so the same value can tie nearly every
pose for one compound and none for another.

`decided_by` in the select manifest reads `"iptm"` only when ipTM actually
promoted a pose above its dE rank, not merely when a tie group existed.

---

## 9. Environments

| field | type | default | meaning |
|---|---|---|---|
| `env_main` | `str` | `"openmm"` | conda environment for refinement |
| `env_score` | `str` | `"orbmol"` | conda environment for scoring |

When `orb_models` cannot be imported in the current kernel and
`cross_env=True`, the score stage shells out to `env_score` via `conda run`,
passing the config as JSON.

---

## 10. Entry points

```python
cs.run_all(cfg, stages=STAGES, cross_env=True, **kw)
```

`STAGES = ("convert", "refine", "validate", "score", "select", "plot")`.
Returns `{stage_name: DataFrame}`. Per-stage keyword arguments are passed by
stage name:

```python
cs.run_all(cfg,
           stages = ("convert", "refine", "validate", "score", "select"),
           refine = dict(only=("cmpd_0", "cmpd_5")),
           score  = dict(only=("cmpd_0", "cmpd_5")))
```

`cs.run_stage(cfg, name, cross_env=True, **kw)` runs one stage.

### Stage arguments

| function | arguments |
|---|---|
| `stage_convert(cfg, force=False)` | `force` re-converts existing PDBs |
| `stage_refine(cfg, force=False, limit=None, only=None, verbose_first=True)` | `only` restricts to compounds in either notation; `limit` caps the structure count |
| `stage_validate(cfg, paths=None, force=False)` | `paths` validates specific files; merges into the manifest rather than replacing it |
| `stage_score(cfg, force=False, only=None)` | per-pose resume; `force` discards the manifest and rescores |
| `stage_select(cfg, validate_candidates=True)` | `validate_candidates` runs PoseBusters on demand for unvalidated candidates |
| `stage_plot(cfg)` | dE-vs-rank small multiples plus the rank1→rank2 gap |

### Construction and persistence

| | |
|---|---|
| `RunConfig.from_env(run_name, **over)` | reads `$CYP_PREDICTIONS` and `$CHALLENGE_CSV` (loading a `.env` if present); any field may be passed in `**over` |
| `cfg.save(path)` / `RunConfig.load(path)` | JSON round-trip; every stage writes `config.json` beside its outputs |
| `cfg.describe()` | human-readable summary including resolved per-arm radii |
| `dataclasses.replace(cfg, ...)` | derive a variant config; `RunConfig` is a plain dataclass, not frozen |

---

## 11. Analysis helpers

| function | purpose |
|---|---|
| `audit_best(cfg, expect=None, quiet=False)` | exactly one best pose per challenge compound? names what went missing and at which stage |
| `merge_best(dest, *cfgs, expect=None, strict=True)` | union the best poses of several runs into one submission directory; a compound selected by two runs is an error |
| `describe_runs(predictions_root, run_name=None, params=None)` | resolve `refined_mixed_<tag>` directory names back to the settings that produced them |
| `sweep(base, grid, prefix="")` | yield configs over a parameter grid, routing tag-changing and variant-changing keys appropriately |
| `sweep_plan(cfgs, strict=True)` | print where every config will write and refuse collisions — run before a sweep |
| `compare_runs(cfgs)` | one row per config: refinement metrics plus PoseBusters pass rate |
| `compare_selections(cfgs, labels=None)` | which compound picks change across a sweep |
| `compare_best_structures(cfgs, labels=None)` | compare the submitted coordinates, not just the picks |
| `preview_de_tol(cfg, values=(...))` | how the picks would change at other `de_tol`, read-only, no re-scoring |

---

## 12. Worked examples

**A mixed run with the settled parameters**

```python
cfg = cs.RunConfig.from_env(
    run_name   = "esmfold2_cyp3a4_run9",
    engine     = "esmfold2",
    arm        = "mixed",
    binders    = BINDERS,
    nonbinders = NONBINDERS,
    strict_arms = True,
    ligand_k   = 15.0,
    ligand_tol = 0.30,
    wall_k     = 1.0e5,
    wall_rc_default = {"nonbinders": 0.26},
    shell      = 3.5,
    de_tol     = 0.20,
    require_pb_pass = True,
    variant    = "shell3p5",
)
results = cs.run_all(cfg)
```

**Re-scoring at a different shell, reusing the refinement**

```python
import dataclasses
alt = dataclasses.replace(cfg, shell=4.5, variant="shell4p5")
cs.run_all(alt, stages=("score", "select"))
```

`tag` is unchanged, so `refine_dir` and the convert/refine/validate manifests
are shared; only score, select and plot are recomputed.

**A single-arm ablation on the wall radius**

```python
alt = dataclasses.replace(cfg, wall_rc_default={"nonbinders": 0.29})
cs.sweep_plan([cfg, alt])     # confirm they write to different directories
cs.run_all(alt)
```

The dict pins only the non-binder arm; binders keep their 0.26 default. Since
`wall_rc_default` is in `REFINE_PARAMS`, `tag` changes and the run gets its own
`refine_dir` automatically.
