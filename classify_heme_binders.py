#!/usr/bin/env python
"""
Split a compound set into heme-iron (Type II) coordinators vs non-coordinators.

(c) Byron DeLaBarre / The Consulting Biochemist 2026

Chemical test, not statistical: does the molecule carry a lone pair that is
  (a) on a donor atom Fe(III) accepts (sp2 aromatic N, basic sp3 N, S, nitrile N),
  (b) NOT delocalised into an amide / sulfonamide / urea / aromatic sextet, and
  (c) not deactivated by ring EWGs or shielded by flanking substituents?

Classes
  BINDER      strong, well-precedented Fe(III) donor -> run WITHOUT iron exclusion
  BORDERLINE  real but weakened donor -> run in BOTH passes, compare
  NON         no usable donor -> run WITH iron exclusion

Usage:  python classify_heme_binders.py compounds.csv [out.csv]
"""
import sys, csv
from rdkit import Chem

# --------------------------------------------------------------------- donors
S = Chem.MolFromSmarts

# free aromatic N: two ring connections, no H, no substituent, neutral.
AZOLE_FREE_N    = [S("[nX2;H0;r5;!$([n+]);!$(n-[O,N])]")]
PYRIDINE_FREE_N = [S("[nX2;H0;r6;!$([n+]);!$(n-[O,N])]")]

# basic sp3 amine: exclude amide/sulfonamide/urea/carbamate N, aniline N,
# N-N / N-O, nitrile, cation, and aromatic N.
BASIC_AMINE = [S(
    "[NX3;!$(N-[C,S,P]=[O,S,N]);!$(N-[C,S]#[N]);!$(N-a);!$(N-[N,O]);"
    "!$([N+]);!$([nX3])]"
)]
NITRILE   = [S("[NX1]#[CX2]")]
THIOETHER = [S("[SX2;!$(S(=O));!$(s)]([#6])[#6]")]
THIOL     = [S("[SX2H1]")]
OXIME     = [S("[CX3]=[NX2]-[OX2H1]"), S("[NX3]-[CX3](=O)-[OX2H1]")]

# 1,2,3-triazole (any substitution) - weaker donor than imidazole/1,2,4-triazole
TZ123 = [S("[n,nH]1[n][n][c][c]1"), S("c1cn[n][n]1")]

# exocyclic heavy substituent on an aromatic carbon (for steric flanking)
AR_SUBST = S("[c;$(c-!@[!#1])]")

# strong ring EWGs that kill pyridine-N basicity
EWG = [S("[c]-!@[F,Cl,Br,I]"), S("[c]-!@[CX2]#[NX1]"),
       S("[c]-!@[N+](=[O])[O-]"), S("[c]-!@C(F)(F)F"),
       S("[c]-!@[SX4](=O)(=O)")]


def _hits(mol, pats):
    out = set()
    for p in pats:
        if p is None:
            continue
        for m in mol.GetSubstructMatches(p):
            out.update(m)
    return out


def _free_ar_n_idx(mol, pats):
    idx = set()
    for p in pats:
        if p is None:
            continue
        for m in mol.GetSubstructMatches(p):
            idx.update(i for i in m
                       if mol.GetAtomWithIdx(i).GetSymbol() == "N")
    return idx


def classify(smiles):
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return "PARSE_FAIL", "rdkit could not parse", ""

    ri = mol.GetRingInfo()
    tz_atoms   = _hits(mol, TZ123)
    ewg_atoms  = _hits(mol, EWG)
    subst_c    = _hits(mol, [AR_SUBST])

    azole_n = _free_ar_n_idx(mol, AZOLE_FREE_N)
    pyr_n   = _free_ar_n_idx(mol, PYRIDINE_FREE_N)

    groups, notes = [], []
    strong = False

    # --- 5-ring free N (imidazole / 1,2,4-triazole / 1,2,3-triazole) ---
    for n in azole_n:
        hindered = all(nb.GetIdx() in subst_c
                       for nb in mol.GetAtomWithIdx(n).GetNeighbors())
        if n in tz_atoms:
            groups.append("1,2,3-triazole sp2 N (weak donor)")
        else:
            groups.append("free azole sp2 N (imidazole / 1,2,4-triazole)")
            if hindered:
                notes.append("azole N flanked by substituted carbons")
            else:
                strong = True

    # --- 6-ring free N (pyridine) ---
    for n in pyr_n:
        ring = next((r for r in ri.AtomRings() if n in r), ())
        deact = bool(set(ring) & ewg_atoms)
        hindered = all(nb.GetIdx() in subst_c
                       for nb in mol.GetAtomWithIdx(n).GetNeighbors())
        groups.append("pyridine-type ring N")
        if deact:
            notes.append("pyridine ring carries a strong EWG (donor weakened)")
        if hindered:
            notes.append("pyridine N flanked by substituted carbons")
        if not deact and not hindered:
            strong = True

    if _hits(mol, BASIC_AMINE):
        groups.append("basic sp3 amine")
        strong = True
    if _hits(mol, NITRILE):
        groups.append("nitrile (weak donor)")
    if _hits(mol, THIOETHER):
        groups.append("thioether S (weak donor)")
    if _hits(mol, THIOL):
        groups.append("thiol S")
        strong = True
    if _hits(mol, OXIME):
        groups.append("oxime / hydroxamate")
        strong = True

    if strong:
        cls = "BINDER"
    elif groups:
        cls = "BORDERLINE"
    else:
        cls = "NON"
    return cls, "; ".join(groups) if groups else "no Fe(III) donor", "; ".join(notes)


if __name__ == "__main__":
    path = sys.argv[1]
    out = sys.argv[2] if len(sys.argv) > 2 else "heme_binder_classes.csv"
    rows = list(csv.DictReader(open(path)))
    id_col = "id" if "id" in rows[0] else list(rows[0])[0]
    smi_col = "SMILES" if "SMILES" in rows[0] else "smiles"

    res = []
    for r in rows:
        if not (r.get(smi_col) or "").strip():
            continue
        cls, why, note = classify(r[smi_col])
        res.append({"id": r[id_col], "SMILES": r[smi_col],
                    "class": cls, "donor_groups": why, "caveats": note})

    with open(out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["id", "SMILES", "class",
                                           "donor_groups", "caveats"])
        w.writeheader(); w.writerows(res)

    for c in ("BINDER", "BORDERLINE", "NON", "PARSE_FAIL"):
        sub = [x for x in res if x["class"] == c]
        if not sub:
            continue
        print(f"\n=== {c}  (n={len(sub)}) ===")
        for x in sub:
            tail = f"   [{x['caveats']}]" if x["caveats"] else ""
            print(f"  {x['id']:<16} {x['donor_groups']}{tail}")
    print(f"\nwrote {out}")
