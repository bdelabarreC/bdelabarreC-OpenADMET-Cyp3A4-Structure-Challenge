# run in the environment where you do the OpenMM refinement
import importlib, sys
print(sys.executable)
for m in ("openmm", "openmmforcefields", "openff.toolkit", "gemmi", "rdkit",
          "posebusters", "orb_models", "torch", "scipy", "pandas"):
    try:
        mod = importlib.import_module(m)
        print(f"  {m:22s} {getattr(mod, '__version__', 'ok')}")
    except Exception as e:
        print(f"  {m:22s} MISSING ({type(e).__name__})")