This will be where I will host detailed information on my participation in OpenADMET's 2026 Cyp Structure Prediction Challenge
Start date:  27Sept2026
Competition end date:  3Nov2026

Competition update for my entries as of this upload:


| Rank | Username | Submitted  | LDDT-PLI | LDDT-LP | BiSyRMSD | Coverage | Run # | Diff to Top | Median | Note |
|-----:|----------|------------|---------:|--------:|---------:|---------:|-------|-------------|--------|------|
| 8 | TCB | 2026-09-29 | ~0.27  | -      | -      | 1 | 1(A) | 0.16 | 0.39  | Direct openfold output best of 4x (7 entries failed QC) |
| 5 | TCB | 2026-09-30 | 0.3506 | 5.5842 | 0.5916 | 1 | 1(C) | 0.075 | 0.39  | Output from run 1 put through light OpenMM minimization |
| 7 | TCB | 2026-10-02 | 0.3917 | 4.2177 | 0.6702 | 1 | 2(C) | 0.034 | 0.39  | Best of 20 per molecule taken and put through light minimization with explicit solvent model |


Initial thoughts: 

This is a larger protein - more vRAM would be helpful but I will stick to workflows suitable for my NVidia 4070 with 16 Gb vRAM.

The PXR top scoring entries used ESMFold2 but that does not fit within the above criterion.  Will likely stick to OpenFold3 here.

Goals:
Focus on opensource and homebrewed methods.
Improve workflow and documentation!  
Focus on co-folding and improve post co-fold selection.  The PXR structure prediction challenge revealed that I had the correct answer 80-90% of the time but wasn't able to fish it out reliably.
Post-cofold selection strategies were explored in the PXR challenge.  Ultimately resorted to internal scoring as it produced the highest 'checkpoint' scores (LDDT-PLI against 50% of dataset).  Ultimately the final score was not a reflection of this checkpoint score!

Post-folding scoring exploration to focus on:
-MD simulation to identify stable binding poses
-shape matching (Roshambo) against existing structures
-pharmacophore matching

QC checks will be as PoseBusters as before (now an official check in the competition).  Visual inspection will be feasible with only 20 structures.


I have decided to open source my work flows during the competition. 

Feel free to use as you wish within the restrictions of the GNU license.

Much of this workflow was developed in the previous (PXR) challenge but not explicitly shared there - consider this the update to that effort as well.

I am also building upon the opensource efforts of many other people - I will generate a complete list eventually.  For now it will be apparent by perusing the dependencies in the Jupyter notebooks included here.

My local information (directories, usernames, etc) is protected through use of a .env file that will not be uploaded, but an example of how to format one will be.

In some cases I have generated 'githubsafe' versions of files where I have removed directory structures so that your own can be placed there.


File descriptions (work in progress):

### Notebooks & Scripts
* `run_openfold_githubsafe.txt` - How I execute an OpenFold run
* `structure_prediction.ipynb` - Jupyter notebook for processing post-OF3 run
* `structure_refinement.ipynb` - Jupyter notebook for further processing via OpenMM

### Data Files
* `cyp3a4_A.npz` - Alignment file for input to OF3
* `cyp3a_PDB_2026.xlsx` - Current list of human CYP3A4 files in the PDB (via UniProt lookup)
* `cyp3a4_pdb_summary.csv` - Detailed information on PDB CYP3A4 structures
* `validation.py` - Presubmission validation file provided by OpenADMET
* `cyp3a4_challenge_vs_pdb.csv` - Comparison of PDB ligands to the current OpenADMET challenge set
* `cyp_challenge_compound.csv` - OpenADMET challenge compounds
* `openadmet20_query_cyp3A4_githubsafe.json` - JSON file for running OF3

### Auxiliary Files
* `dotenv_example_forgithub` - Example of a `.env` file structure
* `LICENSE` - GNU license file
* `README.md` - This file               




Next steps - an expanded run of 20 poses per compound showed that OF3 is predicting ligand positions within a very narrow range.
Although I had hoped that using the latest .pt file from OF3 would give a wider range of possibilities, it may be necessary to select
a set of .pdb structures and do some fine-tuning to minimize the current bias to place all ligands very close to the heme group.
