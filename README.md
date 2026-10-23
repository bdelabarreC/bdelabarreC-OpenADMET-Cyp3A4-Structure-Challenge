This will be where I will host detailed information on my participation in OpenADMET's 2026 Cyp Structure Prediction Challenge
Start date:  27Sept2026
Competition end date:  3Nov2026

Most recently updated competition standings (live table at HF: https://huggingface.co/spaces/openadmet/cyp-challenge)


| Rank | Username | Submitted  | LDDT-PLI | BiSyRMSD | LDDT-LP | Run # | Diff to Top | Note |
|-----:|----------|------------|---------:|--------:|---------:|---------:|--------|----------|
| 1 | TCB | 2026-10-23 | 0.447  | 3.865  | 0.8095 | 5B   | 0.034 | EMSFold2.  OpenMM minimization.  Best of 5 by orbmol/iptm 
| 3 | TCB | 2026-10-20 | 0.3773 | 4.3487 | 0.6521 | 4    | 0.104 | OF3 prediction, Split dataset into heme binders vs non-hemebinders and used different pocket environments.  Scoring by orbmol energy calc only to pick best prediction |
| 5 | TCB | 2026-09-29 | ~0.27  | -      | -      | 1(A) | 0.20 | Direct openfold output best of 4x (7 entries failed QC) |
| 4 | TCB | 2026-09-30 | 0.3506 | 5.5842 | 0.5916 | 1(C) | 0.12 | Output from run 1 put through light OpenMM minimization |
| 2 | TCB | 2026-10-02 | 0.3917 | 4.2177 | 0.6702 | 2(C) | 0.075 | Best of 20 per molecule taken and put through light minimization with explicit solvent model |
| 5 | TCB | 2026-10-22 | 0.275  | 10.2 | 0.507 | 5A | 0.21 | Baseline EMSFold2.  No minimization, scored best of 5 by IPTM (4 failed QC) |


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


A review of the ~120 structures in the PDB of molecules bound to the heme iron showed a huge bias that will be difficult to overcome with weight based tuning.  After removing duplicates, the dataset consists of 94 heme coordinating and only 14 non-heme coordinating.

Finetuning here would be poor returns on the work required to run that.

The bias was made clear by running a modified compound 9 where the pyrdine ring was replaced by a penyl group (ie - coordinating nitrogen removed).  OF3 will place these two compounds into the exact same place, with a slight offset due to the proton on the phenyl ring.  This is a wrong answer but one that is understandable based on the bias in the training set.

So instead - I started using pocket definitions and reworked some of the OF3 code to allow me to create exclusion zones, specifically around the Heme iron atom.

The pocket definitions is based on reviewing all of the experimental structures and picking residues that are within 3.5A of any.  For the 19 non-heme, I had a separate pocket defintion.  I split the challenge compounds into two sets (heme binding (11) and non heme-binding (9).  I did two separate runs generating 20 structures per compound (5 different starting seeds * 4 runs per seed).

I retained the light minimization via OpenMM and added an energy evaluation step via OrbMol.  Unfortunately - the resulting score did not beat the previous efforts!  I used pure OrbMol scoring for selection - it may be worth a revisit to use a combination of OrbMol and IPTM type scoring.

However - I am at least a full tenth below the top scores using this approach, so my next effort will pivot to using ESMFold2 to generate the initial poses.  The difficulty here is that I don't think I can run their best model locally, and have to rely on the default 100 tokens/day or get a waiver to allow bigger runs to compute the poses.  My first effort had me running out of tokens around compound 9 generating just 5 poses per compound.
