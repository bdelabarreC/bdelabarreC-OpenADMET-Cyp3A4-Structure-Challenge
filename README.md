This will be where I will host detailed information on my participation in OpenADMET's 2026 Cyp Structure Prediction Challenge
Start date:  27Sept2026
Competition end date:  3Nov2026

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

Update:

I have decided to open source my work flows during the competition.  
Feel free to use as you wish - I put a GNU general public license on everything.
I am building upon the opensource efforts of many other people - I will generate a complete list eventually.  For now it will be apparent by perusing the dependencies in the Jupyter notebooks.
My directory structures are protected through use of a .env file that will not be uploaded, but an example of how to format one will be.
If you wish to use this - you will have to fill in your own variables or use the dotenv module along with your own .env file
In some cases I have generated 'githubsafe' versions of files where I have removed directory structures so that your own can be placed there.
Much of this workflow was developed in the previous (PXR) challenge but not explicitely shared there - consider this the update to that effort as well.

File description:

Notebooks / scripts
run_openfold_githubsafe.txt <how I execute an openfold run>
structure_prediction.ipynb <jupyter notebook for processing post OF3 run>
structure_refinement.ipynb <jupyter notebook for further processing via OpenMM>

Data files:
cyp3a4_A.npz <alignment file for input to OF3>                
cyp3a_PDB_2026.xlsx <current list of human CYP3A4 files in the pdb - via Uniprot lookup>         
cyp3a4_pdb_summary.csv <more detailed information on PDB cyp3a4 structures>       
validation.py <presubmission validation file provided by OpenADMET>
cyp3a4_challenge_vs_pdb.csv <comparison of PDB ligands to the current OpenADMET challenge set>  
cyp_challenge_compound.csv  <openADMET challenge compounds>
openadmet20_query_cyp3A4_githubsafe.json  <.json file for running OF3>

Auxilliary files: 
dotenv_example_forgithub <example of a .env file>    
LICENSE <GNU license file>                                   
README.md <this file>                    



Results update as of this upload:
2 submissions:
Rank	Username	Submitted   LDDT-PLI	LDDT-LP	BiSyRMSD    Coverage    Run #	Note

8	    TCB	        2026-09-29  ~0.27						    1           1(A)	Direct openfold output best of 4x (7 entries failed QC)
5	    TCB	        2026-09-30  0.3506	    5.5842	0.5916	    1		    1(C)    Output from run 1 put through light OpenMM minimization 


Next steps - an expanded run of 20 poses per compound showed that OF3 is predicting ligand positions within a very narrow range.
Although I had hoped that using the latest .pt file from OF3 would give a wider range of possibilities, it may be necessary to select
a set of .pdb structures and do some fine-tuning to minimize the current bias to place all ligands very close to the heme group.
