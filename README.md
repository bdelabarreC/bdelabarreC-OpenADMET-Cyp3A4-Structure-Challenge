# OpenADMET 2026 CYP3A4 Structure Prediction Challenge

Detailed information on my participation in OpenADMET's 2026 CYP3A4 Structure Prediction Challenge.

- **Start date:** 27 Sept 2026
- **Competition end date:** 3 Nov 2026
- **Challenge background:** [OpenADMET CYP3A4 Challenge](https://huggingface.co/spaces/openadmet/cyp-challenge)
- **View my predictions:** [interactive viewer](https://bdelabarrec.github.io/bdelabarreC-OpenADMET-Cyp3A4-Structure-Challenge/viewer/#9A)

## Results

Ranks are checkpoint standings at time of submission (LDDT-PLI against 50% of the dataset). Follow the links for specific runs.

| Rank | Username | Submitted | LDDT-PLI | BiSyRMSD | LDDT-LP | Run # | Diff to Top | Note |
|-----:|----------|-----------|---------:|---------:|--------:|-------|------------:|------|
| * | TCB | 2026-11-02 | 0.379 | 4.414 | 0.821 | [9A](https://bdelabarrec.github.io/bdelabarreC-OpenADMET-Cyp3A4-Structure-Challenge/viewer/#9A) | 0.102 | ESMFold2, all structures folded with same pocket definition but then refined in 2 different groups (heme binder vs non-binder). Selection by OrbMol energy / PoseBusters pass. IPTM used as tie breaker when dE did not distinguish. Run 9A as final submitted. Scored against all 20 compounds. |
| 1 | TCB | 2026-10-25 | 0.457 | 3.786 | 0.8107 | 5D | 0.023 | ESMFold2, all structures folded with same pocket definition but then refined in 2 different groups (heme binder vs non-binder). |
| 1 | TCB | 2026-10-28 | 0.456 | 3.785 | 0.8105 | [9A](https://bdelabarrec.github.io/bdelabarreC-OpenADMET-Cyp3A4-Structure-Challenge/viewer/#9A) | 0.023 | ESMFold2, full structure with 5 subpockets defined. Same settings as 5D with refactored full-stack run, test for reproducibility. Selection by OpenMM → PoseBusters → OrbMol (tie breaker IPTM). |
| 1 | TCB | 2026-10-28 | 0.453 | 3.793 | 0.808 | 10A | 0.028 | Added pocket 6 to expand coverage. Relaxed heme constraint for non-binders and optimized OrbMol energy scoring by extending from 3.5 to 4.5 Å radius. |
| 1 | TCB | 2026-10-29 | 0.454 | 3.786 | 0.810 | 9B | 0.027 | Rerun of 9A using optimized refinement/selection (0.26 wall_fc / 4.5 Å OrbMol sphere). |
| 2 | TCB | 2026-10-23 | 0.447 | 3.865 | 0.8095 | 5B | 0.034 | ESMFold2. OpenMM minimization. Best of 5 by OrbMol/IPTM. 3 failed QC. |
| 3 | TCB | 2026-10-30 | 0.438 | 3.781 | 0.6734 | [2D](https://bdelabarrec.github.io/bdelabarreC-OpenADMET-Cyp3A4-Structure-Challenge/viewer/#2D) | 0.043 | OF3 — optimized full-stack refinement and selection on OpenFold run 2 (more straightforward and higher-scoring run; see run 2C for more background). |
| 3 | TCB | 2026-10-27 | 0.434 | 4.08 | 0.796 | 8A | 0.05 | ESMFold2 starting at AA 21, run across 5 different subpockets (like 5D). OpenMM minimization. Had to pass PoseBusters, then ranking by OrbMol with tiebreaker by IPTM. All chiralities submitted correct; ESMFold2 scrambles. |
| 4 | TCB | 2026-10-27 | 0.431 | 4.100 | 0.7944 | 7A | 0.05 | ESMFold2 starting at AA 21, run across 10 different subpockets plus full pocket. OpenMM minimization. Had to pass PoseBusters, then best by OrbMol energy/IPTM. NB all chiralities submitted correct; ESMFold2 scrambled some. |
| 5 | TCB | 2026-11-01 | 0.417 | 3.942 | 0.6593 | 50A | 0.06 | Combined: pool of 2 best ESMFold2 and OF3 runs. Attempted objective comparison via OrbMol energy and tiebreaker by normalized IPTM. |
| 6 | TCB | 2026-10-02 | 0.3917 | 4.2177 | 0.6702 | 2(C) | 0.075 | OF3: best of 20 per molecule taken and put through light minimization with implicit solvent model added. |
| 7 | TCB | 2026-10-20 | 0.3773 | 4.3487 | 0.6521 | 4 | 0.104 | OF3 prediction. Split dataset into heme binders vs non-heme binders and used different pocket environments. Scoring by OrbMol energy calculation only to pick best prediction. |
| 8 | TCB | 2026-09-30 | 0.3506 | 5.5842 | 0.5916 | 1(C) | 0.12 | OF3 output from run 1 put through light OpenMM minimization. |
| 9 | TCB | 2026-10-22 | 0.275 | 10.2 | 0.507 | 5A | 0.21 | Baseline ESMFold2. No minimization, scored best of 5 by IPTM (4 failed QC). |
| 10 | TCB | 2026-09-29 | ~0.27 | - | - | 1(A) | 0.20 | OF3: direct OpenFold output, best of 4× (7 entries failed QC). |

## Goals and initial thoughts

Initial thoughts at the start of the competition:

This is a larger protein — more vRAM would be helpful, but I will stick to workflows suitable for my NVIDIA 4070 with 16 GB vRAM. *(Update: all ESMFold2 entries were run on the Biohub platform — many thanks to Biohub for granting my request for more tokens.)*

The PXR top-scoring entries used ESMFold2, but that does not fit within the above criterion. Will likely stick to OpenFold3 here. *(Update: OF3 seemed to have topped out — went to ESMFold2 in the last few weeks of the competition. May revisit OF3 with the updated ranking (OpenMM / PoseBusters / OrbMol) created for ESMFold2 output.)*

Goals:

- Focus on open-source and homebrewed methods.
- Improve workflow and documentation.
- Focus on co-folding and improve post-co-fold selection. The PXR structure prediction challenge revealed that I had the correct answer 80–90% of the time but wasn't able to fish it out reliably.
- Post-co-fold selection strategies were explored in the PXR challenge. Ultimately I resorted to internal scoring as it produced the highest checkpoint scores (LDDT-PLI against 50% of the dataset) — but the final score was not a reflection of that checkpoint score.

Post-folding scoring exploration to focus on:

- MD simulation to identify stable binding poses
- Shape matching (Roshambo) against existing structures
- Pharmacophore matching

QC checks will be PoseBusters as before (now an official check in the competition). Visual inspection will be feasible with only 20 structures.

## Licensing and reuse

I have decided to open source my workflows during the competition. Feel free to use them as you wish within the restrictions of the MIT license.

Much of this workflow was developed in the previous (PXR) challenge but not explicitly shared there — consider this the update to that effort as well.

I am also building upon the open-source efforts of many other people. I will generate a complete list eventually; for now it will be apparent by perusing the dependencies in the Jupyter notebooks included here.

My local information (directories, usernames, etc.) is protected through use of a `.env` file that will not be uploaded, but an example of how to format one will be. In some cases I have generated "githubsafe" versions of files where I have removed directory structures so that your own can be placed there.

## File descriptions

*Work in progress.*

### Code / notebooks

* `openfold3_cyp3a4_setup.ipynb` — how I execute an OpenFold run
* `esmfold2_cyp3a4_setup.ipynb` — Jupyter notebook for processing post-OF3 run
* `cyp_suite_driver_v2.ipynb` — post-co-fold processing of `.cif` files moved into a single folder. Currently set up to do either ESMFold or OpenFold3, but in principle can be easily extended to any other co-folder.
* `cyp_suite_pools.ipynb` — distilled specialized (extremely beta) approach to blending output from multiple co-folders
* `structure_prediction_submit.ipynb` — run final check and send submission to OpenADMET via API
* `pose_in_pocket_check.ipynb` — make sure the co-folder didn't place the small molecule outside of the defined pocket
* `classify_heme_binders.py` — automated sorter to identify heme-Fe binder motifs from SMILES codes
* `cofolder_processer.py` — module used by `cyp_suite_driver_v2.ipynb`. Can be generalized to any co-folding output.
* `OpenMM_environmentcheck.py` — ensure that required modules for OpenMM are present in a conda environment
* `validation.py` — presubmission validation, required by `structure_prediction_submit.ipynb` (created by the OpenADMET team)
* `submit_when_open.py` — sleep-saving utility for the competition: make a timed submission to OpenADMET with the latest prediction

### Data files

* `cyp3a4_A.npz` — alignment file for input to OF3
* `cyp3a_PDB_2026.xlsx` — current list of human CYP3A4 files in the PDB (via UniProt lookup)
* `cyp3a4_pdb_summary.csv` — detailed information on PDB CYP3A4 structures
* `cyp3a4_challenge_vs_pdb.csv` — comparison of PDB ligands to the current OpenADMET challenge set
* `cyp_challenge_compound.csv` — OpenADMET challenge compounds
* `openadmet20_query_cyp3A4_githubsafe.json` — example JSON file for running OF3

### Auxiliary files

* `dotenv_example_forgithub` — example of a `.env` file structure
* `LICENSE` — MIT license file
* `README.md` — this file

## Progress log

### OF3 pose diversity and heme bias

An expanded run of 20 poses per compound showed that OF3 is predicting ligand positions within a very narrow range. Although I had hoped that using the latest `.pt` file from OF3 would give a wider range of possibilities, it may be necessary to select a set of `.pdb` structures and do some fine-tuning to reduce the current bias to place all ligands very close to the heme group.

A review of the ~120 CYP3A4 structures in the PDB with molecules bound to the heme iron showed a huge bias that will be difficult to overcome with weight-based tuning. After removing duplicates, the dataset consists of 94 heme-coordinating and only 19 non-heme-coordinating structures. Fine-tuning here would be poor returns on the work required.

The bias was made clear by running a modified compound 9 where the pyridine ring was replaced by a phenyl group (i.e. coordinating nitrogen removed). OF3 places these two compounds in the exact same place, with a slight offset due to the proton on the phenyl ring. This is a wrong answer, but an understandable one given the bias in the training set.

### Pocket definitions and exclusion zones

Instead, I started using pocket definitions and reworked some of the OF3 code to allow me to create exclusion zones, specifically around the heme iron atom.

The pocket definitions are based on reviewing all of the experimental structures and picking residues within 3.5 Å of any ligand. For the 19 non-heme structures I had a separate pocket definition. I split the challenge compounds into two sets — heme binding (11) and non-heme binding (9) — and did two separate runs generating 20 structures per compound (5 different starting seeds × 4 runs per seed).

I retained the light minimization via OpenMM and added an energy evaluation step via OrbMol. Unfortunately the resulting score did not beat the previous efforts. I used pure OrbMol scoring for selection; it may be worth a revisit to use a combination of OrbMol and IPTM-type scoring.

### Pivot to ESMFold2

I am at least a full tenth below the top scores using this approach, so my next effort will pivot to using ESMFold2 to generate the initial poses. The difficulty here is that I don't think I can run their best model locally, and have to rely on the default 100 tokens/day or get a waiver to allow bigger runs. My first effort had me running out of tokens around compound 9, generating just 5 poses per compound.

Optimizing pays off. I combined the OpenMM, PoseBusters and OrbMol scoring into a single routine (will post as soon as I clean up the code a bit). I used a run on ESMFold2 where I generated 5 different subpockets to simulate the random seed option available in OpenFold3. I had tried to create different pockets for heme-bound vs non-heme-bound compounds (based on ligand chemistry), but these results were poorer than simply giving the entire pocket to ESMFold2. The ligand division was still necessary for the subsequent refinement steps in OpenMM, because atoms near the heme-Fe need different treatment depending on whether they are bound or not.

The unified approach made post-co-folding refinement a one-click operation. Because PoseBusters is used early, all structures passed the submission check, so all 20 got scored. My highest result yet — putting me into 4th at time of writing, just 0.024 points off the top entry.

### Tuning the refinement/selection stack

I have probably pushed ESMFold2 as far as I can go; the limited tokens make it difficult to experiment much, and any academics with freer access definitely have an advantage here. *(Update: the ESMFold2 team granted me an increased daily token limit — thank you.)*

Perhaps it is time to let Boltz2 have a try at prediction. I used it for the PXR challenge with limited success, but my methods were still evolving at that point, so another shot would be a useful benchmark so that OF3, ESMFold2 and Boltz2 can get a proper side-by-side comparison. Boltz2 may have to come post-competition — there is less than a week remaining, so focus will be on tweaking the OpenMM and OrbMol settings to see if I can get better resolution in the energies.

Currently I am testing a method to optimize subpocket selection to get maximum coverage of space within the known binding pocket. Previous efforts at simply maximizing subpocket selection through clustering seemed to slightly degrade the output. *(Update: this didn't seem to make a difference — see run 10A.)* Two things were changed: pocket addition and optimization of some of the selection parameters. I will run the optimized selection on run 9A, which has 5 pockets, to test the additional pocket introduced in run 10.

Post-run selection optimization appears to have topped out — I will revisit post-competition with the full dataset in hand.

With pockets and refinement/selection parameters fully optimized, perhaps it is time to revisit OF3 and then give Boltz2 a shot. Or maybe not — it is getting close to the end of the competition. I may try a few more refinement runs with redefined binders/non-binders. There were a few ambiguous calls, and correct categorization of all of them would likely improve the score, since refinement hinges on proximity to the heme Fe (bond vs no bond).

### Applying the stack back to OF3

In the meantime, a nice result using the integrated OpenMM / PoseBusters / OrbMol refinement and selection stack. I was able to boost the OpenFold3 run 2 effort (a straightforward run, 20 tries per compound) from 0.39 (run 2C) to 0.44 (run 2D). The plot of energies and IPTM shows the benefit of generating multiple poses during co-folding and then letting a combination of energy and IPTM sort things out:

<img width="1483" height="1153" alt="Energy vs IPTM selection plot" src="https://github.com/user-attachments/assets/3849adf7-17e9-4da8-befa-e611e07183a3" />

Useful post-co-folding selection is what was missing from my previous PXR efforts. In that case, comparison to the actual structures showed that with enough co-folding attempts the correct structure was there (or at least within 4 Å) over 90% of the time.

This new refinement/selection stack (working name: sweet-sipper) improves both the OF3 and the ESMFold2 scores. It will be interesting to compare it against the released data after the competition.

## Initial conclusions

Final scoring completed — my standing tumbled quite a bit, from the top 10 at close of competition down to 30th. This is a bit disappointing; I did better in the previous PXR challenge.

Congratulations to "dargason" with a clear top-scoring entry. They went all out and took input from every available co-folding platform, generating hundreds of potential poses per ligand, then triaged down to consensus entries along lines similar to what I used here (heme binder type, energetics in the pocket, IPTM-based scoring). We had the same split of Type I / Type II binders — this is easily done by inspecting the compounds. I think my methods were sensible; I may have simply been underpowered in the number of attempts made per ligand.

Taken from dargason's [interactive presentation of their efforts](https://aetherark.com/articles/admet-cyp-structure/): correct answers came mostly from ESMFold2 and RoseTTAFold — but not always. Boil-the-ocean approaches are not my favorite way to do things, but clearly it worked well here.

<img width="733" height="1093" alt="Source of correct answers by prediction package" src="https://github.com/user-attachments/assets/16dbdc6e-92e4-467d-9add-1806174453f3" />

Final scoring of all entries for the 2026 CYP3A4 prediction challenge:

<img width="720" height="442" alt="Final scoring of all challenge entries" src="https://github.com/user-attachments/assets/fba697f9-2bb2-47e2-9c52-930d27af66ae" />

In retrospect, I should have spent more time tuning some of my OpenMM and OrbMol settings on known structures. In practice I only streamlined my methods towards the end of October; in early October I was jumping between multiple notebook cells and copying files through the terminal. Important lesson: streamline and optimize as soon as you can, so that small changes can be reliably tested. I will revisit this when the final answers come out, using them as my hold-out set against parameters tuned to PDB structures.

I had been bullish on co-folding methods based on the PXR results, where I saw the correct answer being found within fewer than 50 poses generated per structure (still hard to pick out the correct one). That theirs required hundreds of poses and the running of a half-dozen prediction packages tells us that diffusion-based co-folding does not work for all targets.

There were over 100 examples of CYP3A4 in the PDB. There was definitely a bias towards Type II binders (bonding via heme-Fe), but there were ~20 Type I structures. The models could perhaps be improved by selecting only well-supported ligands (i.e. good electron density). Reliably predicting ligand position in a protein pocket would appear to still be an unsolved problem — if we can't make reliable predictions with this much data, there is clearly something amiss in our methods.
