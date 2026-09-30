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
