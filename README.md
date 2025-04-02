# synapse-evolution
Evolution of the Tri-PDZ domain in PSD95 (DLG-4 gene).  
Code used in manuscript by Nilkant, et. al.

**Protein Language Modeling**: ESM PCA analysis 
1. Used ESM2-15B (see repo [here](https://github.com/facebookresearch/esm)) to extract embeddings, which were mean-pooled: `python extract.py esm2_t48_15B_UR50D PDZs.fasta output_dir --include mean`
2. Conducted PCA and applied statistical tests using `plotPCA.py` and `python anovay.py`.

**AF2-Multimer Analysis**  
Used localcolabfold's implementation (repo [here](https://github.com/YoshitakaMo/localcolabfold)) to output structure predictions, pLDDTs, and pAEs.

Ligand Interaction Score (LIS):  
We calculated the LIS scores, as described [here](https://doi.org/10.1101/2024.02.19.580970), using `python lis.py` and `python lis_controls.py` and plotted histograms with `python plot.py`.

Hydrogen Bonding Analysis:  
`python analyze_all_hbonds.py`  
This uses gromacs to process files and then MDAnalysis (in get_hbonds.py) to measure h-bonds. 

**Water Structure Hydrophobicity Analysis**  
The water 3-body angle distribution (also called water triplet distribution) analysis has been described by [Monroe & Shell](10.1063/1.5111545). We used scripts from the implementation in [this repo](https://github.com/samlobe/protein_WaterStructure_Hydrophobicity).  
1) Process a pdb file, e.g. `bash process_with_gromacs.sh PDZ.pdb`; it outputs PDZ_processed.gro
2) Running MD simulations with OpenMM using `python simulate_with_openmm.py ...` with unique random seeds (by using the `--random_seed` flag); it outputs traj.dcd
3) Measuring the water angles around the key h-bond donors, e.g. `python triplet.py ../PDZ_processed.gro ../traj.dcd --selection 'resid 14 and name H'`; it outputs e.g. angles/PDZ_resid_14_and_name_H_angles.txt
4) Process the angles to get the triplet distribution, e.g. `python process_angles.py ../PDZ.pdb --oneAnglesFile angles/PDZ_resid_14_and_name_H_angles.txt -o PDZ_resid_14_and_name_H_triplet_data.csv`
5) Integrate the distribution from 100-120° to get the tetrahedral fraction, which is related to hydrophobicity as described by Monroe & Shell and others.