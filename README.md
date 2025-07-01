# Automated Dynamic e-Pharmacophore Generation Pipeline

# Overview
This Python-based pipeline implements a fully automated workflow for the generation of dynamic e-pharmacophore models derived from time-resolved protein-ligand complexes, typically obtained from molecular dynamics (MD) simulations. The approach integrates key components of the Schrödinger Suite to reflect the dynamic nature of protein-ligand interactions in pharmacophore modeling.

# Methodology
For each selected frame of the input ensemble, the pipeline executes the following steps:

1. Structure preprocessing and pH-dependent protonation using PrepWizard.
2. Separation of the ligand and receptor components.
3. Calculation of the ligand-centered binding site coordinates.
4. Generation of Glide-compatible receptor grid files.
5. Construction of e-pharmacophore hypotheses based on the ligand's in situ conformation.

**This strategy enables the characterization of pharmacophoric features across conformationally diverse states, capturing temporal variations in ligand binding modes and interaction profiles.

# Key Features
- Integration with Schrödinger's Glide, PrepWizard, and e-Pharmacophore utilities.
- Designed for seamless use with MD trajectories, but also applicable to other aligned ligand-receptor ensembles.
- Support for high-throughput batch processing and parallel execution.
- Organized output with traceable intermediate and final data products.

# Software Requirements
- Python **3.x**
- A working installation of the **Schrödinger Suite**
- Environment variable `SCHRODINGER` must point to the root directory of the Schrödinger installation  

# Input File Requirements
- **Aligned protein-ligand complex files** in Maestro `.mae` format

# Usage Example

$SCHRODINGER/run python3 run_dynamic_epharmacophore.py --start 1 --end 5001 --step 1

# Output
All results produced by the pipeline will be automatically organized into a predefined directory structure. This ensures that intermediate files and final pharmacophore models are stored in a clear and traceable manner.

- Intermediate files such as preprocessed .mae structures, receptor-ligand splits, grid generation data, and logs will be saved under:

<working_directory>/DYNOPHORE_ANALYSIS/PROCESSED_FILES/

- Final pharmacophore hypotheses (.phypo files) generated for each MD frame or input structure will be collected under:

<working_directory>/DYNOPHORE_ANALYSIS/saved_HYPOTHESIS/

***This organized output allows for easy post-processing, analysis, and validation of the dynamic pharmacophore models across the entire ensemble.

# Note
Although primarily intended for MD-based analyses, this pipeline can be employed with any ensemble of structurally aligned protein-ligand complexes to investigate the pharmacophoric implications of binding flexibility and structural heterogeneity.

# Citation
If you use this tool in your academic work, please cite:

Computational Drug Design Center (HITMER), Faculty of Pharmacy, Bahçeşehir University, Istanbul, Turkey
