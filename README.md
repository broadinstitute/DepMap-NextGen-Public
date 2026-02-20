# depmap-nextgen

# Contents
- [Setup](#setup)
- [Running the code](#running-the-code)

# Setup

## Python

We use [`poetry`](https://python-poetry.org/docs/) (version 2.1.3) for python package management. Documentation can be found at https://python-poetry.org/docs, with intructions for installing poetry and using some common commands. A human-readable version of package requirements for the environment used in this study can be found in `pyproject.toml`.

To recreate the python environment, navigate to the folder containing `poetry.lock` and `pyproject.toml`, and run `poetry install`. Poetry will install all package dependencies, as specified by their exact versions in the `.lock` file, into an environment which can either be activated for interactive use with `poetry shell` ([documentation](https://github.com/python-poetry/poetry-plugin-shell)) or by prefixing python commands with `poetry run`. The latter case is used in the run instructions below.

## R

The R environment used in Section 3 uses R version 4.4.2, with the following additional packages (and their versions) installed: 
- `scales_1.3.0`
- `gplots_3.2.0`
- `RColorBrewer_1.1-3`
- `dendextend_1.19.0`

# Running the code

## Setup and data acquisition:

Clone this repository, and download the contents of the raw data that will be released in the DepMap Portal into the `data` directory.

The below sections will reproduce the files deposited in `processed`.

## Section 1 (python)

___Main Figures 1b-g, 1i, 2a, 2c-h, 3a, 4a-f, 5a-b, 5d-h; Extended Figures 1-c, 6f, 10e, 11a-c, 12a-d___

To reproduce analytical results for the above figure panels, run the following code (runtime ~20 minutes): ```poetry run python src/generate.py```

Subsequently, run the remaining scripts, in any order, to generate these figure panels (runtime ~5 minutes total):

    poetry run python src/achilles_qc.py
    poetry run python src/celligner.py
    poetry run python src/experimental_validation.py
    poetry run python src/expression_addiction.py
    poetry run python src/figure_legends.py
    poetry run python src/glioblastoma_dependencies.py
    poetry run python src/metaprogram_expression.py
    poetry run python src/minipool.py
    poetry run python src/omics_summaries.py
    poetry run python src/onc_tsg_associations.py
    poetry run python src/oncoplot.py
    poetry run python src/organoid_vs_adherent.py
    poetry run python src/paralog_dependencies.py
    poetry run python src/pdac_classical_dependencies.py
    poetry run python src/tda_class.py

The generation step will populate the `processed` directory with `.csv` files, while the subsequent steps will simultaneously populate the `processed` directory with `.csv` files and `figures` directory with `.pdf` files.

## Section 2 (R)

___Main Figures 3b-d; Extended Figure 6___
To reproduce analytical results for the above figure panels, run the following R scripts:
    Rscript Fig3bc_EDFig6ab_GBM_celligner_transcriptional_analysis.R
    Rscript Fig3d_EDFig6c_GBM_GSEA.R
    

## Section 3 (R)

___Main Figures 3g; Extended Figures 3a, 3d, 8a, 8b, 9b-e, 10a-c___ 

To reproduce analytical results for the above figure panels, run the following R scripts:

    Rscript Fig3g_CDK6_dependency_CN_correlation.R 
    Rscript EDFig3a_CNS_spheroid_adherent_hierarchical_clustering.R
    Rscript EDFig3d_CNS_spheroid_adherent_individual_gene_expression.R
    Rscript EDFig8a_PDAC_classical_assoc_dep_in_RSpondin_Wnt3a_CM_models.R 
    Rscript EDFig8b_Wnt_dependency_mutation_association.R                 
    Rscript EDFig9be_MP_associated_dependencies.R   
    Rscript EDFig10ac_2Dvs3D_oncogene_TSG_dependency.R