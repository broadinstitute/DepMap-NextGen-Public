# depmap-nextgen

# Contents
- [Setup](#setup)
- [Running the code](#running-the-code)

# Setup

## Python

We use [`poetry`](https://python-poetry.org/docs/) (version 2.1.3) for python package management. Documentation can be found at https://python-poetry.org/docs, with intructions for installing poetry and using some common commands. A human-readable version of package requirements for the environment used in this study can be found in `pyproject.toml`.

To recreate the python environment, navigate to the folder containing `poetry.lock` and `pyproject.toml`, and run `poetry install`. Poetry will install all package dependencies, as specified by their exact versions in the `.lock` file, into an environment which can either be activated for interactive use with `poetry shell` ([documentation](https://github.com/python-poetry/poetry-plugin-shell)) or by prefixing python commands with `poetry run`. The latter case is used in the run instructions below.

## R

The R environment used in Section 3 uses R version 4.4.3, with the following additional packages (and their versions) installed: 
- `scales_1.4.0`
- `gplots_3.2.0`
- `RColorBrewer_1.1_3`
- `dendextend_1.19.0`
- `dplyr_1.2.0`
- `readr-2.2.0`

# Running the code

## Setup and data acquisition:

Clone this repository, and download the contents of the raw data that will be released in the DepMap Portal into the `data` directory.

The below sections will reproduce the files deposited in `processed`.

## Section 1 (python)

___Main Figures 1b-g, 1i, 2a-h, 3a, 4a-f, 5a-b, 5d-h; Extended Figures 1-c, 6f, 10e, 11a-c, 12a-d___

To reproduce analytical results for the above figure panels, run the following code (runtime ~20 minutes): ```poetry run python src/generate.py```

Subsequently, run the remaining scripts, in any order, to generate these figure panels (runtime ~5 minutes total):

    poetry run python src/achilles_qc.py
    poetry run python src/celligner.py
    poetry run python src/experimental_validation.py
    poetry run python src/expression_addiction.py
    poetry run python src/figure_legends.py
    poetry run python src/metaprogram_expression.py
    poetry run python src/minipool.py
    poetry run python src/omics_summaries.py
    poetry run python src/onc_tsg_associations.py
    poetry run python src/oncoplot.py
    poetry run python src/organoid_vs_adherent.py
    poetry run python src/paralog_dependencies.py
    poetry run python src/pdac_classical_dependencies.py

The generation step will populate the `processed` directory with `.csv` files, while the subsequent steps will simultaneously populate the `processed` directory with `.csv` files and `figures` directory with `.pdf` files.

## Section 2 (R)

___Main Figures 3b-d; Extended Figure 6___
To reproduce analytical results for the above figure panels, run the following R scripts:

    Rscript Fig3bc_EDFig6ab_GBM_celligner_transcriptional_analysis.R
    Rscript Fig3d_EDFig6c_GBM_GSEA.R
    

## Section 3 (R)

___Main Figures 3d, 3f; Extended Figures 5a, 5d, 6b, 8e-g, 11a, 11c, 12a-c, R9a, R13d___ 

To reproduce analytical results for the above figure panels, run the following R scripts:

    Rscript EDFig5d_cds2d3d_expression_volcano.R
    Rscript EDFig6b_cnsorg2d3d_dependency_volcano.R
    Rscript EDFig8e_glial_nextgen_diff_dependency.R
    Rscript EDFig8f_cdk6dep_mut_association.R
    Rscript EDFig11a_pdacc_associated_dep_wnt3a.R
    Rscript Fig3d_glial_mesenchymal_diff_dependency.R

To produce figures, run the following scripts

    Rscript EDFig5a_cns2d3d_expression_dendrogram.R
    Rscript EDFig8g_cdk6dep_cn_association_glial.R
    Rscript EDFig11c_wnt_dependency_mut_association.R
    Rscript EDFig12ac_2d3d_onctsg_dependency.R
    Rscript Fig3f_cdk6dep_cn_association.R
    Rscript FigR9a_depmap_han_2d3d_screens.R
    Rscript FigR13d_cdk6dep_cdkn2acn_per_lineage.R