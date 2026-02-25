# src

Below is a manifest of files included in this directory.

# Section 1

## Scripts

- `generate.py`: Runs all data processing and analysis for Figures 1-2, 3a, 4-5
- `achilles_qc.py`: Generates figures related to screen quality control
- `celligner.py`: Generates figures related to Celligner 
- `experimental_validation.py`: Generates figures related to experimental validation of adhesion genes
- `expression_addiction.py`: Generates figures related to expression addiction biomarkers
- `figure_legends.py`: Generates figure legends
- `metaprogram_expression.py`: Generates figures related to metaprogram expression (from Gavish et al., 2023) and comparisons
- `minipool.py`: Generates figures related to minipool screens
- `omics_summaries.py`: Generates figures summarizing the data available per model
- `onc_tsg_associations.py`: Generates figures related to oncogene or tumor suppressor alterations as biomarkers
- `oncoplot.py`: Generates oncoprints
- `organoid_vs_adherent.py`: Generates figures comparing dependencies in organoid models to dependencies in adherent models
- `paralog_dependencies.py`: Generates figures related to paralog dependency biomarkers
- `pdac_classical_dependencies.py`: Generates figures related to the PDAC-Classical metaprogram (from Gavish et al., 2023)
- `tda_class.py`: Generates figures related to gene dependency classes

## Helper functions, constants, and assets

- `constants.py`: Constants used commonly throughout the project
- `data_utils.py`: Helper functions for wrangling data and performing statistical tests
- `figure_utils.py`: Helper functions for plot-making
- `gene_utils.py`: Helper functions for standardizing gene names
- `assets/custom_legend_elements.py`: Custom markers for oncoprints
- `assets/stylesheet.mplstyle`: Stylesheet for figures

# Section 2

## Scripts

# Section 3

## Scripts

- `EDFig11a_pdacc_associated_dep_wnt3a.R`
- `EDFig11c_wnt_dependency_mut_association.R`
- `EDFig12ac_2d3d_onctsg_dependency.R`
- `EDFig5a_cns2d3d_expression_dendrogram.R`
- `EDFig5d_cds2d3d_expression_volcano.R`
- `EDFig6b_cnsorg2d3d_dependency_volcano.R`
- `EDFig8e_glial_nextgen_diff_dependency.R`
- `EDFig8f_cdk6dep_mut_association.R`
- `EDFig8g_cdk6dep_cn_association_glial.R`
- `Fig3bc_EDFig6ab_GBM_celligner_transcriptional_analysis.R`
- `Fig3d_EDFig6c_GBM_GSEA.R`
- `Fig3d_glial_mesenchymal_diff_dependency.R`
- `Fig3f_cdk6dep_cn_association.R`
- `FigR13d_cdk6dep_cdkn2acn_per_lineage.R`
- `FigR9a_depmap_han_2d3d_screens.R`
