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

- `Fig3bc_EDFig6ab_GBM_celligner_transcriptional_analysis.R`
- `Fig3d_EDFig6c_GBM_GSEA.R`

# Section 3

## Scripts

- `EDFig11a_pdacc_associated_dep_wnt3a.R`: Computes correlations with PDAC-classical expression among the models grown with Wnt3a
- `EDFig11c_wnt_dependency_mut_association.R`: Computes Wilcoxon rank-sum tests for dependency of Wnt-related genes on genomic alterations
- `EDFig12ac_2d3d_onctsg_dependency.R`: Computes differential dependency across growth formats for oncogenes and tumor suppressor genes
- `EDFig5a_cns2d3d_expression_dendrogram.R`: Clusters model gene expression across growth formats
- `EDFig5d_cds2d3d_expression_volcano.R`: Plots differential model gene expression across growth formats
- `EDFig6b_cnsorg2d3d_dependency_volcano.R`: Computes differential dependencies across growth formats for CNS models
- `EDFig8e_glial_nextgen_diff_dependency.R`: Computes differential dependencies between glial GBMs and other NextGen models
- `EDFig8f_cdk6dep_mut_association.R`: Computes CDK6 dependency associations with genomic alterations
- `EDFig8g_cdk6dep_cn_association_glial.R`: Plots correlations between CDK6 dependency and copy number alterations in glial models
- `Fig3d_glial_mesenchymal_diff_dependency.R`: Computes differential dependency between glial and mesenchymal GBM models
- `Fig3f_cdk6dep_cn_association.R`: Plots correlations between CDK6 dependency and copy number alterations in all GBM models
- `FigR13d_cdk6dep_cdkn2acn_per_lineage.R`: Plots CDK6 dependency in relation to CNS subtype and genomic alterations in CDKN2A
- `FigR9a_depmap_han_2d3d_screens.R`: Clusters the aggregate gene effect scores in different culture conditions DepMap NextGen models and Han et al., 2020
