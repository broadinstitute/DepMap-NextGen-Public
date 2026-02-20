# data

Below is a manifest of files that are expected to be in the data folder, categorized by the section that requires them. They can be accessed via Figshare ([this study](https://figshare.com/s/94af655d3cd71f748be3), [DepMap24Q4](https://plus.figshare.com/articles/dataset/DepMap_24Q4_Public/27993248)) or the [UCSC Xena browser](https://xenabrowser.net/datapages/?cohort=TCGA%20Pan-Cancer%20(PANCAN)&removeHub=https%3A%2F%2Fxena.treehouse.gi.ucsc.edu%3A443) ([Gistic2 copy number](https://tcga-xena-hub.s3.us-east-1.amazonaws.com/download/TCGA.PANCAN.sampleMap%2FGistic2_CopyNumber_Gistic2_all_data_by_genes.gz); [Sample metadata](https://tcga-pancan-atlas-hub.s3.us-east-1.amazonaws.com/download/TCGA_phenotype_denseDataOnlyDownload.tsv.gz))

## Section 1

- Model and screen information:
    - `model_metadata.csv`
    - `screen_metadata.csv`
- Achilles outputs for genome-wide CRISPR screen data:
    - `screen_gene_effect.csv`
    - `screen_gene_dependency.csv`
    - `crispr_control_genes.csv`
    - `crispr_naive_gene_score.csv`
    - `crispr_screen_qc.csv`
- Omics data for NextGen DepMap models published in this manuscript:
    - `next_gen_expression.csv`
    - `next_gen_copy_number.csv`
    - `next_gen_damaging.csv`
    - `next_gen_hotspot.csv`
    - `next_gen_omics_signatures.csv`
- Celligner data:
    - `celligner_coordinates.csv`
    - `celligner_distance_matrix.csv`
    - `tumor_expression.csv.gz`
- Gene-level annotations:
    - `geneset_table.csv`
    - `hgnc_table.csv`
    - `paralog_pair_table.csv`
- Minipool screen data:
    - `minipool_guide_metadata.csv`
    - `minipool_library_composition.csv`
    - `minipool_readcounts.csv`
    - `minipool_screen_metadata.csv`
    - `minipool_sequence_metadata.csv`
- Single-KO viability data:
    - `validation_normalized_viability.csv`
- TCGA PanCancer selected mutations in COAD, READ, PAAD, BRCA, and PRAD:
    - `tcga_select_mutations.csv`
- From [DepMap24Q4](https://plus.figshare.com/articles/dataset/DepMap_24Q4_Public/27993248):
    - `public_24q4_copy_number.csv`
    - `public_24q4_damaging.csv`
    - `public_24q4_expression.csv`
    - `public_24q4_hotspot.csv`
    - `public_24q4_omic_signatures.csv`

## Section 2

## Section 3

- `CDK6CN_GBM_HV2.csv`
- `CN_24q4_GBM_HV2.csv`
- `CNS2D3D_1.csv`
- `CNS2D3D_2.csv`
- `CRISPR_24q4v42_mONCTSG_2DS.csv`
- `CRISPR_24q4v42_mONCTSG_3DO.csv`
- `CRISPR_24q4v42_mONCTSG_mf_2DS.csv`
- `CRISPR_24q4v42_mONCTSG_mf_2DS.csv`
- `CRISPR_24q4v42_mONCTSG_pw_2DS.csv`
- `CRISPR_24q4v42_mONCTSG_pw_3DO.csv`
- `CRISPR_WNT_24q4_top6.csv`
- `MPvsSelectiveDep_Org_v31.csv`
- `OrgSelective_24q4v43_wnt3a.csv`
- `repetitive_mut_3DO.csv`
- `nextgen_cns_2D3D_highvar.csv`