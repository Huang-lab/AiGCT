"""dbNSFP column groups shared by the extraction and ClinVar steps."""

# Variant identity, transcript/protein mapping and population frequency.
COLUMN_INFO = [
    "#chr",
    "pos(1-based)",
    "ref",
    "alt",
    "aaref",
    "aaalt",
    "aapos",
    "rs_dbSNP",
    "hg19_chr",
    "hg19_pos(1-based)",
    "hg18_chr",
    "hg18_pos(1-based)",
    "genename",
    "Ensembl_geneid",
    "Ensembl_transcriptid",
    "Ensembl_proteinid",
    "gnomAD2.1.1_exomes_controls_AF",
    "gnomAD4.1_joint_AF",
    "VEP_canonical",
    "Uniprot_acc",
]

# ClinVar annotations are only pulled when building the ClinVar benchmark.
CLINVAR_COLUMNS = ["clinvar_clnsig", "clinvar_review"]

# Raw and rank score columns for the 37 evaluated VEPs.
VEP_LIST = [
    "SIFT_score", "SIFT_converted_rankscore",
    "SIFT4G_score", "SIFT4G_converted_rankscore",
    "Polyphen2_HDIV_score", "Polyphen2_HDIV_rankscore",
    "Polyphen2_HVAR_score", "Polyphen2_HVAR_rankscore",
    "MutationTaster_score", "MutationTaster_converted_rankscore",
    "MutationAssessor_score", "MutationAssessor_rankscore",
    "PROVEAN_score", "PROVEAN_converted_rankscore",
    "VEST4_score", "VEST4_rankscore",
    "MetaSVM_score", "MetaSVM_rankscore",
    "MetaLR_score", "MetaLR_rankscore",
    "MetaRNN_score", "MetaRNN_rankscore",
    "M-CAP_score", "M-CAP_rankscore",
    "REVEL_score", "REVEL_rankscore",
    "MutPred_score", "MutPred_rankscore",
    "MVP_score", "MVP_rankscore",
    "gMVP_score", "gMVP_rankscore",
    "MPC_score", "MPC_rankscore",
    "PrimateAI_score", "PrimateAI_rankscore",
    "DEOGEN2_score", "DEOGEN2_rankscore",
    "BayesDel_addAF_score", "BayesDel_addAF_rankscore",
    "BayesDel_noAF_score", "BayesDel_noAF_rankscore",
    "ClinPred_score", "ClinPred_rankscore",
    "LIST-S2_score", "LIST-S2_rankscore",
    "VARITY_R_score", "VARITY_R_rankscore",
    "VARITY_ER_score", "VARITY_ER_rankscore",
    "VARITY_R_LOO_score", "VARITY_R_LOO_rankscore",
    "VARITY_ER_LOO_score", "VARITY_ER_LOO_rankscore",
    "ESM1b_score", "ESM1b_rankscore",
    "AlphaMissense_score", "AlphaMissense_rankscore",
    "PHACTboost_score", "PHACTboost_rankscore",
    "MutFormer_score", "MutFormer_rankscore",
    "MutScore_score", "MutScore_rankscore",
    "CADD_raw", "CADD_raw_rankscore",
    "DANN_score", "DANN_rankscore",
    "fathmm-XF_coding_score", "fathmm-XF_coding_rankscore",
    "Eigen-raw_coding", "Eigen-raw_coding_rankscore",
    "Eigen-PC-raw_coding", "Eigen-PC-raw_coding_rankscore",
]

# Three predictors name their raw score column without a "_score" suffix in
# dbNSFP; rename them on output so every VEP follows the same convention.
VEP_RENAME = {
    "CADD_raw": "CADD_raw_score",
    "Eigen-raw_coding": "Eigen-raw_coding_score",
    "Eigen-PC-raw_coding": "Eigen-PC-raw_coding_score",
}

EXTRACT_COLUMNS = COLUMN_INFO + VEP_LIST

# Columns that dbNSFP stores as ';'-delimited lists parallel to
# Ensembl_transcriptid, and which therefore have to be split apart before a
# single transcript can be selected.
MULTI_TRANSCRIPT_COLUMNS = [
    "aapos",
    "genename",
    "Ensembl_geneid",
    "Ensembl_transcriptid",
    "Ensembl_proteinid",
    "VEP_canonical",
    "SIFT_score",
    "SIFT4G_score",
    "Polyphen2_HDIV_score",
    "Polyphen2_HVAR_score",
    "MutationTaster_score",
    "MutationAssessor_score",
    "PROVEAN_score",
    "VEST4_score",
    "MetaRNN_score",
    "REVEL_score",
    "MutPred_score",
    "MVP_score",
    "gMVP_score",
    "MPC_score",
    "DEOGEN2_score",
    "LIST-S2_score",
    "VARITY_R_score",
    "VARITY_ER_score",
    "VARITY_R_LOO_score",
    "VARITY_ER_LOO_score",
    "ESM1b_score",
    "AlphaMissense_score",
    "PHACTboost_score",
    "Uniprot_acc",
]

# Key identifying a variant across all steps of the pipeline.
VARIANT_KEY = ["#chr", "pos(1-based)", "ref", "alt"]
