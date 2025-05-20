"""
Module description here
"""

import os
import pandas as pd
import numpy as np
import re
from ..file_util import create_folder
from ..date_util import now_str_compact
from ..pd_util import filter_dataframe_by_list
from ..repository import (
    DATA_FOLDER,
    TASK_FOLDERS,
    TABLE_DEFS,
    RepoSessionContext
)
from ..util import Config

# DATA_FOLDER = "data1"
# TASK_FOLDERS = [os.path.join(DATA_FOLDER, task) for task in ["CANCER", "ADRD",
#                                                              "CHD", "DDD"]]


COLUMN_NAME_MAP = {
    "chr": "CHROMOSOME",
    "pos": "POSITION",
    "ref": "REFERENCE_NUCLEOTIDE",
    "alt": "ALTERNATE_NUCLEOTIDE",
    "aaref": "REFERENCE_AMINO_ACID",
    "aaalt": "ALTERNATE_AMINO_ACID",
    "aapos": "AMINO_ACID_POSITION",
    "rs_dbSNP": "RS_DBSNP",
    "hg19_chr": "PRIOR_CHROMOSOME",
    "hg19_pos": "PRIOR_POSITION",
    "hg18_chr": "PRIOR_PRIOR_CHROMOSOME",
    "hg18_pos": "PRIOR_PRIOR_POSITION",
    "genename": "GENE_SYMBOL",
    "Ensembl_geneid": "ENSEMBL_GENE_ID",
    "Ensembl_transcriptid": "ENSEMBL_TRANSCRIPT_ID",
    "Ensembl_proteinid": "ENSEMBL_PROTEIN_ID"
}
VEP_COLUMN_LIST = [
    {"CODE": "MAVENAVG", "raw_score": "MAVEN_(average)_score",
     "rank_score": "MAVEN_(average)_score"},
    {"CODE": "MAVEN", "raw_score": "MAVEN_score",
     "rank_score": "MAVEN_score"},
    {"CODE": "SIFT", "raw_score": "SIFT_score",
     "rank_score": "SIFT_converted_rankscore"},
    {"CODE": "SIFT4G", "raw_score": "SIFT4G_score",
     "rank_score": "SIFT4G_converted_rankscore"},
    {"CODE": "POLYP2HDIV", "raw_score": "Polyphen2_HDIV_score",
     "rank_score": "Polyphen2_HDIV_rankscore"},
    {"CODE": "POLYP2HVAR", "raw_score": "Polyphen2_HVAR_score",
     "rank_score": "Polyphen2_HVAR_rankscore"},
    {"CODE": "MUTTASTE", "raw_score": "MutationTaster_score",
     "rank_score": "MutationTaster_converted_rankscore"},
    {"CODE": "MUTASSESS", "raw_score": "MutationAssessor_score",
     "rank_score": "MutationAssessor_rankscore"},
    {"CODE": "PROVEAN", "raw_score": "PROVEAN_score",
     "rank_score": "PROVEAN_converted_rankscore"},
    {"CODE": "VEST4", "raw_score": "VEST4_score",
     "rank_score": "VEST4_rankscore"},
    {"CODE": "METASVM", "raw_score": "MetaSVM_score",
     "rank_score": "MetaSVM_rankscore"},
    {"CODE": "METALR", "raw_score": "MetaLR_score",
     "rank_score": "MetaLR_rankscore"},
    {"CODE": "METARNN", "raw_score": "MetaRNN_score",
     "rank_score": "MetaRNN_rankscore"},
    {"CODE": "MCAP", "raw_score": "M-CAP_score",
     "rank_score": "M-CAP_rankscore"},
    {"CODE": "REVEL", "raw_score": "REVEL_score",
     "rank_score": "REVEL_rankscore"},
    {"CODE": "MUTPRED", "raw_score": "MutPred_score",
     "rank_score": "MutPred_rankscore"},
    {"CODE": "MVP", "raw_score": "MVP_score",
     "rank_score": "MVP_rankscore"},
    {"CODE": "GMVP", "raw_score": "gMVP_score",
     "rank_score": "gMVP_rankscore"},
    {"CODE": "MPC", "raw_score": "MPC_score",
     "rank_score": "MPC_rankscore"},
    {"CODE": "PRIMAI", "raw_score": "PrimateAI_score",
     "rank_score": "PrimateAI_rankscore"},
    {"CODE": "DEOGEN2", "raw_score": "DEOGEN2_score",
     "rank_score": "DEOGEN2_rankscore"},
    {"CODE": "BAYESDAAF", "raw_score": "BayesDel_addAF_score",
     "rank_score": "BayesDel_addAF_rankscore"},
    {"CODE": "BAYESDNAF", "raw_score": "BayesDel_noAF_score",
     "rank_score": "BayesDel_noAF_rankscore"},
    {"CODE": "CLINPRED", "raw_score": "ClinPred_score",
     "rank_score": "ClinPred_rankscore"},
    {"CODE": "LISTS2", "raw_score": "LIST-S2_score",
     "rank_score": "LIST-S2_rankscore"},
    {"CODE": "VAR_R", "raw_score": "VARITY_R_score",
     "rank_score": "VARITY_R_rankscore"},
    {"CODE": "VAR_ER", "raw_score": "VARITY_ER_score",
     "rank_score": "VARITY_ER_rankscore"},
    {"CODE": "VAR_RL", "raw_score": "VARITY_R_LOO_score",
     "rank_score": "VARITY_R_LOO_rankscore"},
    {"CODE": "VAR_ERL", "raw_score": "VARITY_ER_LOO_score",
     "rank_score": "VARITY_ER_LOO_rankscore"},
    {"CODE": "ESM1B", "raw_score": "ESM1b_score",
     "rank_score": "ESM1b_rankscore"},
    {"CODE": "ALPHAM", "raw_score": "AlphaMissense_score",
     "rank_score": "AlphaMissense_rankscore"},
    {"CODE": "PHACTB", "raw_score": "PHACTboost_score",
     "rank_score": "PHACTboost_rankscore"},
    {"CODE": "MUTFORMER", "raw_score": "MutFormer_score",
     "rank_score": "MutFormer_rankscore"},
    {"CODE": "MUTSCORE", "raw_score": "MutScore_score",
     "rank_score": "MutScore_rankscore"},
    {"CODE": "CADD", "raw_score": "CADD_raw_score",
     "rank_score": "CADD_raw_rankscore"},
    {"CODE": "DANN", "raw_score": "DANN_score",
     "rank_score": "DANN_rankscore"},
    {"CODE": "FATHMMXF", "raw_score": "fathmm-XF_coding_score",
     "rank_score": "fathmm-XF_coding_rankscore"},
    {"CODE": "EIGEN", "raw_score": "Eigen-raw_coding_score",
     "rank_score": "Eigen-raw_coding_rankscore"},
    {"CODE": "EIGENPC", "raw_score": "Eigen-PC-raw_coding_score",
     "rank_score": "Eigen-PC-raw_coding_rankscore"},
]
VEP_COLUMN_LIST_OLD = [
    {"CODE": "REVEL", "raw_score": "REVEL_score",
     "rank_score": "REVEL_rankscore"},
    {"CODE": "GMVP", "raw_score": "gMVP_score", 
     "rank_score": "gMVP_rankscore"},
    {"CODE": "VAR_R", "raw_score": "VARITY_R_score",
     "rank_score": "VARITY_R_rankscore"},
    {"CODE": "VAR_ER", "raw_score": "VARITY_ER_score",
     "rank_score": "VARITY_ER_rankscore"},
    {"CODE": "VAR_RL", "raw_score": "VARITY_R_LOO_score",
     "rank_score": "VARITY_R_LOO_rankscore"},
    {"CODE": "VAR_ERL", "raw_score": "VARITY_ER_LOO_score",
     "rank_score": "VARITY_ER_LOO_rankscore"},
    {"CODE": "ESM1B", "raw_score": "ESM1b_score",
     "rank_score": "ESM1b_rankscore"},
    {"CODE": "EVE", "raw_score": "EVE_score",
     "rank_score": "EVE_rankscore"},
    {"CODE": "ALPHAM", "raw_score": "AlphaMissense_score",
     "rank_score": "AlphaMissense_rankscore"},
    {"CODE": "SIFT", "raw_score": "SIFT_score",
     "rank_score": "SIFT_converted_rankscore"},
    {"CODE": "SIFT4G", "raw_score": "SIFT4G_score",
     "rank_score": "SIFT4G_converted_rankscore"},
    {"CODE": "POLYP2HDIV", "raw_score": "Polyphen2_HDIV_score",
     "rank_score": "Polyphen2_HDIV_rankscore"},
    {"CODE": "POLYP2HVAR", "raw_score": "Polyphen2_HVAR_score",
     "rank_score": "Polyphen2_HVAR_rankscore"},
    {"CODE": "LRT", "raw_score": "LRT_score",
     "rank_score": "LRT_converted_rankscore"},
    {"CODE": "MUTTASTE", "raw_score": "MutationTaster_score",
     "rank_score": "MutationTaster_converted_rankscore"},
    {"CODE": "MUTASSESS", "raw_score": "MutationAssessor_score",
     "rank_score": "MutationAssessor_rankscore"},
    {"CODE": "FATHMM", "raw_score": "FATHMM_score",
     "rank_score": "FATHMM_converted_rankscore"},
    {"CODE": "PROVEAN", "raw_score": "PROVEAN_score",
     "rank_score": "PROVEAN_converted_rankscore"},
    {"CODE": "VEST4", "raw_score": "VEST4_score",
     "rank_score": "VEST4_rankscore"},
    {"CODE": "METASVM", "raw_score": "MetaSVM_score",
     "rank_score": "MetaSVM_rankscore"},
    {"CODE": "METALR", "raw_score": "MetaLR_score",
     "rank_score": "MetaLR_rankscore"},
    {"CODE": "METARNN", "raw_score": "MetaRNN_score",
     "rank_score": "MetaRNN_rankscore"},
    {"CODE": "MCAP", "raw_score": "M_CAP_score",
     "rank_score": "M_CAP_rankscore"},
    {"CODE": "MUTPRED", "raw_score": "MutPred_score",
     "rank_score": "MutPred_rankscore"},
    {"CODE": "MVP", "raw_score": "MVP_score",
     "rank_score": "MVP_rankscore"},
    {"CODE": "MPC", "raw_score": "MPC_score",
     "rank_score": "MPC_rankscore"},
    {"CODE": "PRIMAI", "raw_score": "PrimateAI_score",
     "rank_score": "PrimateAI_rankscore"},
    {"CODE": "DEOGEN2", "raw_score": "DEOGEN2_score",
     "rank_score": "DEOGEN2_rankscore"},
    {"CODE": "BAYESDAAF", "raw_score": "BayesDel_addAF",
     "rank_score": "BayesDel_addAF_rankscore"},
    {"CODE": "BAYESDNAF", "raw_score": "BayesDel_noAF_score",
     "rank_score": "BayesDel_noAF_rankscore"},
    {"CODE": "CLINPRED", "raw_score": "ClinPred_score",
     "rank_score": "ClinPred_rankscore"},
    {"CODE": "LISTS2", "raw_score": "LIST_S2_score",
     "rank_score": "LIST_S2_rankscore"},
    {"CODE": "CADD", "raw_score": "CADD_raw",
     "rank_score": "CADD_raw_rankscore"},
    {"CODE": "DANN", "raw_score": "DANN_score",
     "rank_score": "DANN_rankscore"},
    {"CODE": "FATHMMMLK", "raw_score": "fathmm_MKL_coding_score",
     "rank_score": "fathmm_MKL_coding_rankscore"},
    {"CODE": "FATHMMXF", "raw_score": "fathmm_XF_coding_score",
     "rank_score": "fathmm_XF_coding_rankscore"},
    {"CODE": "EIGEN", "raw_score": "Eigen_raw_coding",
     "rank_score": "Eigen_raw_coding_rankscore"},
    {"CODE": "EIGENPC", "raw_score": "Eigen_PC_raw_coding",
     "rank_score": "Eigen_PC_raw_coding_rankscore"},
    {"CODE": "GENCANYON", "raw_score": "GenoCanyon_score",
     "rank_score": "GenoCanyon_rankscore"},
    {"CODE": "INTFITCONS", "raw_score": "Integrated_fitCons_score",
     "rank_score": "Integrated_fitCons_rankscore"},
    {"CODE": "GM12878FC", "raw_score": "GM12878_fitCons_score",
     "rank_score": "GM12878_fitCons_rankscore"},
    {"CODE": "H1HESCFC", "raw_score": "H1_hESC_fitCons_score",
     "rank_score": "H1_hESC_fitCons_rankscore"},
    {"CODE": "HUVECFC", "raw_score": "HUVEC_fitCons_score",
     "rank_score": "HUVEC_fitCons_rankscore"},
    {"CODE": "LINSIGHT", "raw_score": "LINSIGHT",
     "rank_score": "LINSIGHT_rankscore"},
    {"CODE": "GERPRS", "raw_score": "GERP_RS",
     "rank_score": "GERP_RS_rankscore"},
    {"CODE": "PHP100VERT", "raw_score": "phyloP100way_vertebrate",
     "rank_score": "phyloP100way_vertebrate_rankscore"},
    {"CODE": "PHP470MAM", "raw_score": "phyloP470way_mammalian",
     "rank_score": "phyloP470way_mammalian_rankscore"},
    {"CODE": "PHP17PRI", "raw_score": "phyloP17way_primate",
     "rank_score": "phyloP17way_primate_rankscore"},
    {"CODE": "PHC100VERT",
     "raw_score": "phastCons100way_vertebrate",
     "rank_score": "phastCons100way_vertebrate_rankscore"},
    {"CODE": "PHC470MAM",
     "raw_score": "phastCons470way_mammalian",
     "rank_score": "phastCons470way_mammalian_rankscore"},
    {"CODE": "PHC17PRI", "raw_score": "phastCons17way_primate",
     "rank_score": "phastCons17way_primate_rankscore"},
    {"CODE": "SIPHY29LO", "raw_score": "SiPhy_29way_logOdds",
     "rank_score": "SiPhy_29way_logOdds_rankscore"},
    {"CODE": "BSTAT", "raw_score": "bStatistic",
     "rank_score": "bStatistic_converted_rankscore"}
]
VARIANT_EFFECT_SOURCE_DATA = [
    ["REVEL", "REVEL", "VEP", "REVEL"],
    ["GVMP", "gVMP", "VEP", "gVMP"],
    ["VAR_R", "VARITY_R", "VEP", "VARITY_R"],
    ["VAR_ER", "VARITY_ER", "VEP", "VARITY_ER"],
    ["VAR_RL", "VARITY_R_LOO", "VEP", "VARITY_R_LOO"],
    ["VAR_ERL", "VARITY_ER_LOO", "VEP", "VARITY_ER_LOO"],
    ["ESM1b", "ESM1b", "VEP", "ESM1b"],
    ["EVE", "EVE", "VEP", "EVE"],
    ["ALPHAM", "AlphaMissense", "VEP", "AlphaMissense"],
]
VARIANT_DATA_SOURCE_DATA = [
    ["GNOMGE", "GNOMAD_GENOMES", "GNOMAD GENOMES"],
    ["GNOMEX", "GNOMAD_EXOMES", "GNOMAD EXOMES"]
]


class RepositoryLoader:

    def __init__(self, config: Config,
                 repo_context: RepoSessionContext):
        self._log_folder = config.log.dir
        self._repo_context = repo_context

    def _convert_dot_to_nan(self, val):
        if val == '.':
            return np.nan
        return val

    def _derive_variant_effect_source_columns(self, row):
        source_name = re.match("(.+)_rankscore", row["rank_score"]).group(1)
        return [row["CODE"], source_name, "VEP", source_name]

    def _task_full_path_name(self, task: str, file_name: str):
        return os.path.join(DATA_FOLDER, task, file_name)

    def init_variant_task(self):

        create_folder(DATA_FOLDER)
        for task_folder in TASK_FOLDERS:
            create_folder(task_folder)

        variant_effect_task_df = pd.DataFrame(
            data=np.array([['CANCER', 'CANCER', 'Cancer', 'Cancer']]),
            columns=TABLE_DEFS["VARIANT_TASK"].columns
            )
        variant_effect_task_df.to_csv(os.path.join(DATA_FOLDER,
                                      "variant_task.csv"),
                                      index=False)

    def init_variant_effect_source(self):

        variant_effect_source_df = pd.DataFrame(
            data=VEP_COLUMN_LIST)
        """
            columns=['VEP', 'RAW_SCORE', 'RANK_SCORE']
        )
        """
        ves_columns = TABLE_DEFS["VARIANT_EFFECT_SOURCE"].columns
        ves_file_name = TABLE_DEFS["VARIANT_EFFECT_SOURCE"].file_name
        variant_effect_source_df[ves_columns] =\
            variant_effect_source_df.apply(
                                self._derive_variant_effect_source_columns,
                                axis=1, result_type="expand")

        variant_effect_source_df[ves_columns].to_csv(
            os.path.join(DATA_FOLDER,
                         "variant_effect_source.csv"),
            index=False)

        vds_columns = TABLE_DEFS["VARIANT_DATA_SOURCE"].columns
        vds_file_name = TABLE_DEFS["VARIANT_DATA_SOURCE"].file_name
        variant_data_source_df = pd.DataFrame(
            data=np.array(VARIANT_DATA_SOURCE_DATA),
            columns=vds_columns
            )
        variant_data_source_df.to_csv(os.path.join(DATA_FOLDER,
                                      vds_file_name),
                                      index=False)

    def _build_excep_where_clause(self, column_list: list[str],
                                  suffixes: list[str]):
        """
        Builds a where clause to be used in a DataFrame.query method
        where it checks for inequality between any of the columns
        in the dataframe. For each column in column_list it constructs
        a comparison clause where suffixes[0] is appended to the column
        name on the left side of the comparison and suffixes[1] is
        appended to the column name on the right side of the comparison.

        Parameters
        ----------
        column_list : list(str)
            List of column names to compare.
        suffixes : list(str)
            A list of 2 suffixes with first suffix to be appended to each column
            for left side of comparison and second suffix to be appended to
            column name on right side
        """

        where_list = ["((" + col + suffixes[0] + ".notna() or " +
                      col + suffixes[1] + ".notna()) & " +
                      col + suffixes[0] + " != " + col + suffixes[1] + ")"
                      for col in column_list]
        return " or ".join(where_list)

    def _excep_file_full_path_name(self, task: str, repo_file_name: str):
        os.path.join(self._log_folder, task + "_" +
                     repo_file_name.removesuffix(".csv") + "_" +
                     now_str_compact() + ".csv")

    def _upsert_repository_file(self, new_data: pd.DataFrame, task: str,
                                columns: list[str], repo_file_name: str,
                                pk_columns: list[str]):
        """
        General function for updating one of the repository data files
        with new data.

        To update the files we call the _upsert_repository_file method.
        This method first checks if the row already exists in the file.
        If it doesn't exist it adds the row. If it does exist it updates
        the existing row with the new values.

        Parameters
        ----------
        new_data : pd.DataFrame
            DataFrame containing new data to be loaded.
        task : str
        columns : list(str)
            List of columns in new_data dataframe and in repository data
            file. The columns in the data file are inserted into or updated
            from the columns in the new_data dataframe.
        repo_file_name : str
            Name of repository data file to be inserted/updated.
        pk_columns: list(str)
            List of column names in both new_data and repo_file_name that
            uniquely identify a row. We determine if a row in new_data
            already exists in repo_file_name by using the values in this
            combination of columns to look up a row in repo file.
        """

        repo_file = self._task_full_path_name(task, repo_file_name)
        # Create a new dataframe selecting only the required columns
        new_data_df = new_data[columns]
        # If the repo file doesn't exist we simply populate it from
        # the new_data_df and return.
        if not os.path.exists(repo_file):
            new_data_df.to_csv(repo_file, index=False)
            return
        # Read repository file into a dataframe
        repo_df = pd.read_csv(repo_file)

        # Do a merge between data in repo data frame and new data
        # to find what rows already exist in the repo file using the
        # pk_columns
        df_merge = repo_df.merge(new_data_df, how="inner",
                                 on=pk_columns)
        non_pk_columns = list(set(columns) - set(pk_columns))

        # Build a where clause to find existing rows in repo file
        # where the values of the non pk columns in new data differ
        # from the values in the repo file. We then write these
        # records to an exceptions file.
        exception_where = self._build_excep_where_clause(non_pk_columns,
                                                         ["_x", "_y"])
        exception_df = df_merge.query(exception_where)
        exception_df.to_csv(
            self._excep_file_full_path_name(
                task, "upsert_" + repo_file_name + "_exceptions"),
            index=False)

        """
        Another less concise way to do an upsert

        repo_df.update(repo_df[pk_columns]
                        .merge(new_data_df, on=pk_columns, how="left"))
        left_join_df = new_data_df.merge(repo_df, how="left",
                                on=pk_columns,
                                indicator=True, suffixes=[None, "_yyyy"])
        repo_df = pd.concat(repo_df,
                            left_join_df[left_join_df["_merge"] == "left_only"]
                            [columns])
        """

        # The combine_first method will do following:
        # For each row in new_data_df that doesn't exist in repo_df based
        # upon the values in pk_columns, it will insert that row into
        # repo_df. For each row in new_data_df that does exist in repo_df
        # based on matching values in pk_columns, it will update the
        # row in repo_df with the values from new_data_df.
        # The pk_columns must be the index of both data frames.

        print('here')
        print(repo_file_name)
        print(pk_columns)
        reg = repo_df[pk_columns].groupby(pk_columns).size()

        neg = new_data_df.groupby(pk_columns).size()
        print(len(reg[reg>1]))
        print(len(neg[neg>1]))

        new_data_df = new_data_df.drop_duplicates()
        new_data_df_grouped = new_data_df.groupby(pk_columns).size()
        new_data_dup_pks = new_data_df_grouped[new_data_df_grouped > 1]
        if len(new_data_dup_pks) > 0:
            new_data_df = filter_dataframe_by_list(
                new_data_df, new_data_dup_pks.reset_index(),
                pk_columns, in_list=False)

        neg = new_data_df.groupby(pk_columns).size()
        print(len(neg[neg>1]))

        repo_df.set_index(pk_columns, inplace=True)
        new_data_df = new_data_df.set_index(pk_columns)
        repo_df = repo_df.combine_first(new_data_df)

        reg = repo_df.reset_index()[pk_columns].groupby(pk_columns).size()
        print(len(reg[reg>1]))

        repo_df.to_csv(repo_file)

    def load_variant_file(self, genome_assembly: str, task: str,
                          data_file: str,
                          file_folder: str, data_source: str,
                          binary_label: int, prior_genome_assembly: str,
                          prior_prior_genome_assembly: str):

        """
        Function for loading data from a data file containing data as it is
        downloaded from a source data site into our platform repository data
        files. The input data_file is assumed to contain one row per variant
        along with the label. There will be separate column in that row for
        each vep score. For each row in the input data_file we populate
        the following files:

        - variant.csv - We create one row.
        - variant_effect_label.csv - We create one row with the label and
            other informational columns.
        - variant_effect_score.csv - We create one row for each vep score
            column. So if we have 5 vep score columns we would create 5
            rows in this file.

        To update the files we call the _upsert_repository_file method.
        This method first checks if the row already exists in the file.
        If it doesn't exist it adds the row. If it does exist it updates
        the existing row with the new values.

        Parameters
        ----------
        genome_assembly : str
            Genome assembly, typically hg38
        task : str
            task code
        data_file : str
            File containing data to be loaded.
        file_folder : str
            Location of data_file
        data_source: str
            Source of the input data_file. i.e. HOTSPOT
        binary_label: int
            1 or 0. This is the binary label to be assigned to all the
            variants in the data_file. The assumption is that all of the
            variants in the file have the same label.
        prior_genome_assembly : str
            Genome assembly prior to genome_assembly that we have chromosome,
            position data for. typically hg19
        prior_prior_genome_assembly : str
            Genome assembly prior to prior_genome_assembly that we have,
            chromsome position data for. typically hg18
        """

        variant_df = pd.read_csv(os.path.join(file_folder, data_file))
        variant_df['GENOME_ASSEMBLY'] = genome_assembly
        variant_df["RAW_LABEL"] = np.nan
        variant_df["LABEL_SOURCE"] = data_source
        variant_df["BINARY_LABEL"] = binary_label
        if "gnomAD_exomes_AF" in variant_df.columns:
            variant_df[["ALLELE_FREQUENCY", "ALLELE_FREQUENCY_SOURCE"]]\
                = variant_df.apply(
                lambda row: [row["gnomAD_exomes_AF"], "GNOMEX"] if
                not row["gnomAD_exomes_AF"].isnull() else
                [row["gnomAD_genomes_AF"], "GNOMGE"] if
                not row["gnomAD_genomes_AF"].isnull() else
                [np.nan, np.nan], axis=1, result_type="expand")
        else:
            variant_df[["ALLELE_FREQUENCY", "ALLELE_FREQUENCY_SOURCE"]]\
                = [np.nan, np.nan]

        variant_df = variant_df.map(self._convert_dot_to_nan)
        variant_df.rename(columns=COLUMN_NAME_MAP, inplace=True)
        variant_df["PRIOR_GENOME_ASSEMBLY"] = np.where(
            variant_df['PRIOR_CHROMOSOME'].isnull(),
            np.nan, prior_genome_assembly)
        variant_df["PRIOR_PRIOR_GENOME_ASSEMBLY"] = np.where(
            variant_df['PRIOR_PRIOR_CHROMOSOME'].isnull(),
            np.nan, prior_prior_genome_assembly)
        self._upsert_repository_file(variant_df, "",  # task,
                                     TABLE_DEFS["VARIANT"].columns,
                                     "variant.csv",
                                     TABLE_DEFS["VARIANT"].pk_columns)

        """
        Create dataframe for populating variant_effect_score.csv file.
        For each vep column we create a row in the variant_effect_score_df
        dataframe.
        """
        variant_effect_score_df = pd.DataFrame(
            columns=TABLE_DEFS["VARIANT_EFFECT_SCORE"].columns)
        for vep_columns in VEP_COLUMN_LIST:
            if vep_columns["raw_score"] not in variant_df.columns:
                continue
            vep_df = variant_df.query(vep_columns["rank_score"] +
                                      ".notna()")
            if len(vep_df) == 0:
                continue
            vep_df = vep_df[TABLE_DEFS["VARIANT_EFFECT_SCORE"].columns[:5] +
                            [vep_columns["raw_score"],
                             vep_columns["rank_score"]]]
            vep_df["SCORE_SOURCE"] = vep_columns["CODE"]
            vep_df.rename(columns={vep_columns["raw_score"]: "RAW_SCORE",
                                   vep_columns["rank_score"]: "RANK_SCORE"},
                          inplace=True)
            variant_effect_score_df = pd.concat(
                [variant_effect_score_df,
                 vep_df[TABLE_DEFS["VARIANT_EFFECT_SCORE"].columns]])

        self._upsert_repository_file(
            variant_effect_score_df, task,
            TABLE_DEFS["VARIANT_EFFECT_SCORE"].columns,
            "variant_effect_score.csv",
            TABLE_DEFS["VARIANT_EFFECT_SCORE"].pk_columns)
        self._upsert_repository_file(
            variant_df, task,
            TABLE_DEFS["VARIANT_EFFECT_LABEL"].columns,
            "variant_effect_label.csv",
            TABLE_DEFS["VARIANT_EFFECT_LABEL"].pk_columns)
