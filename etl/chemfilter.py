"""
SPARK VERSION

This library contains tools for selecting molecules based
on specific criteria.

Loads the existing data set with QED and feature mapping
already completed.

After import data stores in a spark RDD, each containing:
- ZINC ID
- SMILES code
- tranche
- net charge
- pH
- molecular weight
- logP
- reactive
- feature map
- QED score

Features map to one or several of:
  aggregator
  biogenic
  endogenous
  fda
  in-man
  investigational
  in-vitro
  in-vivo
  metabolite
  world
"""

import argparse
import logging
import os
import pyspark
from pyspark.sql.functions import udf, round
from pyspark.sql.types import StringType
from chemloader import ChemLoader
from features import Features


class ChemFilter:
    """
    This class pulls the overall molecule data, now that it's pretty.
    """
    def __init__(self, load_path):
        """
        Loads the preprocessed data.

        :param load_path: path containing CSV source data
        :type load_path: str
        """
        self.session = (pyspark.sql.SparkSession.builder
                        .appName("molecule-loader")
                        .getOrCreate())
        self.path = load_path
        self.source_files = os.listdir(load_path)
        self.molecule_df = None

    def describe(self):
        """
        Describes statistical information about the
        whole data set.

        :return: dataframe statistics like min, max, average,...
        :rtype: pyspark.dataframe
        """
        return self.molecule_df.describe()

    def features(self):
        """
        Counts the various feature combinations present
        in the data.

        :return: row of feature map and count
        :rtype: pyspark.dataframe
        """
        return self.molecule_df.groupBy('feature_map').count()

    def group_by_qed(self, partitions):
        """
        Creates a QED score grouping with a resolutiion of
        the passed partition.

        :param partitions: number of groupings to create
        :type partitions: int
        :return: dataframe with 'qed_group' column
        :rtype: pyspark.dataframe
        """
        self.molecule_df = (self.molecule_df
                            .withColumn('qed_group',
                                        round(self.molecule_df['qed']
                                              * partitions)
                                        )
                            )
        return self.molecule_df

    def load(self):
        """
        Loads the data from the source_files

        :return: molecules
        :rtype: pyspark.dataframe
        """
        base_df = None
        for file_name in self.source_files:
            full_name = f"{self.path}/{file_name}"
            logging.info(f"Loading data from {full_name}")
            interim_df = (self.session.read.option("sep", "\t")
                          .csv(full_name, header=True))
            if base_df:
                base_df = base_df.union(interim_df)
            else:
                base_df = interim_df
        self.molecule_df = base_df
        return base_df

    def pretty_features(self):
        """
        Makes a spreadsheet-like list of the feature combinations
        present in the data.

        :return: row of feature map, count, and a 1/0
                 pipe-delimited section suitable for Excel split
        :rtype: pyspark.dataframe
        """
        udf_features = udf(self._discrete_features, StringType())
        feature_df = self.features()
        return feature_df.withColumn("|".join(Features.KNOWN_FEATURES),
                                     udf_features("feature_map"))

    def select(self, command):
        """
        Selects using a SQL command.

        :param command: SQL command
        :type command: str
        :return: selected rows
        :rtype: pyspark.dataframe
        """
        self.molecule_df.createOrReplaceTempView("dfTable")
        result = self.session.sql(command)
        return result

    def write_group(self, group_col=None, path="", suffix=""):
        """
        Writes the current data frame, grouped
        by group_col and with the path and suffix.

        :param group_col: column for grouping
        :type group_col: str
        :param path: base path for writing
        :type path: str
        :param suffix: suffix for base path
        :type suffix: str
        """
        (self.molecule_df.repartition(group_col)
         .write.partitionBy(group_col).csv(f"{path}{suffix}",
                                           sep="\t", header=True)
         )

    @staticmethod
    def _discrete_features(feature_map):
        """
        Helper for UDF to get the feature string from an incoming feature map

        :param feature_map: comma separated feature bits
        :type feature_map: str
        :return: feature str
        :rtype: str
        """
        out_string = ""
        for feature in Features.KNOWN_FEATURES:
            if Features.KNOWN_FEATURES[feature] & feature_map:
                out_string += "|1"
            else:
                out_string += '|0'
        return out_string[1:]


if __name__ == "__main__":
    # pylint: disable=invalid-name
    parser = argparse.ArgumentParser("Load molecules from existing CSVs")
    parser.add_argument('source', metavar="source", type=str,
                        help="Source path")
    parser.add_argument('dest', metavar="dest", type=str,
                        help="Destination base filename")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    mols = ChemFilter(args.source)
    mols.load()
    mols.group_by_qed(2)
    mols.write_group('qed_group', args.dest, "_qed")
    # mol_class.describe().show()
    # filtered = mol_class.select("SELECT * FROM dfTable "
    #                             "WHERE qed BETWEEN 0.39999 AND 0.4")
    # filtered.show()
    # for x in range(0, delta - step, step):
    #     fx = float(x)
    #     filtered = mol_class.select(f"SELECT * FROM dfTable "
    #                                 f"WHERE qed BETWEEN "
    #                                 f"{fx / fd} "
    #                                 f"AND {(fx + fs) / fd}")
    #     (filtered.repartition(1).coalesce(1)
    #      .write.csv(f"{args.dest}_qed_{x}", sep="\t",
    #      header=True))
    # stats = ChemStats(mols)
    # print(stats.count())
    # print(stats.describe().show())
    # print(stats.pretty_features().show(n=10000, truncate=False))
