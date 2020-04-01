"""
SPARK VERSION

This library contains tools for generating stats from the ZINC15 dataset.

Use ChemLoader to load the ZINC15 dataset.

After import data stores in a spark RDD, each containing:
- ZINC ID
- SMILES code
- tranche
- net charge
- pH
- molecular weight
- logP
- reactive
- feature set

Features are one or several of:
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

import logging
import pyspark
from pyspark.sql.functions import udf
from pyspark.sql.types import StringType
from features import Features


class ChemStats:
    """
    This class generates statistical information about the molecules.
    """
    def __init__(self, molecules):
        """
        Pass the chemical data from chemloader.

        :param molecules: molecule data
        :type molecules: pyspark.dataframe
        """
        self.session = (pyspark.sql.SparkSession.builder
                        .appName("molecule-loader").getOrCreate())
        self.molecule_df = molecules

    def count(self):
        """
        Counts the number of elements.

        :return: number of elements
        :rtype: int
        """
        return self.molecule_df.count()

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
    from chemloader import ChemLoader
    logging.basicConfig(level=logging.INFO)
    # molecules = ChemLoader("/mnt/ssd2/small-molecules.tsv")
    mols = ChemLoader("/mnt/ssd2/molecules.tsv").load()
    stats = ChemStats(mols)
    # print(stats.count())
    # print(stats.describe().show())
    print(stats.pretty_features().show(n=10000, truncate=False))
