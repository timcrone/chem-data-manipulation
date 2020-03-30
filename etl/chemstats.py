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
from pyspark.sql.functions import col, when, count
from features import Features


class ChemStats:
    """
    This class generates statistical information about the molecules.
    """
    def __init__(self, molecules):
        """
        Pass the chemical data from chemloader.

        :param molecules: molecule data
        :type molecules: ChemLoader
        """
        self.session = pyspark.sql.SparkSession.builder.appName("molecule-loader").getOrCreate()
        self.molecule_df = molecules.molecule_df

    def count(self):
        return self.molecule_df.count()

    def describe(self):
        return self.molecule_df.describe()

    def features(self):
        return self.molecule_df.groupBy('feature_map').count()


if __name__ == "__main__":
    from chemloader import ChemLoader
    logging.basicConfig(level=logging.INFO)
    # molecules = ChemLoader("/mnt/ssd2/small-molecules.tsv")
    molecules = ChemLoader("/mnt/ssd2/molecules.tsv")
    stats = ChemStats(molecules)
    print(stats.count())
    print(stats.describe().show())
    print(stats.features().show())
