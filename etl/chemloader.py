"""
SPARK VERSION

This library contains tools for loading raw data from the ZINC15 dataset.

It is expected that the ZINC15 3D data set is in a file,
based on the data available at ../files.docking.org/3D/*
in the standard prefix/suffix/prefixsuffix.txt format.
For example the file
../files.docking.org/3D/AA/AAHL/AAAAHL.txt should exist
and be in the standard tab-separated format that ZINC
uses.  These files need to be collated into a single
file without headers (i.e. find -exec tail -n +1 {}) as
a prerequisite to calling this loader, for performance
reasons this seems to be the fastest way by far even when
comparing to HDFS or asynchronous file loading from the
filesystem.

All other files should be ignored.

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
import time
import pyspark
from pyspark.sql.functions import udf
from pyspark.sql.types import IntegerType
from features import Features


class ChemLoader:
    """
    This class loads data from an individual file and returns it in the
    specified format.
    """
    COLUMNS = ["smiles", "zinc_id", "prot_id", "files.db2",
               "substance.inchikey", "net_charge", "ph_mod_fk",
               "substance.mwt", "substance.logp", "purchasable",
               "reactive", "features", "tranche_name"]

    def __init__(self, source_file):
        """
        Pass the full path of a file (files.docking.org/3D) to
        start the loader.

        :param source_file: full path of a subtree
        :type source_file: str
        """
        self.session = (pyspark.sql.SparkSession.builder
                        .appName("molecule-loader").getOrCreate())
        self.source = source_file
        self.molecule_df = None

    def load(self):
        """
        Loads the compound data from the instance source file.

        :return: compound data frame
        :rtype: pyspark.dataframe
        """
        udf_features = udf(self._features, IntegerType())
        interim_df = self.session.read.option("sep", "\t").csv(self.source)
        self.molecule_df = (interim_df.toDF(*self.COLUMNS)
                            .withColumn("feature_map",
                                        udf_features("features")))
        basetime = time.time()
        logging.info("Initialization time: %s seconds",
                     time.time() - basetime)
        return self.molecule_df

    @staticmethod
    def _features(feature_string):
        """
        Helper for UDF to get the feature map from an incoming feature string

        :param feature_string: feature string, comma separated
        :type feature_string: str
        :return: feature map
        :rtype: int
        """
        if feature_string is None:
            return 0
        return Features(feature_string).map()


if __name__ == "__main__":
    # pylint: disable=invalid-name
    logging.basicConfig(level=logging.INFO)
    test = ChemLoader("/mnt/ssd2/molecules.tsv").load()
