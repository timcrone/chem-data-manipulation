"""
SPARK VERSION

This library contains tools for loading raw data from the ZINC15 dataset.

It is expected that the ZINC15 3D data set is available at ../files.docking.org/3D/*
in the standard prefix/suffix/prefixsuffix.txt format, so for example the file
../files.docking.org/3D/AA/AAHL/AAAAHL.txt should exist and be in the standard
tab-separated format that ZINC uses.

All other files are ignored.

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

import concurrent.futures
import logging
import os
import re
import sys
import time
import pyspark
from pyspark.sql.functions import udf
from pyspark.sql.types import *


class ChemLoader:
    """
    This class loads data from an individual subtree and returns it in the
    specified format.
    """
    COLUMNS = ["smiles", "zinc_id", "prot_id", "files.db2",
               "substance.inchikey", "net_charge", "ph_mod_fk",
               "substance.mwt", "substance.logp", "purchasable",
               "reactive", "features", "tranche_name"]

    def __init__(self, tree):
        """
        Pass the full path of a tree (files.docking.org/3D) to
        start the loader.  All files matching the pattern within the
        tree will be loaded, and control will return synchronously
        to the caller./

        This loader will not work on Windows due to double slashes.

        :param tree: full path of a subtree
        :type tree: str
        """
        self.session = pyspark.sql.SparkSession.builder.appName("molecule-loader").getOrCreate()
        # https://stackoverflow.com/questions/33681487/how-do-i-add-a-new-column-to-a-spark-dataframe-using-pyspark/33683462
        udf_features = udf(self.features, StringType())
        # self.molecule_df = self.session.read.option("sep", "\t").csv(f"{tree}")
        interim_df = self.session.read.option("sep", "\t").csv(f"/mnt/ssd2/molecules.tsv")
        self.molecule_df = interim_df.toDF(*self.COLUMNS).withColumn("feature_map", udf_features("features"))
        basetime = time.time()
        print(self.molecule_df.count())
        print(f"That took {time.time() - basetime}")
        # print(self.features(feature_map=1023))

    @staticmethod
    def features(feature_list=None, feature_map=None):
        """
        Map features to bits or vice versa.  If feature_list is passed
        then the bitmap is returned, if feature_map is passed then the
        string form is returned.
        :param feature_list: list of features in comma-delimited form
        :type feature_list: str
        :param feature_map: bitmap of features
        :type feature_map: int
        :return: bitmap or list of features
        :rtype: int|str
        """
        KNOWN_FEATURES = {
            'aggregator': 1,
            'biogenic': 2,
            'endogenous': 4,
            'fda': 8,
            'in-man': 16,
            'investigational': 32,
            'in-vitro': 64,
            'in-vivo': 128,
            'metabolite': 256,
            'world': 512
        }
        if feature_list:
            _result = 0
            for feature in feature_list.split(','):
                _result |= KNOWN_FEATURES[feature]
        elif feature_map:
            _result = ""
            for feature in KNOWN_FEATURES:
                if feature_map ^ KNOWN_FEATURES[feature]:
                    _result += f"{feature},"
            _result = _result[:-1]
        else:
            return 0
        return _result


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    test = ChemLoader("../../../files.docking.org/3D/")
    # test = ChemLoader("hdfs://molecules")
    # test = ChemLoader("/mnt/ssd2/molecules.tsv")
    # print(test.tranche)
