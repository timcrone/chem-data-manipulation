"""
MEMORY/SWAP VERSION

This library contains tools for loading raw data from the ZINC15 dataset.

It is expected that the ZINC15 3D data set is available at ../files.docking.org/3D/*
in the standard prefix/suffix/prefixsuffix.txt format, so for example the file
../files.docking.org/3D/AA/AAHL/AAAAHL.txt should exist and be in the standard
tab-separated format that ZINC uses.

All other files are ignored.

After import data stores in a list of dicts, each containing:
- ZINC ID
- SMILES code
- tranche
- net charge
- pH
- molecular weight
- logP
- reactive
- feature bitmap

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
    def __init__(self, tree):
        """
        Pass the full path of a tree (files.docking.org/3D) to
        start the loader.  All files matching the pattern within the
        tree will be loaded, and control will return synchronously
        to the caller.

        This loader will not work on Windows due to double slashes.

        :param tree: full path of a subtree
        :type tree: str
        """
        basetime = time.time()  # Let's figure out how long it takes to load
        self.molecule_df = None
        _tranche_list = []
        _tranche_files = []
        root_spec = re.compile("[A-Z]{2}")
        path_spec = re.compile("[A-Z]{4}")
        file_spec = re.compile("([A-Z]{6}).txt")
        for root_path in os.listdir(tree):
            if not self.path_match(root_path, root_spec):
                continue
            for mid_path in os.listdir(f"{tree}/{root_path}"):
                if not self.path_match(mid_path, path_spec):
                    continue
                for file_name in os.listdir(f"{tree}/{root_path}/{mid_path}"):
                    if not self.path_match(file_name, file_spec):
                        continue
                    _tranche_files += [f"{tree}/{root_path}/{mid_path}/{file_name}"]
        print(f"Took me {time.time() - basetime} to get the file list together")
        basetime = time.time()
        session = pyspark.sql.SparkSession.builder.appName("molecule-loader").getOrCreate()
        udf_features = udf(self.features, StringType())

        with concurrent.futures.ThreadPoolExecutor(max_workers=7) as executor:
            futures = {executor.submit(self.load_chem, session, file_name, udf_features)
                       for file_name in _tranche_files}
            for future in concurrent.futures.as_completed(futures):
                _tranche_data = future.result()
                if self.molecule_df:
                    self.molecule_df = self.molecule_df.unionByName(_tranche_data)
                else:
                    self.molecule_df = _tranche_data
        print(f"Took me {time.time() - basetime} to get the DFs loaded")
        basetime = time.time()
        print(self.molecule_df.count())
        print(f"That took {time.time() - basetime}")

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

    @staticmethod
    def load_chem(session, filename, udf_features):
        """
        Loads the chemical list from filename and returns the appropriate
        dictionary.  Best to keep under the memory limit.
        :param session: Spark session context
        :type session: SparkSession
        :param filename: file name to read
        :type filename: str
        :param udf_features: UDF for mapping features
        :type udf_features: pyspark.sql udf
        :return: list of chemicals as defined above
        :rtype: dict
        """
        # session = pyspark.sql.SparkSession.builder.appName(filename).getOrCreate()
        _result = session.read.option("sep", "\t").csv(filename)
        return _result.withColumn("feature_map", udf_features("features"))
        # session.stop()

    def path_match(self, filename, expression):
        """
        Checks whether the last element of the filename matches
        the expression.
        :param filename: file or path name, slash-delimited
        :type filename: str
        :param expression: compiled expression that might match
        :type expression: re.compile
        :return: the filename matches the expression
        :rtype: bool
        """
        return expression.match(filename.split('/')[-1])


if __name__ == "__main__":
    logging.basicConfig(level=logging.INFO)
    test = ChemLoader("../../../files.docking.org/3D/")
