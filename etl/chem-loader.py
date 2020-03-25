"""
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

import logging
import os
import re
import subprocess
import sys


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

        TODO use subprocess.run to load async

        This loader will not work on Windows due to double slashes.

        :param tree: full path of a subtree
        :type tree: str
        """
        self.tranche = {}
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
                    _tranche = file_name.split('/')[-1].split('.txt')[0]
                    logging.info("Starting %s", _tranche)
                    self.tranche[_tranche] = self.load_chem(f"{tree}/{root_path}/{mid_path}/{file_name}")
                    logging.info("Completed %s", _tranche)

    @staticmethod
    def feature_map(feature_list=None, feature_map=None):
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
        _result = 0
        if feature_list:
            for feature in feature_list.split(','):
                _result |= KNOWN_FEATURES[feature]
        return _result

    def load_chem(self, filename):
        """
        Loads the chemical list from filename and returns the appropriate
        dictionary.  Best to keep under the memory limit.
        :param filename: file name to read
        :type filename: str
        :return: list of chemicals as defined above
        :rtype: dict
        """
        _return = []
        with open(filename, 'r') as chem_file:
            chem_file.readline()  # throw out the header row
            for line in chem_file:
                chunks = line.split("\t")
                _return += [
                    {
                        'zinc': chunks[1],
                        'smile': chunks[0],
                        'tranche': chunks[12],
                        'charge': chunks[5],
                        'ph': chunks[6],
                        'molwt': chunks[7],
                        'logp': chunks[8],
                        'reactive': chunks[10],
                        'features': self.feature_map(chunks[11])

                    }
                ]
        return _return

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
    test = ChemLoader("../../files.docking.org/3D/")
    print(test.tranche)
