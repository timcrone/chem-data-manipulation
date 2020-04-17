"""
This file does some basic manipulations on SMILES codes
using RDKit.  Handy for generating meta-information for
reporting.
"""

import argparse
import logging
import os
import pyspark
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, QED
from numba import jit, cuda


class ChemTools:
    @staticmethod
    @jit(forceobj=True)
    def fingerprint(a):
        """
        Returns the Morgan fingerprint of a molecule

        :param a: molecule to calculate
        :type a: rdkit.molecule
        :return: Morgan fingerprint
        :rtype: rdkit.bitvector
        """
        return AllChem.GetMorganFingerprintAsBitVect(a, 2,
                                                     nBits=2048,
                                                     useChirality=False)

    @classmethod
    @jit(forceobj=True)
    def similarity(cls, a, b):
        """
        Calculate Tanimoto similarity between two molecules

        :param a: first molecule
        :type a: rdkit.molecule
        :param b: second molecule
        :type b: rdkit.molecule
        :return: Tanimoto similarity
        :rtype: float
        """
        fp1 = cls.fingerprint(a)
        return cls.similarity_known(fp1, b)

    @classmethod
    @jit(forceobj=True)
    def similarity_known(cls, fp, b):
        """
        Calculate Tanimoto similarity between two molecules, one
        of which already has a calculated fingerprint.

        :param fp: fingerprint of first molecule
        :type fp: rdkit.bitvector
        :param b: second molecule
        :type b: rdkit.molecule
        :return: Tanimoto similarity
        :rtype: fload
        """
        fp2 = cls.fingerprint(b)
        return DataStructs.TanimotoSimilarity(fp, fp2)


if __name__ == "__main__":
    # pylint: disable=invalid-name
    parser = argparse.ArgumentParser("Maniuplate SMILES codes")
    parser.add_argument('initial', metavar="smiles", type=str,
                        help="SMILES code")
    parser.add_argument('target', metavar="smiles2", type=str,
                        help="Second SMILES code if applicable",
                        default='', nargs='?')
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    mol1 = Chem.MolFromSmiles(args.initial)
    if args.target:
        mol2 = Chem.MolFromSmiles(args.target)
    else:
        mol2 = None
    logging.info(f"QED parameters: {QED.properties(mol1)}")
    logging.info(f"QED value: {QED.qed(mol1)}")
    if mol2:
        logging.info(f"Tanimoto similarity: "
                     f"{ChemTools.similarity(mol1, mol2)}")
