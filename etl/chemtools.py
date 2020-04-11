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

class ChemTools:
    @staticmethod
    def similarity(a, b):
        """
        Calculate Tanimoto similarity between two molecules

        :param a: first molecule
        :type a: rdkit.molecule
        :param b: second molecule
        :type b: rdkit.molecule
        :return: Tanimoto similarity
        :rtype: float
        """
        fp1 = AllChem.GetMorganFingerprintAsBitVect(a, 2,
                                                    nBits=2048,
                                                    useChirality=False)
        fp2 = AllChem.GetMorganFingerprintAsBitVect(b, 2,
                                                    nBits=2048,
                                                    useChirality=False)
        return DataStructs.TanimotoSimilarity(fp1, fp2) 


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

