"""
This file does some basic manipulations on SMILES codes
using RDKit.  Handy for generating meta-information for
reporting.
"""

import argparse
import logging
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, QED


class ChemTools:
    """
    Module containing various methods for manipulating molecules
    """

    @staticmethod
    def fingerprint(mol1):
        """
        Returns the Morgan fingerprint of a molecule

        :param mol1: molecule to calculate
        :type mol1: rdkit.molecule
        :return: Morgan fingerprint
        :rtype: rdkit.bitvector
        """
        return AllChem.GetMorganFingerprintAsBitVect(mol1, 2,
                                                     nBits=2048,
                                                     useChirality=False)

    @classmethod
    def similarity(cls, mol1, mol2):
        """
        Calculate Tanimoto similarity between two molecules

        :param mol1: first molecule
        :type mol1: rdkit.molecule
        :param mol2: second molecule
        :type mol2: rdkit.molecule
        :return: Tanimoto similarity
        :rtype: float
        """
        fprint1 = cls.fingerprint(mol1)
        return cls.similarity_known(fprint1, mol2)

    @classmethod
    def similarity_known(cls, fprint, mol2):
        """
        Calculate Tanimoto similarity between two molecules, one
        of which already has a calculated fingerprint.

        :param fprint: fingerprint of first molecule
        :type fprint: rdkit.bitvector
        :param mol2: second molecule
        :type mol2: rdkit.molecule
        :return: Tanimoto similarity
        :rtype: fload
        """
        fprint2 = cls.fingerprint(mol2)
        return DataStructs.TanimotoSimilarity(fprint, fprint2)


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
    molecule1 = Chem.MolFromSmiles(args.initial)
    if args.target:
        molecule2 = Chem.MolFromSmiles(args.target)
    else:
        molecule2 = None
    logging.info("QED parameters: %s", QED.properties(molecule1))
    logging.info("QED value: %s", QED.qed(molecule1))
    if molecule2:
        logging.info("Tanimoto similarity: %s",
                     ChemTools.similarity(molecule1, molecule2))
