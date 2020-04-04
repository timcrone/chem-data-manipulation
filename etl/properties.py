"""
This library contains tools that help characterize
molecules.
"""

import logging
from rdkit import Chem
from rdkit.Chem import Lipinski
from rdkit.Chem import QED


class Properties:
    """
    Defines methods for working with compound properties.
    """
    def __init__(self, smiles):
        """
        Pass the SMILES string to initialize.

        :param smiles: SMILES-format string for the target
                       compound
        :type smiles: str
        """
        self.smiles = smiles
        self.molecule = Chem.MolFromSmiles(smiles)

    def h_bond_donors(self):
        """
        :return: number of Hydrogen bond donors
        :rtype: int
        """
        return Chem.Lipinski.NumHDonors(self.molecule)

    def h_bond_acceptors(self):
        """
        :return: number of Hydrogen bond acceptors
        :rtype: int
        """
        return Chem.Lipinski.NumHAcceptors(self.molecule)

    def rotatable_bonds(self):
        """
        :return: number of rotatable bonds
        :rtype: int
        """
        return Chem.Lipinski.NumRotatableBonds(self.molecule)

    def qed(self):
        """
        :return: QED using default parameters
        :rtype: double
        """
        return Chem.QED.default(self.molecule)


if __name__ == "__main__":
    # pylint: disable=invalid-name
    logging.basicConfig(level=logging.INFO)
    # Randomly selected ZINC000226860047 from FE/CEHL/FECEHL.txt
    test_string = "O=C1[N-]C(=O)c2cc(S(=O)(=O)[N-]c3ccc(O)c4ccccc34)ccc21"
    test = Properties(test_string)
    assert int(test.qed() * 1000) == 710
    assert test.h_bond_acceptors() == 5
    assert test.h_bond_donors() == 1
    assert test.rotatable_bonds() == 3
