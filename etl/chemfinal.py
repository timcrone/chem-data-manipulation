"""
This loads relatively similar molecules, collates them,
and outputs files in the form
mol1smiles mol2smiles qed1 qed2 similarity

Takes as input a file in the form
zincid smiles: similarity
This should have been stream-filtered by similarity to
the COVID treatments using the scripts in tools/
and thus should fit into memory (hundreds of thousands)

The records are randomly separated and recombined to
generate the output, and similarity and QED scores are
recalculated.
"""

import logging
import random
from rdkit import Chem
from chemtools import ChemTools
from properties import Properties


class ChemFinal:
    """
    Final processing to generate test/train data.
    """
    def __init__(self, source_file):
        """
        Pass the full path of a file (files.docking.org/3D) to
        start the loader.

        :param source_file: full path of a subtree
        :type source_file: str
        """
        self.source = source_file
        self.mols = {}
        self.mol_pairs = []

    def load(self):
        """
        Loads the compound data from the instance source file.

        :return: compound dictionary
        :rtype: dict
        """
        _mols = {}
        with open(self.source, "r") as infile:
            line = infile.readline()
            while line:
                _zinc, _smiles, _ = line.split()
                _smiles = str(_smiles.strip(":"))
                if _zinc.startswith("ZINC") and _smiles:
                    _mols[_zinc] = {
                        'smiles': _smiles,
                        'qed': Properties(_smiles).qed()
                    }
                line = infile.readline()
        self.mols = _mols
        return self.mols

    def pick(self, count=None):
        """
        Makes pairs from the molecules.  Note pairs may not
        be distinct nor unique.

        :param count: number of pairs to create, if None then
                      half the number of self.mols
        :type count: int|None
        :return: molecule pairs as [(smiles1, smiles2, qed1, qed2, sim)]
        :rtype: list
        """
        if not count:
            count = len(self.mols) / 2
        _mol_pairs = []
        while count > 0:
            _new_pair = random.sample(list(self.mols), 2)
            _smile1 = self.mols[_new_pair[0]]['smiles']
            _smile2 = self.mols[_new_pair[1]]['smiles']
            _new_sim = ChemTools.similarity(Chem.MolFromSmiles(_smile1),
                                            Chem.MolFromSmiles(_smile2))
            _new = (_smile1,
                    _smile2,
                    self.mols[_new_pair[0]]['qed'],
                    self.mols[_new_pair[1]]['qed'],
                    _new_sim)
            print(_new[0], _new[1], _new[2], _new[3], _new[4])
            _mol_pairs += [_new]
            count -= 1
        self.mol_pairs = _mol_pairs
        return _mol_pairs


if __name__ == "__main__":
    # pylint: disable=invalid-name
    logging.basicConfig(level=logging.INFO)
    test = ChemFinal("/mnt/ssd-data/mols/covid.list")
    test.load()
    test.pick()
