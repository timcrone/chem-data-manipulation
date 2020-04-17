"""
SPARK VERSION

This calculates the Tanimoto similarity of a known
compound with the rest of the ZINC database.

Loads the existing data set with QED and feature mapping
already completed.

After import data stores in a spark RDD, each containing:
- ZINC ID
- SMILES code
- tranche
- net charge
- pH
- molecular weight
- logP
- reactive
- feature map
- QED score

Features map to one or several of:
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

import argparse
import logging
import os
import pyspark
from pyspark.sql.functions import udf, lit
from pyspark.sql.types import FloatType
from rdkit import Chem
from chemfilter import ChemFilter
from chemloader import ChemLoader
from chemtools import ChemTools
from features import Features


class ChemSimilarity(ChemFilter):
    """
    This class pulls the overall molecule data, now that it's pretty.
    """
    def __init__(self, load_path, base_molecule):
        """
        Loads the preprocessed data.

        :param load_path: path containing CSV source data
        :type load_path: str
        :param base_molecule: known compound for comparison
        :type base_molecule: rdkit.molecule
        """
        ChemFilter.__init__(self, load_path)
        self.base_smiles = base_molecule
        mol = Chem.MolFromSmiles(base_molecule)
        self.base_fp = ChemTools.fingerprint(mol)

    def similarity(self):
        """
        Adds a similarity column to the molecules.

        Notably there is not a single similarity, each
        set of similarities will be distinct.

        :return: DataFrame with new column
        :rtype: pyspark.dataframe
        """

        udf_sim = udf(self._similarity, FloatType())
        self.molecule_df = (self.molecule_df
                            .withColumn('similarity_basis', lit(self.base_smiles)))
        self.molecule_df = (self.molecule_df
                            .withColumn('similarity',
                                        udf_sim(self.molecule_df['similarity_basis'],
                                                self.molecule_df['smiles'])
                                        )
                            )
        return self.molecule_df

    @staticmethod
    def _similarity(base_fp, new_smiles):
        """
        Helper for UDF to get the Tanimoto similarity

        :param base_fp: known base fingerprint
        :type base_fp: rdkit.bitvector
        :param new_smiles: incoming SMILES string
        :type new_smiles: str
        :return: Tanimoto similarity
        :rtype: float
        """
        return ChemTools.similarity(Chem.MolFromSmiles(base_fp),
                                    Chem.MolFromSmiles(new_smiles))


if __name__ == "__main__":
    # pylint: disable=invalid-name
    parser = argparse.ArgumentParser("Load molecules and calculate similarity")
    parser.add_argument('smiles', type=str, help="Known SMILES string")
    parser.add_argument('source', metavar="source", type=str,
                        help="Source path")
    parser.add_argument('dest', metavar="dest", type=str,
                        help="Destination base filename")
    args = parser.parse_args()
    logging.basicConfig(level=logging.INFO)
    mols = ChemSimilarity(args.source, args.smiles)
    mols.load()
    mols.similarity()
    mols.write_group(f'similarity', args.dest, "_sim")
    # filtered = mol_class.select("SELECT * FROM dfTable "
    #                             "WHERE qed BETWEEN 0.39999 AND 0.4")
    # filtered.show()
    # for x in range(0, delta - step, step):
    #     fx = float(x)
    #     filtered = mol_class.select(f"SELECT * FROM dfTable "
    #                                 f"WHERE qed BETWEEN "
    #                                 f"{fx / fd} "
    #                                 f"AND {(fx + fs) / fd}")
    #     (filtered.repartition(1).coalesce(1)
    #      .write.csv(f"{args.dest}_qed_{x}", sep="\t",
    #      header=True))
    # stats = ChemStats(mols)
    # print(stats.count())
    # print(stats.describe().show())
    # print(stats.pretty_features().show(n=10000, truncate=False))
