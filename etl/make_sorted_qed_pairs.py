import pandas as pd
import numpy as np
from collections import Counter
import sys

import rdkit
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem
from rdkit.Chem import Descriptors
import rdkit.Chem.QED as QED
import networkx as nx

# Generates sorted pairs training sets
# In: python etl/make_sorted_qed_pairs.py < data/pairs
# Out: data/covid_train_pairs.txt
# Out: data/covid_valid.txt
# Out: data/covid_test.txt
# Out: data/covid_mols.txt 

VAL_SIZE = 400
TEST_SIZE = 800

def similarity(a, b):
    if a is None or b is None: 
        return 0.0
    amol = Chem.MolFromSmiles(a)
    bmol = Chem.MolFromSmiles(b)
    if amol is None or bmol is None:
        return 0.0

    fp1 = AllChem.GetMorganFingerprintAsBitVect(amol, 2, nBits=2048, useChirality=False)
    fp2 = AllChem.GetMorganFingerprintAsBitVect(bmol, 2, nBits=2048, useChirality=False)
    return DataStructs.TanimotoSimilarity(fp1, fp2) 

def qed(s):
    if s is None: return 0.0
    mol = Chem.MolFromSmiles(s)
    if mol is None: return 0.0
    return QED.qed(mol)
bad_values = 0
sorted_pairs = []
if sys.stdin:
    for line in sys.stdin:
        values = line.split(" ")
        if len(values) == 5:
            x,y,qed1,qed2,sim2D = values
        else:
            bad_values +=1
            if bad_values < 5:
                print(line)
            continue
        if y == "None": y = None
        try: 
            sim2D = float(sim2D)
        except:
            sim2D = similarity(x, y)
        try: 
            qed1 = float(qed1)
        except:
            qed1 = qed(x)
        try: 
            qed2 = float(qed2)
        except:
            qed2 = qed(y)

        if qed1 < qed2:
            sorted_pairs.append((x,y,sim2D,qed1,qed2))
        else:
            sorted_pairs.append((y,x,sim2D,qed2,qed1))

    print(f"{bad_values} lines are bad")
    df = pd.DataFrame(sorted_pairs, columns=['mol1','mol2','sim','qed1','qed2'])
    df = df.drop_duplicates()
    df = df.query('sim < 1.0')
    df['qed_delta'] = df['qed2'] - df['qed1']
    # remove duplicated 
    df.to_pickle("data/df_covid_clean-04192125.pickle")
    print('-'*50)
    print(df['sim'].describe())
    print('-'*50)
    print(df['qed_delta'].describe())
    print('-'*50)
    for i in np.arange(0.1, 0.99, 0.1):
        print(f"{len(df[df['sim']>=i])} training pairs with sim >= {round(i, 2)}")
    print('-'*50)
    # The ideal training set consists of sim >= 0.4 pairs
    train_pairs = df.query("sim > 0.5 and qed_delta > 0.2")
    leftovers = df.query("sim < 0.5")
    valid = leftovers[:VAL_SIZE]
    test = leftovers[VAL_SIZE:VAL_SIZE+TEST_SIZE]

    train_pairs[['mol1','mol2']].to_csv("data/covid_train_pairs.txt", sep=' ', header=False, index=False)

    valid[['mol1']].to_csv("data/covid_valid.txt", sep=' ', header=False, index=False)
    test[['mol1']].to_csv("data/covid_test.txt", sep=' ', header=False, index=False)
    # TODO(vg) May need to handle out of vocab values in train set
    mols = list(train_pairs[['mol1']].values) + list(train_pairs[['mol2']].values) + list(valid[['mol1']].values) + list(test[['mol1']].values)
    pd.DataFrame(mols).drop_duplicates().to_csv("data/covid_mols.txt", sep=' ', header=False, index=False)

