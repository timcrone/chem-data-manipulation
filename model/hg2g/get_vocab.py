import sys
from hgraph import *
from rdkit import Chem
from multiprocessing import Pool

def process(data):
    vocab = set()
    for line in data:
        s = line.strip("\r\n ")
        hmol = MolGraph(s)
        for node,attr in hmol.mol_tree.nodes(data=True):
            smiles = attr['smiles']
            vocab.add( attr['label'] )
            for i,s in attr['inter_label']:
                vocab.add( (smiles, s) )
    return vocab

if __name__ == "__main__":
    with open("data/qed/covid_mols.txt", 'r') as fh:
        data = [mol for line in fh.readlines() for mol in line.split()[:2]]
    #data = [mol for line in sys.stdin for mol in line.split()[:2]]
    data = list(set(data))
    print(len(data))
    #ncpu = 4
    #batch_size = len(data) // ncpu + 1
    #batches = [data[i : i + batch_size] for i in range(0, len(data), batch_size)]

    #pool = Pool(ncpu)
    #vocab_list = pool.map(process, batches)
    vocab_list = process(data)
    #vocab = [(x,y) for vocab in vocab_list for x,y in vocab]
    vocab = list(set(vocab_list))
    with open("data/qed/covid_vocab.txt", 'w') as fh:
        for x,y in sorted(vocab):
            fh.write(f"{x} {y} \n")
