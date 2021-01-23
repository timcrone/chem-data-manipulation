# NOTE THIS DOES NOT WORK
# For some reason there is a flipped bit that gets stuck between
# quick reads and writes; after readAll the read queue reports cleared,
# but the next read is corrupted.  Using the non-file-based version
# this does not appear because QT and Python reinitialize, and we
# were not able to fix it faster than just running the
# stream-based version.
# Including for future update and repair if a faster solution is
# required.

from PyQt5 import QtCore, QtNetwork

import random
from rdkit import Chem
from rdkit.Chem import QED
from gpusim_utils import smiles_to_fingerprint_bin


def parse_args():
    import argparse
    parser = argparse.ArgumentParser(description="Sample GPUSim Server - "
            "run an HTTP server that loads fingerprint data onto GPU and " #noqa
            "responds to queries to find most similar fingperints.") #noqa
    parser.add_argument('dbname', help=".fsim file containing fingerprint "
                        "data to be searched")
    parser.add_argument('dbkey', default="", help="Key for fsim file")
    parser.add_argument('smi_file', default="covid.smi", help="SMILES file")
    return parser.parse_args()


def main():
    args = parse_args()
    app = QtCore.QCoreApplication([])

    socket = QtNetwork.QLocalSocket(app)
    dbcount = 1
    dbname = args.dbname
    dbkey = args.dbkey
    socket.connectToServer('gpusimilarity')

    # precalculate QED scores
    in_smiles = {}
    with open(args.smi_file) as smifile:
        lines = smifile.readlines()
        for line in lines:
            _smi, _id = line.split()
            _mol = Chem.MolFromSmiles(_smi)
            in_smiles[_smi] = {
                    'id': _id,
                    'qed': QED.qed(_mol)
                    }

    #while smiles and smiles.lower() not in ('quit', 'exit'):
    for i_smiles in in_smiles:
        return_count = 1000
        similarity_cutoff = 0.4

        fp_binary, _ = smiles_to_fingerprint_bin(i_smiles)
        fp_qba = QtCore.QByteArray(fp_binary)

        output_qba = QtCore.QByteArray()
        output_qds = QtCore.QDataStream(output_qba, QtCore.QIODevice.WriteOnly)

        output_qds.writeInt(dbcount)
        output_qds.writeString(dbname.encode())
        output_qds.writeString(dbkey.encode())

        request_num = random.randint(0, 2**31)
        output_qds.writeInt(request_num)
        output_qds.writeInt(return_count)
        output_qds.writeFloat(similarity_cutoff)
        output_qds << fp_qba

        socket.write(output_qba)
        socket.flush()
        socket.waitForReadyRead(30000)
        output_qba = socket.readAll()

        smiles = []
        scores = []
        ids = []

        data_reader = QtCore.QDataStream(output_qba)
        returned_request = data_reader.readInt()
        if request_num != returned_request:
            raise RuntimeError("Incorrect result ID returned! %s %s", returned_request, request_num)

        return_count = data_reader.readInt()

        for i in range(return_count):
            smiles.append(data_reader.readString())
        for i in range(return_count):
            ids.append(data_reader.readString())
        for i in range(return_count):
            scores.append(data_reader.readFloat())

        for cid, smi, score in zip(ids, smiles, scores):
            if score < similarity_cutoff:
                continue
            f_smiles = smi.decode("utf-8")
            print("{} {} {} {} {}".format(i_smiles,
                                          f_smiles,
                                          in_smiles[i_smiles],
                                          in_smiles[f_smiles],
                                          score))
        #smiles = input("Smiles: ")


if __name__ == '__main__':
    main()
