import sys
import os
import gzip
import h5py
import numpy as np

def main():
    filepath = sys.argv[1]
    if not os.path.isfile(filepath):
        print("File path {} does not exist. Exiting...".format(filepath))
        sys.exit()
  
    chrmset = []
    chrm = []
    pos = []
    with gzip.open(filepath,"rt") as fh:
        for i, line in enumerate(fh):
            fields = line.strip().split("\t")
            if fields[0] not in chrmset:
                chrmset.append(fields[0])
            chrm.append(chrmset.index(fields[0]))
            pos.append(int(fields[1])+1)

    with h5py.File(sys.argv[2],"w") as f:
        f.create_dataset("rowData/chrmset", data=chrmset, compression="gzip")
        f.create_dataset("rowData/chrm", data=chrm, dtype=np.uint8, compression="gzip")
        f.create_dataset("rowData/pos", data=pos, dtype=np.uint64, compression="gzip")

if __name__ == "__main__":
    main()
