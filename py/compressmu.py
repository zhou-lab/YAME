import sys
import os
import gzip
import h5py
import numpy as np

def opengz(fn,m='r'):

    if fn == "-":
        fh = sys.stdin
    elif fn.endswith('.gz'):
        import gzip
        fh = gzip.open(fn,m)
    else:
        fh = open(fn,m)

    return fh

def main():
    filepath = sys.argv[1]
    if filepath != "-" and not os.path.isfile(filepath):
        print("File path {} does not exist. Exiting...".format(filepath))
        sys.exit()
  
    MU = []
    with opengz(filepath) as fh:
        for i, line in enumerate(fh):
            fields = line.strip().split("\t")
            MU.append((int(fields[0]), int(fields[1])))

    dt = np.dtype([("M", np.uint32), ("U", np.uint32)])
    with h5py.File(sys.argv[2],"w") as f:
        f.create_dataset("MU/%s" % sys.argv[3], data=MU, dtype=dt,
                         compression="gzip", compression_opts=9)

if __name__ == "__main__":
    main()
