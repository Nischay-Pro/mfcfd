import os
import argparse
import glob
import shutil

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--folder", help="Phi Split Folder", type=str, default="grids")
    parser.add_argument("-n", "--normalize", help="Normalize", type=str, default="false")
    parser.add_argument("-b", "--bounds", help="Max and Min bounds", type=str, default="false")
    args = parser.parse_args()
    phiPath = args.folder
    normalize = str2bool(args.normalize)
    bounds = str2bool(args.bounds)
    if not os.path.isdir(phiPath):
        print("Phi Split Folder does not exist")
        exit()

    mergedPhi = {}
    
    phiFiles = glob.glob(os.path.join(phiPath, "phid-*"))
    for phiFile in phiFiles:
        data = open(phiFile).read().split("\n")
        for itm in data[1:-1]:
            itm = " ".join(itm.split()).split(" ")
            mergedPhi[int(itm[0])] = tuple(map(float,itm[1:]))
    mergedPhi = dict(sorted(mergedPhi.items()))
    maxPhi = maxDictValue(mergedPhi)
    if not bounds:
        with open("phi_vector_merge.dat", "w+") as the_file:
            for key in mergedPhi.keys():
                if normalize:
                    the_file.write("{}\n".format(" ".join(map(str, normalizeList(mergedPhi[key], maxPhi)))))
                else:
                    the_file.write("{}\n".format(" ".join(map(str, mergedPhi[key]))))
    else:
        minPhi = minDictValue(mergedPhi)
        print("Max Value: {}, Min Value: {}".format(maxPhi, minPhi))


def str2bool(v):
    if isinstance(v, bool):
       return v
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

def maxDictValue(data):
    return max(abs(i) for v in data.values() for i in v)

def minDictValue(data):
    return min(abs(i) for v in data.values() for i in v)

def normalizeList(data, value):
    return [(itm / value) for itm in data]


if __name__ == "__main__":
    main()