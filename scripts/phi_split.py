import os
import argparse
import glob
import shutil

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--folder", help="Grid File Folder", type=str, default="grids")
    parser.add_argument("-p", "--phi_file", help="Phi Vector File", type=str, default="grids")
    parser.add_argument("-b", "--bounds", help="Max and Min bounds", type=str, default="false")
    args = parser.parse_args()
    gridPath = args.folder
    phiFile = args.phi_file
    bounds = str2bool(args.bounds)
    if not os.path.isdir(gridPath):
        print("Grid Folder does not exist")
        exit()

    if not os.path.isfile(phiFile):
        print("Phi Vector File does not exist")
        exit()

    phiValueArray = [None]
    phiValues = open(phiFile).read().split("\n")
    
    for itm in phiValues[:-1]:
        itm = " ".join(itm.split()).split(" ")
        phiValueArray.append(tuple(map(float,itm)))

    try:
        shutil.rmtree("phi_values_split")
    except FileNotFoundError:
        pass

    os.mkdir("phi_values_split")
    
    gridFiles = glob.glob(os.path.join(gridPath, "*"))
    for gridFile in gridFiles:
        phiFileName = "{}{}".format(os.path.basename(gridFile).replace("partGrid","phi-"), ".dat")
        with open(os.path.join("phi_values_split", phiFileName), "w+") as the_file:
            data = open(gridFile).read().split("\n")
            readLimit = 0
            for idx, itm in enumerate(data[:-1]):
                if idx == 0:
                    readLimit = int(itm.split(" ")[2])
                if idx > 0 and idx <= readLimit:
                    pointIdx = int(itm.split(" ")[0])
                    phiSelectArray = phiValueArray[pointIdx]
                    the_file.write("{} {}\n".format(pointIdx, " ".join(map(str, phiSelectArray))))                 

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

if __name__ == "__main__":
    main()