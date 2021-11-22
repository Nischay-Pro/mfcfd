import os
import argparse
import glob
import shutil

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-f", "--folder", help="Grid File Folder", type=str, default="grids")
    parser.add_argument("-v", "--voronoi", help="Voronoi File", type=str, default="grids")
    args = parser.parse_args()
    gridPath = args.folder
    voronoiFile = args.voronoi
    if not os.path.isdir(gridPath):
        print("Grid Folder does not exist")
        exit()

    if not os.path.isfile(voronoiFile):
        print("Voronoi File File does not exist")
        exit()

    voronoiValueArray = [None]
    voronoiValues = open(voronoiFile).read().split("\n")
    
    for itm in voronoiValues[:-1]:
        itm = " ".join(itm.split()).split(" ")
        voronoiValueArray.append(tuple(map(float,itm)))

    try:
        shutil.rmtree("voronoi_values_split")
    except FileNotFoundError:
        pass

    os.mkdir("voronoi_values_split")
    
    gridFiles = glob.glob(os.path.join(gridPath, "*"))
    for gridFile in gridFiles:
        voronoiFileName = "{}{}".format(os.path.basename(gridFile).replace("partGrid","voronoi"), ".dat")
        with open(os.path.join("voronoi_values_split", voronoiFileName), "w+") as the_file:
            data = open(gridFile).read().split("\n")
            readLimit = 0
            for idx, itm in enumerate(data[:-1]):
                if idx == 0:
                    readLimit = int(itm.split(" ")[2])
                if idx > 0 and idx <= readLimit:
                    pointIdx = int(itm.split(" ")[0])
                    voronoiSelectArray = voronoiValueArray[pointIdx]
                    the_file.write("{} {}\n".format(pointIdx, " ".join(map(str, voronoiSelectArray))))                 

if __name__ == "__main__":
    main()