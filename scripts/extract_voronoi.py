import argparse
import os

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="Input File", type=str, default="grids")
    args = parser.parse_args()
    inputFile = args.input
    if not os.path.isfile(inputFile):
        print("Input File does not exist")
        exit()
    
    data = open(inputFile, 'r').read().split("\n")
    data.pop(-1)

    points = int(data[0])

    voronoi_data = {}

    for idx, pt in enumerate(data[1:]):
        pt = " ".join(pt.split()).split(" ")
        voronoi_data[idx + 1] = float(pt[7])

    with open("voronoi", "w+") as the_file:
        for pt in voronoi_data.keys():
            the_file.write("{0}\n".format(voronoi_data[pt]))
        

if __name__ == "__main__":
    main()