import os
import argparse

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--base", help="Phi Base", type=str, default="grids")
    args = parser.parse_args()

    phibase = args.base

    if not os.path.isfile(phibase):
        print("Phi Base File does not exist")
        exit()

    phiValueArray = []
    phiValues = open(phibase).read().split("\n")
    
    for itm in phiValues[:-1]:
        itm = " ".join(itm.split()).split(" ")
        phiValueArray.append(list(map(float,itm)))

    with open("phi_vector_dup.dat", "w+") as the_file:
        for itm in phiValueArray:
            data = " ".join(map(str,itm))
            the_file.write("{} {}\n".format(data, data))

if __name__ == "__main__":
    main()

if __name__ == "__main__":
    main()