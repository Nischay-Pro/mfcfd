import argparse
import os
import numpy as np
import math

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="Input File (Output.dat)", type=str, default="dummy")
    parser.add_argument("-g", "--grid", help="Grid File", type=str, default="dummy")
    parser.add_argument("-m", "--mode", help="Plot Grid Mode", type=int, default=0)
    parser.add_argument("-v", "--voronoi", help="Voronoi Area File", type=str, default="dummy")
    args = parser.parse_args()
    inputFile = args.input
    gridFile = args.grid
    mode = args.mode
    voronoiFile = args.voronoi
    if not os.path.isfile(inputFile):
        print("Output file not found or invalid Output file")
        exit()

    if not os.path.isfile(gridFile):
        print("Grid file not found or invalid Grid file")
        exit()

    if not os.path.isfile(voronoiFile):
        print("Voronoi file not found or invalid Voronoi file")
        exit()

    # data = " ".join(open(inputfile).read().split())
    data = open(inputFile).read().split("\n")
    data.pop(-1)
    data.pop(0)
    data2 = open(gridFile).read().split("\n")
    data2.pop(0)
    data2.pop(-1)
    data3 = open(voronoiFile).read().split("\n")
    data3.pop(-1)
    pointData = {}
    for i in range(1, len(data) + 1):
        temp = " ".join(data[i - 1].split()).split(" ")
        temp2 = " ".join(data2[i - 1].split()).split(" ")
        num_nbhs = int(temp2[7])
        nbhs = tuple(map(int, temp2[8:]))
        assert num_nbhs == len(nbhs)
        pointData[i] = {"x": float(temp2[0]), "y" : float(temp2[1]), "rho": float(temp[4]), "u1": float(temp[5]), "u2": float(temp[6]), "p": float(temp[7]), "nbh": num_nbhs, "conn": nbhs, "voronoi_area": float(data3[i - 1])}
    
    baroclinicVort(pointData, 0)

    for index in pointData.keys():

        rhox = pointData[index]["drho"][0][0]
        rhoy = pointData[index]["drho"][1][0]

        px = pointData[index]["dp"][0][0]
        py = pointData[index]["dp"][1][0]

        ux = pointData[index]["du"][0][0]
        uy = pointData[index]["du"][1][0]

        vx = pointData[index]["dv"][0][0]
        vy = pointData[index]["dv"][1][0]

        bt = (rhox * py - rhoy * px) * (1 / (pointData[index]["rho"] ** 2))
        vor = ux * vy - uy * vx

        pointData[index].update({"bt": bt, "vor": vor})

    if mode == 0:

        with open("special_function_tec.dat", "w+") as the_file:
            the_file.write('TITLE="Special functions"\n')
            the_file.write('VARIABLES = "X","Y","BT","Vorticity", "Vorticity Square", "Vorticity Square x Voronoi Area"\n')
            the_file.write('Zone I = 321 J = 120 F=POINT\n')
            cycle = 1
            for index in pointData.keys():
                if index % 320 == 0:
                    the_file.write("{} {} {} {} {} {}\n".format(pointData[index]["x"], pointData[index]["y"], pointData[index]["bt"], pointData[index]["vor"], pointData[index]["vor"] ** 2, (pointData[index]["vor"] ** 2) * pointData[index]["voronoi_area"]))
                    the_file.write("{} {} {} {} {} {}\n".format(pointData[cycle]["x"], pointData[cycle]["y"], pointData[cycle]["bt"], pointData[cycle]["vor"], pointData[cycle]["vor"] ** 2, (pointData[cycle]["vor"] ** 2) * pointData[cycle]["voronoi_area"]))
                    cycle = index + 1
                elif index == len(pointData):
                    the_file.write("{} {} {} {} {} {}\n".format(pointData[cycle]["x"], pointData[cycle]["y"], pointData[cycle]["bt"], pointData[cycle]["vor"], pointData[cycle]["vor"] ** 2, (pointData[cycle]["vor"] ** 2) * pointData[cycle]["voronoi_area"]))
                    the_file.write("{} {} {} {} {} {}\n".format(pointData[index]["x"], pointData[index]["y"], pointData[index]["bt"], pointData[index]["vor"], pointData[index]["vor"] ** 2, (pointData[index]["vor"] ** 2) * pointData[index]["voronoi_area"]))
                else:
                    the_file.write("{} {} {} {} {} {}\n".format(pointData[index]["x"], pointData[index]["y"], pointData[index]["bt"], pointData[index]["vor"], pointData[index]["vor"] ** 2, (pointData[index]["vor"] ** 2) * pointData[index]["voronoi_area"]))
                
    elif mode == 1:
        with open("special_function_tec.dat", "w+") as the_file:
            for index in pointData.keys():
                the_file.write("{} {} {} {} {} {}\n".format(pointData[index]["x"], pointData[index]["y"], pointData[index]["bt"], pointData[index]["vor"], pointData[index]["vor"] ** 2, (pointData[index]["vor"] ** 2) * pointData[index]["voronoi_area"]))


    elif mode == 2:
        None
        # for j in range(1, 120):
        #     bti = None
        #     vori = None
        #     for i in range(1, 320):
        #         index = i*j

        #         rhox = pointData[index]["drho"][0][0]
        #         rhoy = pointData[index]["drho"][1][0]

        #         px = pointData[index]["dp"][0][0]
        #         py = pointData[index]["dp"][1][0]

        #         ux = pointData[index]["du"][0][0]
        #         uy = pointData[index]["du"][1][0]

        #         vx = pointData[index]["dv"][0][0]
        #         vy = pointData[index]["dv"][1][0]

        #         bt = rhox * py - rhoy * px
        #         vor = ux * vy - uy * vx

        #         if j == 1:
        #             bti = bt
        #             vori = vor

        #         the_file.write("{} {} {} {}\n".format(pointData[index]["x"], pointData[index]["y"], bt, vor))
        #     the_file.write("{} {} {} {}\n".format(pointData[j]["x"], pointData[j]["y"], bti, vori))


def baroclinicVort(pointData, power):
    for i in pointData.keys():
        x_i = pointData[i]["x"]
        y_i = pointData[i]["y"]
        sum_delx_sqr = 0
        sum_dely_sqr = 0
        sum_delx_dely = 0

        sum_delx_delp = np.zeros(1, dtype=np.float64)
        sum_dely_delp = np.zeros(1, dtype=np.float64)

        sum_delx_delr = np.zeros(1, dtype=np.float64)
        sum_dely_delr = np.zeros(1, dtype=np.float64)

        sum_delx_delu = np.zeros(1, dtype=np.float64)
        sum_dely_delu = np.zeros(1, dtype=np.float64)

        sum_delx_delv = np.zeros(1, dtype=np.float64)
        sum_dely_delv = np.zeros(1, dtype=np.float64)

        for conn in pointData[i]["conn"]:
            x_k = pointData[conn]["x"]
            y_k = pointData[conn]["y"]

            delx = x_k - x_i
            dely = y_k - y_i
            dist = math.sqrt(delx*delx + dely*dely)
            weights = dist ** power

            sum_delx_sqr = sum_delx_sqr + ((delx * delx) * weights)
            sum_dely_sqr = sum_dely_sqr + ((dely * dely) * weights)

            sum_delx_dely = sum_delx_dely + ((delx * dely) * weights)

            sum_delx_delp = sum_delx_delp + (weights * delx * (pointData[conn]["p"] - pointData[i]["p"]))
            sum_dely_delp = sum_dely_delp + (weights * dely * (pointData[conn]["p"] - pointData[i]["p"]))

            sum_delx_delr = sum_delx_delr + (weights * delx * (pointData[conn]["rho"] - pointData[i]["rho"]))
            sum_dely_delr = sum_dely_delr + (weights * dely * (pointData[conn]["rho"] - pointData[i]["rho"]))

            sum_delx_delu = sum_delx_delu + (weights * delx * (pointData[conn]["u1"] - pointData[i]["u1"]))
            sum_dely_delu = sum_dely_delu + (weights * dely * (pointData[conn]["u1"] - pointData[i]["u1"]))

            sum_delx_delv = sum_delx_delv + (weights * delx * (pointData[conn]["u2"] - pointData[i]["u2"]))
            sum_dely_delv = sum_dely_delv + (weights * dely * (pointData[conn]["u2"] - pointData[i]["u2"]))

        det = (sum_delx_sqr * sum_dely_sqr) - (sum_delx_dely * sum_delx_dely)
        one_by_det = 1 / det
        
        # For P

        sum_delx_delp1 = sum_delx_delp * sum_dely_sqr
        sum_dely_delp1 = sum_dely_delp * sum_delx_dely

        tempsumpx = one_by_det * (sum_delx_delp1 - sum_dely_delp1)

        sum_dely_delp2 = sum_dely_delp * sum_delx_sqr
        sum_delx_delp2 = sum_delx_delp * sum_delx_dely

        tempsumpy = one_by_det * (sum_dely_delp2 - sum_delx_delp2)

        tempdp = np.array([tempsumpx, tempsumpy], dtype=np.float64)

        # For Rho

        sum_delx_delr1 = sum_delx_delr * sum_dely_sqr
        sum_dely_delr1 = sum_dely_delr * sum_delx_dely

        tempsumrx = one_by_det * (sum_delx_delr1 - sum_dely_delr1)

        sum_dely_delr2 = sum_dely_delr * sum_delx_sqr
        sum_delx_delr2 = sum_delx_delr * sum_delx_dely

        tempsumry = one_by_det * (sum_dely_delr2 - sum_delx_delr2)

        # For U

        sum_delx_delu1 = sum_delx_delu * sum_dely_sqr
        sum_dely_delu1 = sum_dely_delu * sum_delx_dely

        tempsumux = one_by_det * (sum_delx_delu1 - sum_dely_delu1)

        sum_dely_delu2 = sum_dely_delu * sum_delx_sqr
        sum_delx_delu2 = sum_delx_delu * sum_delx_dely

        tempsumuy = one_by_det * (sum_dely_delu2 - sum_delx_delu2)

        # For V

        sum_delx_delv1 = sum_delx_delv * sum_dely_sqr
        sum_dely_delv1 = sum_dely_delv * sum_delx_dely

        tempsumvx = one_by_det * (sum_delx_delv1 - sum_dely_delv1)

        sum_dely_delv2 = sum_dely_delv * sum_delx_sqr
        sum_delx_delv2 = sum_delx_delv * sum_delx_dely

        tempsumvy = one_by_det * (sum_dely_delv2 - sum_delx_delv2)

        tempdp = np.array([tempsumpx, tempsumpy], dtype=np.float64)
        tempdr = np.array([tempsumrx, tempsumry], dtype=np.float64)
        tempdu = np.array([tempsumux, tempsumuy], dtype=np.float64)
        tempdv = np.array([tempsumvx, tempsumvy], dtype=np.float64)

        pointData[i].update({"dp": tempdp, "drho": tempdr, "du": tempdu, "dv": tempdv})

if __name__ == "__main__":
    main()