import argparse
import os
import numpy as np
from scipy.spatial import Delaunay
import math
import triangle as tr
import matplotlib.pyplot as plt
from shapely.geometry import Point, Polygon

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("-i", "--input", help="Input File (Output.dat)", type=str, default="dummy")
    parser.add_argument("-g", "--grid", help="Grid File", type=str, default="dummy")
    parser.add_argument("-m", "--mode", help="Plot Grid Mode", type=int, default=0)
    parser.add_argument("-v", "--voronoi", help="Voronoi Area File", type=str, default="dummy")
    parser.add_argument("-b", "--base-phi", help="Base Phi File", type=str, default="dummy")
    parser.add_argument("-o", "--optimal-phi", help="Optimal Phi File", type=str, default="dummy")
    args = parser.parse_args()
    inputFile = args.input
    gridFile = args.grid
    mode = args.mode
    voronoiFile = args.voronoi
    basePhiFile = args.base_phi
    optimalPhiFile = args.optimal_phi
    if not os.path.isfile(inputFile):
        print("Output file not found or invalid Output file")
        exit()

    if not os.path.isfile(gridFile):
        print("Grid file not found or invalid Grid file")
        exit()

    if not os.path.isfile(voronoiFile):
        print("Voronoi file not found or invalid Voronoi file")
        exit()

    if not os.path.isfile(basePhiFile):
        print("Base Phi file not found or invalid Phi file")
        exit()
    
    if not os.path.isfile(optimalPhiFile):
        print("Optimal Phi file not found or invalid Phi file")
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
    data4 = open(basePhiFile).read().split("\n")
    data4.pop(-1)
    data5 = open(optimalPhiFile).read().split("\n")
    data5.pop(-1)
    pointData = {}
    generalData = {"geometry_pts": [], "geometry": {}}
    points = []
    wallPoints = []
    outerPoints = []
    wallSegments = []
    wallSegmentsCords = []
    outerSegments = []
    for i in range(1, len(data) + 1):
        temp = " ".join(data[i - 1].split()).split(" ")
        if len(temp) == 11:
            temp.pop(2)
            temp.pop(2)
        elif len(temp) == 10:
            temp.pop(2)
        temp2 = " ".join(data2[i - 1].split()).split(" ")
        temp4 = " ".join(data4[i - 1].split()).split(" ")
        temp5 = " ".join(data5[i - 1].split()).split(" ")
        if int(temp2[5]) > 0:
            wallPoints.append(i)
        if int(temp2[4]) == 2:
            outerPoints.append(i)
        num_nbhs = int(temp2[7])
        nbhs = tuple(map(int, temp2[8:]))
        assert num_nbhs == len(nbhs)
        pointData[i] = {"x": float(temp2[0]), "y": float(temp2[1]), "left": int(temp2[2]), "right": int(temp2[3]), "type": int(temp2[4]), "geometry": int(temp2[5]), "rho": float(temp[2]), "u1": float(temp[3]), "u2": float(temp[4]), "p": float(temp[5]), "nbh": num_nbhs, "conn": nbhs, "voronoi_area": float(data3[i - 1])}
        pointData[i].update({"base_phi": { "phi11": float(temp4[0]), "phi12": float(temp4[1]), "phi13": float(temp4[2]), "phi14": float(temp4[3]), "phi21": float(temp4[4]), "phi22": float(temp4[5]), "phi23": float(temp4[6]), "phi24": float(temp4[7])}})
        pointData[i].update({"optimal_phi": { "phi11": float(temp5[0]), "phi12": float(temp5[1]), "phi13": float(temp5[2]), "phi14": float(temp5[3]), "phi21": float(temp5[4]), "phi22": float(temp5[5]), "phi23": float(temp5[6]), "phi24": float(temp5[7])}})
        points.append([float(temp[0]), float(temp[1])])

    calculateBaroclinicVort(pointData)
    calculateCirculation(pointData, generalData)

    tempSegments = []
    tempCordSegments = []
    baseGeometry = 1
    for wallpt in wallPoints:
        currGeometry = pointData[wallpt]["geometry"]
        if currGeometry != baseGeometry:
            wallSegments.append(tempSegments)
            wallSegmentsCords.append(tempCordSegments)
            tempSegments = []
            tempCordSegments = []
            baseGeometry = currGeometry
        leftIdx = pointData[wallpt]["left"] - 1
        rightIdx = pointData[wallpt]["right"] - 1
        tempSegments.append([wallpt - 1, rightIdx])
        tempCordSegments.append([pointData[wallpt]["x"], pointData[wallpt]["y"]])
    wallSegments.append(tempSegments)
    wallSegmentsCords.append(tempCordSegments)

    for outerpt in outerPoints:
        leftIdx = pointData[outerpt]["left"] - 1
        rightIdx = pointData[outerpt]["right"] - 1
        outerSegments.append([outerpt - 1, rightIdx])

    segments = []
    for segment in wallSegments:
        segments += segment
    
    segments += outerSegments

    # holes = [(0.367, 0), (1.14, -0.06)]
    holes = [(0.31, -0.012)]
    # for wall in wallSegmentsCords:
    #     holes.append(random_points_within(wall, 1))

    # segments = np.vstack((wallSegments, outerSegments))
    points = np.array(points)
    tri_data = {"vertices": points, "segments": segments, "holes": holes}

    # segment_idx = 1
    # with open("triangle.poly", "w+") as f:
    #     f.write("# Generated\n")
    #     f.write("{} 2 0 1\n".format(len(pointData)))
    #     f.write("# Shape\n")
    #     for wallpt in wallPoints:
    #         f.write("{} {} {} {}\n".format(wallpt, pointData[wallpt]["x"], pointData[wallpt]["y"], 1))

    #     for pt in pointData.keys():
    #         if pt not in wallPoints:
    #             f.write("{} {} {} {}\n".format(pt, pointData[pt]["x"], pointData[pt]["y"], pointData[pt]["type"] + 100))
        
    #     f.write("# Line Segments\n")
    #     f.write("{} 1\n".format(len(segments)))
    #     for wallpt in wallPoints:
    #         f.write("{} {} {} {}\n".format(segment_idx, wallpt, pointData[wallpt]["right"], 1))
    #         segment_idx += 1
        
    #     for outerpt in outerPoints:
    #         f.write("{} {} {} {}\n".format(segment_idx, outerpt, pointData[outerpt]["right"], 102))
    #         segment_idx += 1
        
    #     f.write("1\n")
    #     f.write("# Holes\n")
    #     f.write("1 0.08 0.01\n")

    triangles = tr.triangulate(tri_data, "p")

    triangle_set = triangles["triangles"]
    print(len(points))
    print(len(triangles["triangles"]))
    print(len(triangles["vertices"]))
#    tr.compare(plt, tri_data, triangles)
#    plt.show()
#    exit()
    print(len(triangles["vertices"]))

    with open("special_function_tec.dat", "w+") as the_file:
        the_file.write('TITLE="Special functions"\n')
        the_file.write('VARIABLES = "X", "Y", "rho", "u1", "u2", "pressure", "mach", "baroclinic torque", "vorticity"\n')
        the_file.write('ZONE T="Triangulation", N={0}, E={1}, DATAPACKING=POINT, ZONETYPE=FETRIANGLE\n'.format(len(pointData), len(triangle_set)))
        for i in pointData.keys():
            p = pointData[i]
            the_file.write("{0} {1} {2} {3} {4} {5} {6} {7} {8}\n".format(p["x"], p["y"], p["rho"], p["u1"], p["u2"], p["p"], p["mach"], p["baroclinic"], p["vorticity"]))
        for i in triangle_set:
            the_file.write("{0} {1} {2}\n".format(i[0] + 1, i[1] + 1, i[2] + 1))
    
    with open("phi_vector_1_tec.dat", "w+") as the_file:
        the_file.write('TITLE="Phi vector 1"\n')
        the_file.write('VARIABLES = "X", "Y", "phi 11", "phi 12", "phi 13", "phi 14", "delta phi 11", "delta phi 12", "delta phi 13", "delta phi 14"\n')
        the_file.write('ZONE T="Triangulation", N={0}, E={1}, DATAPACKING=POINT, ZONETYPE=FETRIANGLE\n'.format(len(pointData), len(triangle_set)))
        for i in pointData.keys():
            p = pointData[i]
            the_file.write("{0} {1} {2} {3} {4} {5} {6} {7} {8} {9}\n".format(p["x"], p["y"], p["base_phi"]["phi11"], p["base_phi"]["phi12"], p["base_phi"]["phi13"], p["base_phi"]["phi14"], p["optimal_phi"]["phi11"] - p["base_phi"]["phi11"], p["optimal_phi"]["phi12"] - p["base_phi"]["phi12"], p["optimal_phi"]["phi13"] - p["base_phi"]["phi13"], p["optimal_phi"]["phi14"] - p["base_phi"]["phi14"]))
        for i in triangle_set:
            the_file.write("{0} {1} {2}\n".format(i[0] + 1, i[1] + 1, i[2] + 1))
    
    with open("phi_vector_2_tec.dat", "w+") as the_file:
        the_file.write('TITLE="Phi vector 2"\n')
        the_file.write('VARIABLES = "X", "Y", "phi 21", "phi 22", "phi 23", "phi 24", "delta phi 21", "delta phi 22", "delta phi 23", "delta phi 24"\n')
        the_file.write('ZONE T="Triangulation", N={0}, E={1}, DATAPACKING=POINT, ZONETYPE=FETRIANGLE\n'.format(len(pointData), len(triangle_set)))
        for i in pointData.keys():
            p = pointData[i]
            the_file.write("{0} {1} {2} {3} {4} {5} {6} {7} {8} {9}\n".format(p["x"], p["y"], p["base_phi"]["phi21"], p["base_phi"]["phi22"], p["base_phi"]["phi23"], p["base_phi"]["phi24"], p["optimal_phi"]["phi21"] - p["base_phi"]["phi21"], p["optimal_phi"]["phi22"] - p["base_phi"]["phi22"], p["optimal_phi"]["phi23"] - p["base_phi"]["phi23"], p["optimal_phi"]["phi24"] - p["base_phi"]["phi24"]))
        for i in triangle_set:
            the_file.write("{0} {1} {2}\n".format(i[0] + 1, i[1] + 1, i[2] + 1))

    with open("enstrophy.dat", "w+") as the_file:
        the_file.write('TITLE="Enstrophy"\n')
        the_file.write('VARIABLES = "X", "Y", "enstrophy"\n')
        the_file.write('ZONE T="Triangulation", N={0}, E={1}, DATAPACKING=POINT, ZONETYPE=FETRIANGLE\n'.format(len(pointData), len(triangle_set)))
        for i in pointData.keys():
            p = pointData[i]
            the_file.write("{0} {1} {2}\n".format(p["x"], p["y"], p["enstrophy"]))
        for i in triangle_set:
            the_file.write("{0} {1} {2}\n".format(i[0] + 1, i[1] + 1, i[2] + 1))

    with open("circulation", "w+") as the_file:
        for ptidx in generalData["geometry_pts"]:
            pt = pointData[ptidx]
            the_file.write("{0} {1} {2}\n".format(pt["geometry"], pt["x"], pt["circulation"]))
    
    for idx, itm in enumerate(generalData["maxCirculation"]):
        print("Total Circulation {0}: {1}".format(idx, itm))

def calculateCirculation(pointData, generalData):
    maxGeo = 0
    # max_circulation = [0, 0]
    max_circulation = [0]
    for ptidx in pointData.keys():
        pt = pointData[ptidx]
        if pt["geometry"] > 0:
            generalData["geometry_pts"].append(ptidx)
            maxGeo = max(maxGeo, pt["geometry"])
            mx = pt["x"]
            my = pt["y"]

            leftidx = pt["left"]
            rightidx = pt["right"]

            lx = pointData[leftidx]["x"]
            ly = pointData[leftidx]["y"]

            rx = pointData[rightidx]["x"]
            ry = pointData[rightidx]["y"]

            dx1 = mx - lx
            dy1 = my - ly

            dx2 = rx - mx
            dy2 = ry - my

            dx = (dx1 + dx2) / 2
            dy = (dy1 + dy2) / 2

            u1 = pt["u1"]
            u2 = pt["u2"]

            pointData[ptidx].update({"circulation": (u1 * dx + u2 * dy)})
            max_circulation[pt["geometry"] - 1] += (u1 * dx + u2 * dy)

            try:
                generalData["geometry"][pt["geometry"]] += pt["circulation"]
            except KeyError:
                generalData["geometry"][pt["geometry"]] = pt["circulation"]

    generalData.update({"maxGeo": maxGeo, "maxCirculation": max_circulation})


def calculateBaroclinicVort(pointData):
    total_enstrophy = 0
    for i in pointData.keys():
        x_i = pointData[i]["x"]
        y_i = pointData[i]["y"]

        sum_delx_sqr = 0
        sum_dely_sqr = 0
        sum_delx_dely = 0

        sum_delx_delu1 = 0
        sum_delx_delu2 = 0
        sum_dely_delu1 = 0
        sum_dely_delu2 = 0

        sum_delx_delP = 0
        sum_delx_delrho = 0
        sum_dely_delP = 0
        sum_dely_delrho = 0

        for conn in pointData[i]["conn"]:
            x_k = pointData[conn]["x"]
            y_k = pointData[conn]["y"]

            delx = x_k - x_i
            dely = y_k - y_i

            sum_delx_sqr += delx ** 2
            sum_dely_sqr += dely ** 2

            sum_delx_dely += delx * dely

            sum_delx_delu1 += delx * (pointData[conn]["u1"] - pointData[i]["u1"])
            sum_delx_delu2 += delx * (pointData[conn]["u2"] - pointData[i]["u2"])

            sum_dely_delu1 += dely * (pointData[conn]["u1"] - pointData[i]["u1"])
            sum_dely_delu2 += dely * (pointData[conn]["u2"] - pointData[i]["u2"])

            sum_delx_delP += delx * (pointData[conn]["p"] - pointData[i]["p"])
            sum_delx_delrho += delx * (pointData[conn]["rho"] - pointData[i]["rho"])

            sum_dely_delP += dely * (pointData[conn]["p"] - pointData[i]["p"])
            sum_dely_delrho += dely * (pointData[conn]["rho"] - pointData[i]["rho"])

        det = sum_delx_sqr * sum_dely_sqr - sum_delx_dely ** 2
        one_by_det = 1 / det

        du2_dx = (sum_delx_delu2 * sum_dely_sqr - sum_dely_delu2 * sum_delx_dely) * one_by_det
        du1_dy = (sum_dely_delu1 * sum_delx_sqr - sum_delx_delu1 * sum_delx_dely) * one_by_det

        du2_dy = (sum_dely_delu2 * sum_delx_sqr - sum_delx_delu2 * sum_delx_dely) * one_by_det
        du1_dx = (sum_delx_delu1 * sum_dely_sqr - sum_dely_delu1 * sum_delx_dely) * one_by_det

        dP_dx = (sum_delx_delP * sum_dely_sqr - sum_dely_delP * sum_delx_dely) * one_by_det
        drho_dy = (sum_dely_delrho * sum_delx_sqr - sum_delx_delrho * sum_delx_dely) * one_by_det

        dP_dy = (sum_dely_delP * sum_delx_sqr - sum_delx_delP * sum_delx_dely) * one_by_det
        drho_dx = (sum_delx_delrho * sum_dely_sqr - sum_dely_delrho * sum_delx_dely) * one_by_det

        vorticity = du2_dx - du1_dy

        baroclinicTorque = (drho_dx * dP_dy - dP_dx * drho_dy) * (1 / pointData[i]["rho"] ** 2)
        
        enstrophy = ((du1_dx * du1_dx) + (du2_dx * du2_dx) + (du1_dy * du1_dy) + (du2_dy * du2_dy)) * pointData[i]["voronoi_area"]

        total_enstrophy += enstrophy

        gamma = 1.4

        sos = math.sqrt(gamma * pointData[i]["p"] / pointData[i]["rho"])
        vel_mag = math.sqrt(pointData[i]["u1"] ** 2 + pointData[i]["u2"] ** 2)
        Mach = vel_mag / sos

        pointData[i].update({"vorticity": vorticity, "baroclinic": baroclinicTorque, "mach": Mach, "enstrophy": enstrophy})
        

    print("Total Enstrophy is {0}".format(total_enstrophy))

def random_points_within(poly, num_points):
    poly = Polygon(poly)
    min_x, min_y, max_x, max_y = poly.bounds

    points = []
    while len(points) < num_points:
        random_point = Point([np.random.uniform(min_x, max_x), np.random.uniform(min_y, max_y)])
        if (random_point.within(poly)):
            points.append(random_point)

    return points

if __name__ == "__main__":
    main()
