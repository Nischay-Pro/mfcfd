import geopandas as gpd
from shapely.ops import voronoi_diagram as svd
import shapely.ops as so
from shapely.geometry import Polygon, MultiPolygon, MultiPoint, Point, LineString
import os
import argparse
import math
import matplotlib.pyplot as plt
import tqdm
from scipy.spatial import Voronoi, voronoi_plot_2d
import matplotlib.pyplot as plt
import numpy as np

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

    globalPoints = []
    wallPoints = []
    outerPoints = []

    min_geometry = 10000000000
    globaldata = {}
    for idx, pt in enumerate(data[1:]):
        pt = " ".join(pt.split())
        tmp = {
            "index": idx + 1,
            "x": float(pt.split(" ")[0]),
            "y": float(pt.split(" ")[1]),
            "left": int(pt.split(" ")[2]),
            "right": int(pt.split(" ")[3]),
            "type": int(pt.split(" ")[4]),
            "geometry": int(pt.split(" ")[5]),
            # "nx": float(pt.split(" ")[6]),
            # "ny": float(pt.split(" ")[7]),
            "min_dist": float(pt.split(" ")[6]),
            "nbhs": tuple(map(int, pt.split(" ")[7:]))
        }
        if tmp["type"] == 0:
            wallPoints.append(idx + 1)
            min_geometry = min(min_geometry, tmp["geometry"])
        elif tmp["type"] == 1 and tmp["geometry"] > 0:
            wallPoints.append(idx + 1)
        elif tmp["type"] == 2:
            outerPoints.append(idx + 1)
        globalPoints.append((tmp["x"], tmp["y"]))
        globaldata[idx + 1] = tmp

    outer_radius = get_outer_radius(globaldata, outerPoints)
    outer_pts, _ = circle(len(outerPoints) * 2,  outer_radius * 1.01)

    for pt in outer_pts:
        globalPoints.append(pt)
    
    assert len(globaldata.keys()) == points

    globalPoints = MultiPoint(globalPoints)

    basePoints = []
    current_geometry = 1

    for pt in wallPoints:
        if globaldata[pt]["geometry"] == current_geometry:
            basePoints.append(pt)
            current_geometry += 1

    wallGeometry = []
    multiWalls = []

    for pt in basePoints:
        tmp = []
        tmp.append(pt)
        current_pt = pt
        while True:
            if current_pt == tmp[0] and len(tmp) > 1:
                tmp.pop(-1)
                break
            else:
                current_pt = globaldata[current_pt]["left"]
                tmp.append(current_pt)
        wallGeometry.append(tmp)
    
    for wallPoints in wallGeometry:
        tmp = []
        for pt in wallPoints:
            tmp.append((globaldata[pt]["x"], globaldata[pt]["y"]))
        multiWalls.append(Polygon(tmp))

    # vor = Voronoi(globalPoints)
    # fig = voronoi_plot_2d(vor)

    # plt.show()
    # exit()
    voronoi = svd(globalPoints)
    
    # assert len(voronoi.geoms) == points

    search_list = list(globaldata.keys())

    # voronoi_list = []

    # for pt in tqdm.tqdm(search_list):
    #     check_pt = Point(globaldata[pt]["x"], globaldata[pt]["y"])
    #     for vor_idx, vor_poly in enumerate(voronoi.geoms):
    #         if vor_idx not in voronoi_list:
    #             if vor_poly.contains(check_pt):
    #                 globaldata[pt]["voronoi"] = vor_poly
    #                 globaldata[pt]["voronoi_area"] = vor_poly.area
    #                 voronoi_list.append(vor_idx)
    #                 search_list.remove(pt)
    #                 break

    for vor_poly in tqdm.tqdm(voronoi.geoms):
        # checks = 0
        for pt in search_list:
            # check_line = LineString([(globaldata[pt]["x"], globaldata[pt]["y"]), (globaldata[pt]["x"] + outer_radius, globaldata[pt]["y"])])
            check_pt = Point(globaldata[pt]["x"], globaldata[pt]["y"])
            if vor_poly.contains(check_pt):
                globaldata[pt]["voronoi"] = vor_poly
                globaldata[pt]["voronoi_area"] = vor_poly.area
                search_list.remove(pt)
                # print(checks)
                # checks = 0
                break

    for wallPt in wallPoints:
        wallCheck = multiWalls[globaldata[wallPt]["geometry"] - 1]
        new_vor = globaldata[wallPt]["voronoi"].difference(wallCheck)
        globaldata[wallPt]["voronoi"] = new_vor
        globaldata[wallPt]["voronoi_area"] = new_vor.area
    
    for outerPt in outerPoints:
        nbhs = globaldata[outerPt]["nbhs"]
        interior_area = 0
        interior_count = 0
        for nbh in nbhs:
            if globaldata[nbh]["type"] == 1:
                interior_count += 1
                interior_area += globaldata[nbh]["voronoi_area"]
        if interior_count > 0:
            globaldata[outerPt]["voronoi"] = None
            globaldata[outerPt]["voronoi_area"] = interior_area / interior_count
    
    # for idx, vor_poly in enumerate(tqdm.tqdm(voronoi.geoms)):
    #     globaldata[idx + 1]["voronoi"] = vor_poly
    #     globaldata[idx + 1]["voronoi_area"] = vor_poly.area

    with open("partGrid-voronoi", "w+") as the_file:
        the_file.write("{0}\n".format(len(globaldata)))
        for idx in globaldata.keys():
            pt = globaldata[idx]
            nbhs_str = " ".join(map(str, pt["nbhs"]))
            the_file.write("{0} {1} {2} {3} {4} {5} {6} {7} {8} {9}\n".format(pt["x"], pt["y"], pt["left"], pt["right"], pt["type"], pt["geometry"], pt["min_dist"], pt["voronoi_area"], len(pt["nbhs"]), nbhs_str))

    # # Set (current) axis to be equal before showing plot
    # for geom in voronoi.geoms:
    #     plt.plot(*geom.exterior.xy)
        
    # plt.gca().axis("equal")
    # plt.savefig("temp.png")
    # plt.show()
    

    # voronoi = voronoiDiagram4plg(gpd.GeoDataFrame(index=[0], crs=None, geometry=[globalPoints]), multiWalls)
    # voronoi.plot()
    # plt.savefig("temp.png")

def voronoiDiagram4plg(gdf, mask):
	'''Create Voronoi diagram / Thiessen polygons based on polygons.
	
	Parameters:
		gdf: polygons to be used to create Voronoi diagram
			Type: geopandas.GeoDataFrame
		mask: polygon vector used to clip the created Voronoi diagram
			Type: GeoDataFrame, GeoSeries, (Multi)Polygon
	Returns:
		gdf_vd: Thiessen polygons
			Type: geopandas.geodataframe.GeoDataFrame
	'''
	gdf.reset_index(drop=True)
	#convert to shapely.geometry.MultiPolygon
	smp = gdf.unary_union
	#create primary voronoi diagram by invoking shapely.ops.voronoi_diagram (new in Shapely 1.8.dev0)
	smp_vd = svd(smp)
	#convert to GeoSeries and explode to single polygons
	#note that it is NOT supported to GeoDataFrame directly
	gs = gpd.GeoSeries([smp_vd]).explode()
	#convert to GeoDataFrame
	#note that if gdf was shapely.geometry.MultiPolygon, it has no attribute 'crs'
	gdf_vd_primary = gpd.geodataframe.GeoDataFrame(geometry=gs, crs=gdf.crs)
	
	#reset index
	gdf_vd_primary.reset_index(drop=True)	#append(gdf)
	#spatial join by intersecting and dissolve by `index_right`
	gdf_temp = ( gpd.sjoin(gdf_vd_primary, gdf, how='inner', op='intersects')
		.dissolve(by='index_right').reset_index(drop=True) )
	gdf_vd = gpd.clip(gdf_temp, mask)
	gdf_vd = dropHoles(gdf_vd)
	return gdf_vd

def dropHolesBase(plg):
	'''Basic function to remove / drop / fill the holes.
	
	Parameters:
		plg: plg who has holes / empties
			Type: shapely.geometry.MultiPolygon or shapely.geometry.Polygon
	Returns:
		a shapely.geometry.MultiPolygon or shapely.geometry.Polygon object
	'''
	if isinstance(plg, MultiPolygon):
		return MultiPolygon(Polygon(p.exterior) for p in plg)
	elif isinstance(plg, Polygon):
		return Polygon(plg.exterior)

def dropHoles(gdf):
	'''Remove / drop / fill the holes / empties for iterms in GeoDataFrame.
	
	Parameters:
		gdf:
			Type: geopandas.GeoDataFrame
	Returns:
		gdf_nohole: GeoDataFrame without holes
			Type: geopandas.GeoDataFrame
	'''
	gdf_nohole = gpd.GeoDataFrame()
	for g in gdf['geometry']:
		geo = gpd.GeoDataFrame(geometry=gpd.GeoSeries(dropHolesBase(g)))
		gdf_nohole=gdf_nohole.append(geo,ignore_index=True)
	gdf_nohole.rename(columns={gdf_nohole.columns[0]:'geometry'}, inplace=True)
	gdf_nohole.crs = gdf.crs
	return gdf_nohole

def circle(N, R):
    i = np.arange(N)
    theta = i * 2 * np.pi / N
    pts = np.stack([np.cos(theta), np.sin(theta)], axis=1) * R
    seg = np.stack([i, i + 1], axis=1) % N
    return pts, seg

def get_outer_radius(gridData, outer_points):
    outer_radius = 0
    p1x = gridData[outer_points[0]]["x"]
    p1y = gridData[outer_points[0]]["y"]
    for pt in outer_points[1:]:
        p2x = gridData[pt]["x"]
        p2y = gridData[pt]["y"]
        dist = math.sqrt((p1x - p2x)**2 + (p1y - p2y)**2)
        if outer_radius == 0:
            outer_radius = dist
        else:
            outer_radius = max(outer_radius, dist)
    return outer_radius

if __name__ == "__main__":
    main()