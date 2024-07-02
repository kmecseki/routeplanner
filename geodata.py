#import geopandas as gpd
#from shapely.geometry import Point
#import json
#import pandas as pd
#from pyproj import Geod
#from pyrosm import OSM, get_data
import osmnx as ox
import os.path
#import networkx as nx
#from datetime import timedelta


class GeoData:
    
    def __init__(self, area="Csongrad megye", fname="geodata.graphml"):
        """Downloads geodata or loads it from file if it already exits."""
        
        if os.path.isfile("csm.graphml"):
            self.G = self.load_graph(fname)
        else:
            self.G = self.download_data(area, fname)

    def download_data(self, area, fname):
        """Download geodata around area, example "Csongrad megye"."""

        # Create the graph of the area from OSM data. It will download the
        # data and create a graph
        G = ox.graph_from_place(area, network_type='drive')
        # (For a better accuracy, create a graph with lot more nodes:)
        #G = ox.graph_from_place(graph_area, network_type='drive', simplify=False)

        # OSM data are sometime incomplete so we use the speed module of osmnx
        # to add missing edge speeds and travel times
        G = ox.add_edge_speeds(G)
        G = ox.add_edge_travel_times(G)

        # Save graph to disk if you want to reuse it
        ox.save_graphml(G, fname)
        return G

    def load_graph(self, fname):
        """Loads graph data from disk."""
    
        G = ox.load_graphml(fname)
        return G
        
G = GeoData("Csongrad megye").G
fig, ax = ox.plot_graph(G, figsize=(10, 10), node_size=0, edge_color='y', edge_linewidth=0.2)

#get_route("Szeged", "Szentes")



class geodata:

    def __init__(self, cities_datafile: str) -> None:

        try:
            with open(cities_datafile, encoding = "utf-8") as f:

                #all_geodata_df = pd.read_json(cities_datafile)['features']
                #features_df = pd.json_normalize(all_geodata_df)
                #cities_df = features_df[features_df['geometry.type']=='Point']
                #print(cities_df)
                #self.cities = cities_df

                gdf = gpd.read_file(cities_datafile)
                points_gdf = gdf[gdf.geometry.type == 'Point'] # EPSG:4326
                #print(points_gdf)
                self.cities = points_gdf
            
        except Exception as e:
            print(f"Error opening geospatial datafile: {e}")
        
    


    def calculate_distance(self, city1: str , city2: str) -> float:
        #print(self.cities.columns.tolist())
        point1 = self.cities[self.cities.name == city1]['geometry'].values[0]#).to_crs(epsg=3857) # EPSG:3857 in meters
        #print(point1.x)
        #point10 = gpd.GeoDataFrame(geometry=point1)#.to_crs(epsg=3857) # EPSG:3857 in meters
        #print(point10)
        point2 = self.cities[self.cities.name == city2]['geometry'].values[0]#).to_crs(epsg=3857)
        #point20 = gpd.GeoDataFrame(geometry=point2).to_crs(epsg=3857) # EPSG:3857 in meters
        #print(point1, point10)
        geod = Geod(ellps='WGS84')
        dist_km = geod.inv(point1.y, point1.x, point2.y, point2.x)[2]
        print(dist_km/1000)
        #print(point1.distance(point2))
        
        #return distance

"""
def calculate_distance(city1: str , city2: str) -> float:
    point1 = cities_df.loc[cities_df['city'] == city1]['geometry'].values[0]
    point2 = cities_df.loc[cities_df['city'] == city2]['geometry'].values[0]
    return point1.distance(point2)



def create_distance_matrix(city_pairs, cities_datafile):
    cities = list(set([city for pair in city_pairs for city in pair]))
    n = len(cities)
 

    cities_df = gpd.GeoDataFrame(cities_data, crs='EPSG:4326', geometry=[Point(xy) for xy in cities_data['coordinates']])
    distances = np.zeros((n, n))
    for i in range(n):
        for j in range(i+1, n):
            city1, city2 = cities[i], cities[j]
            if (city1, city2) in city_pairs:
                distance = calculate_distance(city1, city2, cities_df)
                distances[i][j] = distance
                distances[j][i] = distance
    return distances, cities
"""

def clarke_and_wright(num_cars, capacity, distances):

    n = len(distances)
    savings = np.zeros((n, n))
    for i in range(n):
        for j in range(i+1, n):
            savings[i][j] = distances[i][0] + distances[0][j] - distances[i][j]
    
    indices = [(i, j) for i in range(1, n) for j in range(i+1, n)]
    sorted_indices = sorted(indices, key=lambda x: savings[x[0]][x[1]], reverse=True)

    routes = [[0] * (capacity + 1) for _ in range(num_cars)]
    remaining_capacity = [capacity for _ in range(num_cars)]
    penalties = [0 for _ in range(num_cars)]
    distances_traveled = [0 for _ in range(num_cars)]
    empty_distances_traveled = [0 for _ in range(num_cars)]
    
    for (i, j) in sorted_indices:
        for k in range(num_cars):
            if (i in routes[k]) and (j not in routes[k]) and (remaining_capacity[k] >= 1):
                index = routes[k].index(i)
                routes[k].insert(index+1, j)
                distance = distances[i][j]
                distances_traveled[k] += distance
                remaining_capacity[k] -= 1
                if remaining_capacity[k] > 0:
                    penalties[k] = remaining_capacity[k]
                    empty_distances_traveled[k] += distance
                break

            elif (j in routes[k]) and (i not in routes[k]) and (remaining_capacity[k] >= 1):

                index = routes[k].index(j)
                routes[k].insert(index, i)
                distance = distances[i][j]
                distances_traveled[k] += distance
                remaining_capacity[k] -= 1
                if remaining_capacity[k] > 0:
                    penalties[k] = remaining_capacity[k]
                    empty_distances_traveled[k] += distance
                break
    
    for k in range(num_cars):
        for i in range(1, n):
            if i not in routes[k]:
                min_penalty = min(penalties)
                index = penalties.index(min_penalty)
                routes[index].append(i)
                distance = distances[routes[index][-2]][i]
                distances_traveled[index] += distance
                remaining_capacity[index] -= 1
                if remaining_capacity[index] > 0:
                    penalties[index] = remaining_capacity[index]
                    empty_distances_traveled[index] += distance
                else:
                    penalties[index] = 0
    
    return routes, distances_traveled, empty_distances_traveled

