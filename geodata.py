#import geopandas as gpd
#from shapely.geometry import Point
#import json
#import pandas as pd
#from pyproj import Geod
#from pyrosm import OSM, get_data
import osmnx as ox
import os.path

class GeoData:
    
    def __init__(self, area="Csongrad megye", fname="geodata.graphml"):
        """Downloads geodata or loads it from file if it already exits."""
        
        if os.path.isfile(fname):
            self.G = self.load_graph(fname)
            print("Re-using data.")
        else:
            self.G = self.download_data(area, fname)
            print("Downloading data.")

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