import geodata
import route
import networkx as nx
import osmnx as ox
from datetime import timedelta
import itertools

city_pairs = [['Derekegyhaz', 'Szentes'], ['Szentes', 'Csongrad'], ['Szentes', 'Derekegyhaz'], ['Szentes', 'Szeged'], ['Szentes', 'Szeged'], ['Szeged', 'Szentes']]
numer_of_drivers = 2
drivers = range(numer_of_drivers)
capacity = 2

def assign_drivers(city_pairs, drivers):
    assignment = list(itertools.product(drivers, repeat=len(city_pairs)))
    assignments = []
    for assig in assignment:
        assignments.append(assig)
    return assignments

assign_drivers(city_pairs, drivers)

assignment = (1, 0, 1, 0, 1, 0)

routes = route.Route(assignment)



driver = 1
points = [x for i, x in enumerate(city_pairs) if assignment[i] == driver]

print(points)

def assign_numbers(n):
    result = []
    for i in range(n):
        start = 2 * i
        end = start + 1
        result.append([start,end])
        print(result)
    return result

points2 = assign_numbers(len(points))
flat_list = [item for sublist in points2 for item in sublist]

# Function to check if a permutation maintains the order within sublists
def is_valid_permutation(permutation, sublists):
    indices = {item: i for i, item in enumerate(permutation)}
    for sublist in sublists:
        for i in range(len(sublist) - 1):
            if indices[sublist[i]] > indices[sublist[i + 1]]:
                return False
    return True

# Generate all permutations of the flat list
all_permutations = itertools.permutations(flat_list)

# Filter permutations to keep only valid ones
valid_permutations = [perm for perm in all_permutations if is_valid_permutation(perm, points2)]

routes = []
flat_cities = [x for xs in points for x in xs]

# Print all valid permutations
for perm in valid_permutations:
    print(perm)
    single_route = []
    for k in perm:
        single_route.append(flat_cities[k])
    routes.append(single_route)
print(routes)
exit(0)

# Print the total number of valid permutations
print(f"Total number of valid permutations: {len(valid_permutations)}")


exit(0)

G = geodata.GeoData("Csongrad megye").G
#fig, ax = ox.plot_graph(G, figsize=(10, 10), node_size=0, edge_color='y', edge_linewidth=0.2)

origin_coordinates = ox.geocode("Derekegyhaz")
destination_coordinates = ox.geocode("Szeged")

# In the graph, get the nodes closest to the points
origin_node = ox.nearest_nodes(G, Y=origin_coordinates[0], X=origin_coordinates[1])
destination_node = ox.nearest_nodes(G, Y=destination_coordinates[0], X=destination_coordinates[1])

# Get the shortest route by distance
shortest_route_by_distance = ox.shortest_path(G, origin_node, destination_node, weight='length')

# Plot the shortest route by distance
fig, ax = ox.plot_graph_route(G, shortest_route_by_distance, route_color='y', route_linewidth=6, node_size=0)

# Get the shortest route by travel time
shortest_route_by_travel_time = ox.shortest_path(G, origin_node, destination_node, weight='travel_time')

# Plot the shortest route by travel time
fig, ax = ox.plot_graph_route(G, shortest_route_by_travel_time, route_color='y', route_linewidth=6, node_size=0)

# Plot the 2 routes
fig, ax = ox.plot_graph_routes(G, routes=[shortest_route_by_distance, shortest_route_by_travel_time], route_colors=['r', 'y'], route_linewidth=6, node_size=0)
# Get the travel time, in seconds
# Note here that we use "nx" (networkx), not "ox" (osmnx)
travel_time_in_seconds = nx.shortest_path_length(G, origin_node, destination_node, weight='travel_time')
print("travel time in seconds", travel_time_in_seconds)

#The travel time in "HOURS:MINUTES:SECONDS" format
travel_time_in_hours_minutes_seconds = str(timedelta(seconds=travel_time_in_seconds))
print("travel time in hours minutes seconds", travel_time_in_hours_minutes_seconds)

# Get the distance in meters
distance_in_meters = nx.shortest_path_length(G, origin_node, destination_node, weight='length')
print("distance in meters", distance_in_meters)
# Distance in kilometers
distance_in_kilometers = distance_in_meters / 1000
print("distance in kilometers", distance_in_kilometers)
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

