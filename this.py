import geopandas as gpd
from shapely.geometry import Point
import json

def calculate_distance(city1, city2, cities_df):
    point1 = cities_df.loc[cities_df['city'] == city1]['geometry'].values[0]
    point2 = cities_df.loc[cities_df['city'] == city2]['geometry'].values[0]
    return point1.distance(point2)

def create_distance_matrix(city_pairs, cities_datafile):
    cities = list(set([city for pair in city_pairs for city in pair]))
    n = len(cities)
    with open(cities_datafile) as f:
        cities_data = json.load(f)
    print(cities_data.columns())
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

