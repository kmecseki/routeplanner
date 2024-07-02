import geodata

#city_pairs = [['Derekegyhaz', 'Szentes'], ['Szentes', 'Csongrad'], ['Szentes', 'Derekegyhaz'], ['Szentes', 'Szeged'], ['Szentes', 'Szeged']]


city1 = "Csongrád"
city2 = "Szeged"

cities_data = geodata.geodata("export.geojson")

print(cities_data.calculate_distance(city1, city2))



#distances, cities = this.create_distance_matrix(city_pairs, 'export.geojson')

#num_cars = 2
#capacity = 2

#routes, distances_traveled, empty_distances_traveled = this.clarke_and_wright

