import this

city_pairs = [['Derekegyhaz', 'Szentes'], ['Szentes', 'Csongrad'], ['Szentes', 'Derekegyhaz'], ['Szentes', 'Szeged'], ['Szentes', 'Szeged']]

distances, cities = this.create_distance_matrix(city_pairs, 'export.geojson')

num_cars = 2
capacity = 2

routes, distances_traveled, empty_distances_traveled = this.clarke_and_wright

