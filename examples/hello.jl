# Load necessary packages and modules
using CSV
using DataFrames

# Include the topological analysis tools
using TopologicalAnalysis

# Load the CSV data
data = CSV.read("C:/Users/aur_m/Downloads/centers_of_mass (2).csv", DataFrame)

# Extract the points (x, y, z coordinates)
points = [data.x data.y data.z]

# Check the size of points to ensure consistency
println("Points array size: ", size(points))

# Adjust the 'r' value (you may need to experiment with this value)
r = 5  # or any value based on your understanding of the dataset

# Compute the flip graph
flip_graph = compute_flip_graph(points; r=r)

# Print flip graph details
println("Flip Graph:")
println(flip_graph)
# Now, calculate the distance matrix for the FlipGraph
# Check if optimal_transport or JS options are needed depending on your use case
distance_matrix = calculate_distance_matrix(flip_graph, points)

# Print the distance matrix
println("Distance Matrix:")
println(distance_matrix)

# Optionally, perform Delaunay triangulation as well
delaunay_network = find_delaunay_network(points)

# Print Delaunay network details
println("Delaunay Triangulation Network:")
println(delaunay_network)
