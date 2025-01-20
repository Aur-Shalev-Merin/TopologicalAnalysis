# Import necessary packages
using TopologicalAnalysis
using Random
using LinearAlgebra
using Plots
using CSV
using DataFrames

using Distributed  # For parallel computing

addprocs(12)  # Add 12 worker processes
println(nworkers())  # Confirm the number of workers
@everywhere using TopologicalAnalysis  # Load package on all workers

save_directory = "./"

# Function to read points from CSV and generate random positions
function generate_pts_from_csv(file_path::String)
    # Read data from CSV file
    df = CSV.read(file_path, DataFrame)

    pos = [collect(row) for row in eachrow(df)]

    return pos
end

# Function to generate random points within the convex hull without uniform distribution
function generate_points_within_convex_hull_random(simplices, positions, N)
    random_points = Vector{Vector{Float64}}(undef, N)
    num_tetrahedra = size(simplices, 1)
    for i in 1:N
        # Randomly select a tetrahedron uniformly (without volume weighting)
        tet_index = rand(1:num_tetrahedra)
        idx = simplices[tet_index, :]

        # Get the vertices of the selected tetrahedron
        v0 = positions[idx[1], :]
        v1 = positions[idx[2], :]
        v2 = positions[idx[3], :]
        v3 = positions[idx[4], :]

        # Generate a random point within the tetrahedron using barycentric coordinates
        s, t, u = rand(3)
        if s + t + u > 1
            s = 1 - s
            t = 1 - t
            u = 1 - u
        end
        w = 1 - s - t - u

        random_point = s * v0 + t * v1 + u * v2 + w * v3
        random_points[i] = random_point
    end
    return random_points
end

# Function to convert vector of vectors to matrix
function p_vals(P::Vector)
    return permutedims(hcat(P...))
end

function p_vals(P::Matrix)
    return P
end

function get_random_points_convex_random(random_csv_path::String)
    # Set random seed for reproducibility
    Random.seed!(1234)

    # Generate positions from CSV
    pos = generate_pts_from_csv(random_csv_path)

    # Compute the Delaunay triangulation
    delaunay_info = find_delaunay_network(pos, periodic=false, alpha=0)

    # Extract simplices and positions
    simplices = delaunay_info.simplices  # N x 4 array of indices
    positions = p_vals(pos)  # Convert positions to Nx3 array

    # Number of random points to generate within the convex hull
    N = size(positions, 1)

    # Generate random points within the convex hull without uniform distribution
    random_points = generate_points_within_convex_hull_random(simplices, positions, N)

    return random_points  # This is now an array of arrays
end

# Main script
random_csv_path = "examples/centers_of_mass.csv" 
random_points = get_random_points_convex_random(random_csv_path)

# Convert positions to arrays for plotting (if needed)
positions = generate_pts_from_csv(random_csv_path)

# Prepare data for visualization
positions_matrix = p_vals(positions)
random_points_matrix = p_vals(random_points)

# Visualize the original points and the random points
scatter3d(positions_matrix[:, 1], positions_matrix[:, 2], positions_matrix[:, 3], label="Original Points", markersize=2)
scatter3d!(random_points_matrix[:, 1], random_points_matrix[:, 2], random_points_matrix[:, 3], label="Random Points", markersize=1, alpha=0.5)
xlabel!("X")
ylabel!("Y")
zlabel!("Z")
title!("Random Points Within Convex Hull (Non-Uniform Distribution)")

display(current())

# Check the type
println(typeof(random_points))  # Output: Vector{Vector{Float64}}

# Inspect the first few points
for i in 1:5
    println("Point $(i): ", random_points[i])
end
