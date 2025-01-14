# Example script showing how to use the package TopologicalAnalysis in 3D.
using TopologicalAnalysis
using Random
using LinearAlgebra


using StatsBase
using CSV
using DataFrames
using Plots

using Distributed # for addprocs

addprocs(12) # add, say, 3 processes
println(nworkers()) # check they have been added

@everywhere using TopologicalAnalysis # make sure all processes are aware of TopologicalAnalysis 

save_directory = "./"


using Random

# Function to read points from CSV and generate random positions
function generate_pts_from_csv(file_path::String)
    # Read data from CSV file
    df = CSV.read(file_path, DataFrame)
    
    pos = [collect(row) for row in eachrow(df)]

    return pos
end

# Function to compute volumes of tetrahedra
function compute_tetrahedron_volumes(simplices, positions)
    num_tetrahedra = size(simplices, 1)
    volumes = zeros(Float64, num_tetrahedra)
    for i in 1:num_tetrahedra
        idx = simplices[i, :]
        v0 = positions[idx[1], :]
        v1 = positions[idx[2], :]
        v2 = positions[idx[3], :]
        v3 = positions[idx[4], :]
        # Compute the volume using the determinant formula
        volumes[i] = abs(det(hcat(v1 - v0, v2 - v0, v3 - v0))) / 6.0
    end
    return volumes
end

# Function to generate random points within the convex hull and return an array of arrays
function generate_points_within_convex_hull(simplices, positions, volumes, N)
    
    random_points = Vector{Vector{Float64}}(undef, N)  # Initialize as a vector of vectors
    num_tetrahedra = size(simplices, 1)
    for i in 1:N
  
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
        random_points[i] = random_point  # Store as a vector
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





function get_random_points_convex(random_csv_path::String)

# Set random seed for reproducibility
#Random.seed!(3828)

# Generate positions from CSV
pos = generate_pts_from_csv(random_csv_path)

# Compute the Delaunay triangulation
delaunay_info = find_delaunay_network(pos, periodic=false, alpha=0)

# Extract simplices and positions
simplices = delaunay_info.simplices  # N x 4 array of indices
positions = p_vals(pos)  # Convert positions to Nx3 array

# Compute volumes of tetrahedra
volumes = compute_tetrahedron_volumes(simplices, positions)

# Number of random points to generate within the convex hull
N = size(positions, 1)

# Generate random points within the convex hull
random_points = generate_points_within_convex_hull(simplices, positions, volumes, N)



delaunay_info = find_delaunay_network(random_points, periodic=false, alpha = 0)
return compute_motifs(delaunay_info)
end







# Path to your CSV file
csv_path = "examples/centers_of_mass.csv"  # Provide the correct path

csv_path_h5 = "examples/centers_of_mass.csv"

#random_points = get_random_points_convex(csv_path)

# First we wish to compute the Delaunay. The function below returns a structure containing
# all information relevant to the delaunay triangulation.
# We tell it that the data is not periodic, and also, by setting alpha=0 we use the default
# value of alpha. We can change this as needed.

#delaunay_info = find_delaunay_network(random_points, periodic=false, alpha = 0)

# delaunay_info is a structure of (custom) type TopologicalNetwork.
# delaunay_info.simplices is a N×4 list of vertices making up the N tetrahedrons of the delauany.
# delaunay_info.not_edge are the vertices identified as interior points (can be useful to plot
#   to make sure the alpha value has been correctly chosen).
# There are other fields, which we can ignore for now.

# Now let us do a more realistic use case.
# We will take 5 realizations of a random distribution of points, and 5 realizations of a perturbed grid

function generate_pts(file_path::String)
    # Read data from CSV file
    df = CSV.read(file_path, DataFrame)

    # Convert DataFrame to array of arrays
    pos = [collect(row) for row in eachrow(df)]
    # Generate motifs for the given data
    delaunay_info = find_delaunay_network(pos, periodic=false, alpha = 0)
    return compute_motifs(delaunay_info)
end



function generate_perturbed_grid()
    # generates a MotifArray for a perturbed grid of points
    pos = [[x,y,z] for x in range(0,stop=1,length=15), y in range(0,stop=1,length=15), z in range(0,stop=1,length=15)][:]
    pos .+= [0.01*randn(3) for i = 1:length(pos)]
    delaunay_inf = find_delaunay_network(pos, periodic=false, alpha = 0)
    return compute_motifs(delaunay_inf)
end


# Function to load networks from .h5 files and convert them to MotifArray
function load_networks_from_h5(h5_files::Vector{String})
    motif_arrays = Vector{MotifArray}()
    for file in h5_files
        topological_network = load(file)
        motif_array = compute_motifs(topological_network)
        push!(motif_arrays, motif_array)
    end
    print("Loaded ", length(motif_arrays), " networks from .h5 files.")
    return motif_arrays
end

# Path to your .h5 files
h5_files = [
 "C:/Users/aur_m/OneDrive - UCB-O365/School/Bee Swarm Research/Code/TopologicalAnalysis2/TopologicalAnalysis_release/data/network/Biofilm/salmonella/Pos11/Pos11_ch2_frame000003_Nz100_data.alpha_120.h5",
 "C:/Users/aur_m/OneDrive - UCB-O365/School/Bee Swarm Research/Code/TopologicalAnalysis2/TopologicalAnalysis_release/data/network/Biofilm/ecoli/Pos6-kde2011/Pos6-kde2011_ch2_frame000015_Nz69_data.alpha_60.h5",
 "C:/Users/aur_m/OneDrive - UCB-O365/School/Bee Swarm Research/Code/TopologicalAnalysis2/TopologicalAnalysis_release/data/network/Yeast/Yeast_03.h5"
 ]  # Add your HDF5 files

# Load the networks from .h5 files
existing_motif_arrays = load_networks_from_h5(h5_files)

# take 5 realizations of the random points and 5 of the perturbed grid in one total array
total_motif_array = vcat([generate_pts(csv_path) for i = 1:5])


# Combine with current motif arrays
total_motif_array = vcat(total_motif_array, existing_motif_arrays...)

# compute the flip graph for these. N.B. writing the ... after total_motif_array is like calling
# compute_flip(total_motif_array[1],total_motif_array[2],...,total_motif_array[10];restrict = 0,thresh=0.5)
# See "splat" operation at https://docs.julialang.org/en/v1/base/base/
flip_graph_combined = compute_flip(total_motif_array...; restrict = 0,thresh=0.5)

# Now compute the combined distance matrix
d = calculate_distance_matrix(flip_graph_combined, total_motif_array, optimal_transport=false)

# To save, create a record of which inputs were random points and which were perturbed grids;
#total_string_record = vcat(["Random_point_" * string(i) for i in 1:5], ["Real Points" * string(i) for i in 1:5], ["Network_" * string(i) for i in 1:length(h5_files)])
total_string_record = vcat(["Bee Swarm" * string(i) for i in 1:5], ["Network_" * string(i) for i in 1:length(h5_files)])



using CSV, DataFrames

# Save the result to a CSV file
CSV.write(save_directory * "result.csv", DataFrame(d, total_string_record))


using Plots, MultivariateStats

# Pass the distance matrix into an MDS algorithm (see MultivariateStats package)
MDS_coords = MultivariateStats.transform(MultivariateStats.fit(MDS, d, maxoutdim=3, distances=true))

# Create both scatter plots with adjusted parameters
p = scatter(MDS_coords[1, 1:5], MDS_coords[2, 1:5],
    label="Bee Swarm",
    xlabel="MDS PC 1", 
    ylabel="MDS PC 2",
    grid=false,
    aspect_ratio=true,
    title="MDS PC1 vs PC2",
    margin=10Plots.mm,  
    fontfamily="Computer Modern",
    tickfontsize=10,
    labelfontsize=12,
    legendfontsize=10,
    titlefontsize=14
)

#scatter!(p, MDS_coords[1, 11:end], MDS_coords[2, 11:end], label="Networks")

scatter!(p, [MDS_coords[1, 6]], [MDS_coords[2, 6]],
    marker=:circle, label="Salmonella"
)
scatter!(p, [MDS_coords[1, 7]], [MDS_coords[2, 7]],
    marker=:circle, label="Ecoli"
)
scatter!(p, [MDS_coords[1, 8]], [MDS_coords[2, 8]],
    marker=:circle, label="Yeast"
)

l = scatter(MDS_coords[2, 1:5], MDS_coords[3, 1:5],
    label="Bee Swarm",
    xlabel="MDS PC 2", 
    ylabel="MDS PC 3",
    grid=false,
    aspect_ratio=true,
    title="MDS PC2 vs PC3",
    margin=10Plots.mm,  
    fontfamily="Computer Modern",
    tickfontsize=10,
    labelfontsize=12,
    legendfontsize=10,
    titlefontsize=14
)

#scatter!(l, MDS_coords[2, 11:end], MDS_coords[3, 11:end], label="Networks")

scatter!(l, [MDS_coords[2, 6]], [MDS_coords[3, 6]],
    marker=:circle, label="Salmonella"
)
scatter!(l, [MDS_coords[2, 7]], [MDS_coords[3, 7]],
    marker=:circle, label="Ecoli"
)
scatter!(l, [MDS_coords[2,8]], [MDS_coords[3, 8]],
    marker=:circle, label="Yeast"
)


combined_plot = plot(p, l, 
    layout=(1, 2), 
    size=(1200, 500),  
    left_margin=20Plots.mm,
    bottom_margin=15Plots.mm,
    plot_title="MDS Embedding Visualization",
    plot_titlefontsize=16
)
display(combined_plot)