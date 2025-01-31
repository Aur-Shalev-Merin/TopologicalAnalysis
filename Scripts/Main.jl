using TopologicalAnalysis
using Random
using LinearAlgebra
using StatsBase
using CSV
using DataFrames
using Plots
using Distributed
using Glob
using FilePathsBase: basename, dirname
using MultivariateStats  # For MDS
using Dates  # For timing

println("Starting Main.jl script...")

# Stop all current workers (if they exist)
if !isempty(workers())
    rmprocs(workers())
end

addprocs(12)
println("Number of workers: ", nworkers())

@everywhere using TopologicalAnalysis

save_directory = "./"

println("Defining functions...")

function generate_pts(file_path::String)
    println("Generating points from CSV at: ", file_path)
    df = CSV.read(file_path, DataFrame)
    pos = [collect(row) for row in eachrow(df)]
    delaunay_info = find_delaunay_network(pos, periodic=false, alpha=0)
    return compute_motifs(delaunay_info)
end

@everywhere function load_network(file::String)
    println("Loading file on worker ", myid(), ": ", file)
    net = load(file)
    motifs = compute_motifs(net)
    println("Finished loading and computing motifs for file on worker ", myid(), ": ", file)
    return motifs
end

function load_networks_from_h5(h5_files::Vector{String})
    println("Starting parallel loading of .h5 files on ", nworkers(), " workers...")
    start_time = now()
    motif_arrays = pmap(load_network, h5_files)
    println("Finished parallel loading of .h5 files in ", now() - start_time)
    return motif_arrays
end

println("Main code starts...")

# 1) Path to CSV file (for the single Bee Points dataset)
csv_path = joinpath(@__DIR__, "../BeeData/centers_of_mass.csv")
println("Real points CSV path: ", csv_path)

# 2) Recursively collect all .h5 files
network_dir = joinpath(@__DIR__, "../network")
println("Collecting .h5 files from subdirectories, but only the first directory per each of the first 5 directories...")

all_h5_files = collect(glob("**/*.h5", network_dir))

# Group files by their immediate parent directory
files_by_dir = Dict{String, Vector{String}}()
for f in all_h5_files
    parent_dir_name = basename(dirname(f))
    push!(get!(files_by_dir, parent_dir_name, String[]), f)
end

# Sort directory names, then keep only first 5
sorted_dirs = sort(collect(keys(files_by_dir)))
first_5_dirs = sorted_dirs[1:min(2, length(sorted_dirs))]

# Collect only the first file from each of these 5 directories
selected_h5_files = String[]
for d in first_5_dirs
    file_list = files_by_dir[d]
    sort!(file_list)  # sort by filename if desired
    push!(selected_h5_files, first(file_list))
end

h5_files = selected_h5_files
println("Directories found: ", length(files_by_dir))
println("Selected first 5 directories: ", first_5_dirs)
println("Total selected files: ", length(h5_files))

# 3) Extract directory-based labels for each file
println("Extracting directory-based labels from selected files...")
function dir_label(file_path::String)
    return basename(dirname(file_path))  # immediate directory name
end
network_labels = [dir_label(path) for path in h5_files]
println("Selected labels: ", network_labels)

# 4) Load MotifArrays
start_time = now()
existing_motif_arrays = load_networks_from_h5(h5_files)
println("Loaded MotifArrays in ", now() - start_time)

# 5) Create a single Real Points motif
println("Generating Bee Points motif...")
start_time = now()
real_pts = generate_pts(csv_path)
println("Generated Bee Points motif in ", now() - start_time)

# 6) Combine real points + network arrays
println("Combining Bee swarm points with existing motif arrays...")
start_time = now()
total_motif_array = vcat(real_pts, existing_motif_arrays...)
println("Combined motifs in ", now() - start_time)

# 7) Compute the flip graph and distance matrix
println("Computing flip graph...")
start_time = now()
@time flip_graph_combined = compute_flip(total_motif_array...; restrict=0, thresh=0.5)
println("Flip graph computation complete in ", now() - start_time)

println("Computing distance matrix...")
start_time = now()
d = calculate_distance_matrix(flip_graph_combined, total_motif_array, optimal_transport=false)
println("Computed distance matrix in ", now() - start_time)

# 8) Build labels: first "Real Points," then each directory label
all_labels = vcat(["Bee Swarm"], network_labels)
println("Total labels: ", length(all_labels))

# 9) Save results to CSV
println("Saving results to CSV: result.csv")
start_time = now()
df = DataFrame(d, :auto)
rename!(df, Symbol.(all_labels); makeunique=true)
CSV.write(save_directory * "result.csv", df)
println("Saved results to CSV in ", now() - start_time)

println("Moving on to MDS visualization...")

using MultivariateStats

# 10) MDS for visualization
start_time = now()
MDS_coords = MultivariateStats.transform(MultivariateStats.fit(MDS, d, maxoutdim=3, distances=true))
println("MDS computation complete in ", now() - start_time)

# 11) Plot in 2D using MDS PC2 vs PC3
p = scatter(
    [MDS_coords[1,1]], [MDS_coords[2,1]],
    label=all_labels[1],
    xlabel="MDS PC 1",
    ylabel="MDS PC 2",
    grid=false,
    aspect_ratio=:equal,
    title="MDS PC1 vs PC2"
)
for i in 2:length(all_labels)
    scatter!(
        p,
        [MDS_coords[1, i]],
        [MDS_coords[2, i]],
        marker=:circle,
        label=all_labels[i]
    )
end
p2 = scatter(
    [MDS_coords[2,1]], [MDS_coords[3,1]],
    label=all_labels[1],
    xlabel="MDS PC 2",
    ylabel="MDS PC 3",
    grid=false,
    aspect_ratio=:equal,
    title="MDS PC2 vs PC3"
)
for i in 2:length(all_labels)
    scatter!(
        p2,
        [MDS_coords[2, i]],
        [MDS_coords[3, i]],
        marker=:circle,
        label=all_labels[i]
    )
end

# 12) Combine both plots side by side
combined_plot = plot(
    p, p2,
    layout = (1, 2),
    size = (1200, 500),
    left_margin = 10Plots.mm,
    bottom_margin = 10Plots.mm,
    plot_title = "Comparison of MDS Projections",
    plot_titlefontsize = 14
)
display(combined_plot)

println("Main.jl script complete.")