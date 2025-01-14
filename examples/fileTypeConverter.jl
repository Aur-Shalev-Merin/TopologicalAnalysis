using HDF5
using CSV
using DataFrames

function h5_to_csv(h5_filename::String, output_dir::String=".")
    mkpath(output_dir)
    h5_file = h5open(h5_filename, "r")
    base_name = splitext(basename(h5_filename))[1]
    
    for dataset_path in ["Type", "dim", "idx", "r", "regions", "tvec"]
        try
            data = read(h5_file[dataset_path])
            
            # Convert to DataFrame based on data type
            df = if isa(data, String)
                DataFrame(value = [data])
            elseif isa(data, Number)
                DataFrame(value = [data])
            elseif isa(data, Vector)
                DataFrame(value = data)
            elseif isa(data, Matrix)
                DataFrame(data, :auto)
            else
                throw(ArgumentError("Unsupported data type: $(typeof(data))"))
            end
            
            csv_filename = joinpath(output_dir, "$(base_name)_$(dataset_path).csv")
            CSV.write(csv_filename, df)
            println("Converted dataset '$(dataset_path)' to $(csv_filename)")
        catch e
            println("Warning: Could not convert dataset '$(dataset_path)': $(e)")
        end
    end
    
    close(h5_file)
    println("Conversion complete!")
end

# Example usage
input_file = "C:/Users/aur_m/OneDrive - UCB-O365/School/Bee Swarm Research/Code/TopologicalAnalysis2/TopologicalAnalysis_release/data/motif/Cheng2019/fish_nuclei_pos_2_region.h5"
h5_to_csv(input_file, "output")