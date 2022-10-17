using DrWatson
@quickactivate "granular_media"

using Revise
using speckles
import YAML
using JLD2
using HDF5
using ThreadSafeDicts 
using ProgressBars
import speckles: shift_mat, coarsen

# Set the location of the folder and import all files.
paths = YAML.load_file(projectdir("paths.yml"))
full_matrices  = paths["full_matrices"]
analysis_path  = paths["analysis"]

# Set coarse dimensions
conf = YAML.load_file(projectdir("conf.yml"))
coarse = conf["coarse"]
x = coarse["x"]
y = coarse["y"]
times = coarse["times"]

# Set filenames
coarse_file = joinpath(analysis_path,"coarsen$(x)x$(y).jld2") 
filename = split(coarse_file,".")[1]
temporal_file(root,t) = joinpath(analysis_path,"$(root)_$t.jld2")
temporal_file(filename,1)
##

## %% md Spatial coarse grain of each matrix. Metapixels have size 10x10 pixels, from variance analysis (slow ∼ 5 mins)
@info "Coarse grain and compute time-averaged heatmaps"
data, _ = produce_or_load(full_matrices, filename=coarse_file) do full_matrices
    files = [joinpath(full_matrices,f)  for f in filter(endswith(".h5"),readdir(full_matrices))]
    lk = Threads.ReentrantLock()
    matrices = Vector{Array{Float64,3}}(undef,length(files))
    speeds = Vector{String}(undef,length(files))
    # Threads.@threads  
    for n in eachindex(files)
        file = files[n]
        @info "$file to $( Threads.threadid())"
        _fid = lock(lk) do
            return  read(h5open(joinpath(file),"r"))
        end
        @info "$file loaded"
        matrices[n] = coarsen(shift_mat(_fid["matrix"]), x=x, y=y)
        speeds[n] = string(_fid["speed"])
        @info "$file coarsened"
    end
    return @strdict matrices speeds x y
end


## Compute normalization of the beam

## Temporal coarsegrain
data = tosymboldict(load(coarse_file)) |> dict2ntuple;
@unpack x,y, speeds, matrices= data
for t in ProgressBar(times)
    file = temporal_file(filename,t)
    temp_coarse =Vector{Array{Int16,3}}(undef,length(speeds))
    data, _ = produce_or_load(matrices, filename=file) do matrices 
        for k in eachindex(speeds)
            @debug "time: $t, speed: $(speeds[k])"
            m = speckles.smoothen_t(round.(Int16,matrices[k]), t= t)
            temp_coarse[k] = m
        end
        return @strdict matrices=temp_coarse speeds=speeds times=t
    end
end;


##
#====================================
        Correlate signal
====================================#
corr_file = joinpath(analysis_path,"$(filename)_correlations.jld2")
data, _ = produce_or_load(times,filename=corr_file) do times
    correlations = Matrix{Array{Float64,3}}(undef,length(times), length(speeds))
    iter = ProgressBar(eachindex(times))
    for t in iter
        temporal_file = joinpath(analysis_path,"$(filename)_$(times[t]).jld2")
        set_postfix(iter,Time="$(times[t])")
        temp_matrix =load(temporal_file)["matrices"]
        temp_coarse =Vector{Array{Int16,3}}(undef,length(speeds))
        for k in eachindex(speeds)
            @debug  speeds[k], times[t]
            m = speckles.compute_correlation(temp_matrix[k])
            correlations[t,k] =  m
        end
    end
    return @strdict correlations speeds=speeds times=times
end;
#   #

data["correlations"][1,1]