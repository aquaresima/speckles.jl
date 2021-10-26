include("Speckles.jl")
using .Speckles
using HDF5
using JLD
using LinearAlgebra
using StatsBase
using RollingFunctions
using ProgressBars
using FFTW

#%%md
# Import matrices from npy format.
# Set the location of the folder and import all files.
##
full_matrices = "/run/media/cocconat/data/granulari/full_matrices"
analysis_path = "/run/media/cocconat/data/granulari/analysis/"

(root,dirs,files) =  first(walkdir(full_matrices))
files

pyplot()
default(guidefontsize=18, tickfontsize=13, titlefontsize=15)

#%%md
# Analysis of the temporal variance in order to select the best coarse-grain values. (slow ∼ 5 mins)

file = joinpath(full_matrix,"V_0.00.h5")
vars = measure_temporal_variance(file)
h5open(joinpath(analysis_path,"variance.h5"),"w") do fid
    fid["variance"] = vars
end
plot(var, ylabel="Avg luminiosity StdDev ", xlabel="Coarse grain size (px)")

#
# %% md Spatial coarse grain of each matrix. Metapixels have size 10x10 pixels, from variance analysis (slow ∼ 5 mins)

function shift_mat(mat)
	a = Array{Int16}(mat)
	a[a.<0] .= abs.(a[a.<0]) .+127
	return a
end

coarse10_file = joinpath(analysis_path,"coarsen10.h5")

@printf "Coarse grain and compute time-averaged heatmaps"

isfile(coarse10_file) && rm(coarse10_file)
fid = h5open(coarse10_file,"w")
for _file in files
    !endswith(_file,"h5") && continue
    file = joinpath(full_matrix,_file)
    data = read(h5open(joinpath(file),"r"))
    @show _file, data["speed"]
    matrix = shift_mat(data["matrix"])
    m = Speckles.coarsen(matrix, x=10,y=10)
    fid[string(data["speed"])] = m
end

## Compute normalization of the beam

@printf "Compute the normalization factor"
(root,dirs,files) =  first(walkdir(full_matrix))
norm_matrix = zeros(420,250)
for file in files
    full = mean(read(h5open(joinpath(full_matrix,file),"r"))["matrix"], dims=3)[:,:,1]
    average = copy(full)
    target_mean = mean(full)
    average[full .< (target_mean -50)]
    height = mean(average, dims=2)[:,1]
    width = mean(average, dims=1)[1,:]
    height_norm = height/ maximum(height)|> x -> 1 ./x
    width_norm =  width/maximum(width) |> x-> 1 ./x
    global norm_matrix+= width_norm' .* ones(size(full)) .*height_norm
end
norm_matrix ./= length(files)
_norm = Speckles.coarsen(norm_matrix, x=10,y=10)

norm_file = joinpath(analysis_path,"norm.h5")
h5open(norm_file,"w") do fid
    fid["norm"] = _norm
end




## Temporal coarsegrain

fid = h5open(coarse10_file,"r")
all_matrix = read(fid)
speeds = collect(keys(all_matrix))
temporal_file = joinpath(analysis_path,"coarsen10_time1.h5")
times = ["1", "2","5", "10", "20", "40"];

total = []
for t in times
    temporal_file = joinpath(analysis_path,"coarsen10_time$t.h5")
    isfile(temporal_file) && rm(temporal_file)
    fid = h5open(temporal_file,"w")
    dict = []
    for k in speeds
        @show k,t
        matrix = round.(Int16,all_matrix[k])
        m = Speckles.smoothen_t(matrix, t=parse(Int, t))
        # m = Speckles.coarsen(matrix, t=parse(Int, t))
        fid[k] = m
    end
    close(fid)
end

##
#====================================
        Correlate signal
====================================#

corr_mats = Dict(t=>[] for t in times)
corr_file = joinpath(analysis_path,"coarsen10_correlation.h5")
isfile(corr_file) && rm(corr_file)
fid =h5open(corr_file, "w")

for t in times
    temporal_file = joinpath(analysis_path,"coarsen10_time$t.h5")
    all_matrix =read(h5open(temporal_file, "r"))
    fid_t = create_group(fid, t)
    for k in speeds
        @show  k, t
        m = Speckles.compute_correlation(all_matrix[k])
        push!(corr_mats[t], m)
        fid_t[k] = m
    end
end
close(fid)

## Fit exponentials
using LsqFit
using HypothesisTests
@. gauss_exp(x,p) =exp(-(x/p[1])^2)
@. stretch_exp(x,p) = exp(-(x/abs(p[1]))^abs(p[2]))
@. fit_cosh(x,p) = 2-cosh(x/p[1])
function halfheight(data)
        x = findfirst(x->x<data[1]/2., data)
        return  isnothing(x) ? [length(data)] : [x]
end

"""
Fit decays correlation with three curves:
	- Gaussian_Exponential
	- Hyperbolic Cosine
	- Stretched exponential
and half-width time of the function.

Parameters:
----------
	matrix: correlation matrix
	y: pixel depth
	px: number of points per fit

Return
------

"""
function fit_correlations(corr_mat, no_points=10)
    tt = no_points
    speed_fits = Dict()
    corr_depths = []
    for s in speeds
        speed_fits[s] = []
        mysample = mean(corr_mat[s], dims=2)[:,1,:]
        for y in 1:size(mysample)[1] # number of rows
            ## run all fits
            fits = (
            stretched = curve_fit(stretch_exp, 3tt:4tt, mysample[y,tt:2tt], [tt,1.]).param,
            gaussian  = curve_fit(gauss_exp, 1:tt, mysample[y,1:tt], [1.]).param,
            cosh = curve_fit(fit_cosh, 1:tt, mysample[y,1:tt], [1.]).param,
            half      = halfheight(mysample[y,:])
            )
            push!(speed_fits[s],fits)
        end
    end
    return speed_fits
end

corr_mats["1"]

times_fits = Dict()
for t in times
	@show t
    speed_fits = fit_correlations(corr_mats[t])
    push!(times_fits, (t=>speed_fits))
end

times_fits

## Find maxima of inverse tau

maxima = Dict()
times
inverse_fits = zeros(42,4,8,6) # rows, function, speed, coarse
maxima = zeros(4,8,6)
for (k,t) in enumerate(times)
    for (j,s) in enumerate(speeds[2:end])
        data = times_fits[t][s]
        for (m,f) in enumerate([:stretched, :gaussian, :cosh, :half] )
            values = Vector{Real}()
            for (n,d) in enumerate(data)
                push!(values,1. /getfield(d,f)[1])
                inverse_fits[n,m,j,k] = 1. / getfield(d,f)[1]
            end
            # values = runmean(values,10)
            maxima[m,j,k] = argmax(values)
        end
    end
end



func = 2
mobilmean = 15
coarse = 5

banda = zeros(Int,8,3)
for s in 1:8
    new_val = inverse_fits[:,func,s,coarse]
    deriv = diff(reverse(inverse_fits[:,func,s,coarse])*10)
    mymin = argmin(runmean(deriv,mobilmean))
    mymax = argmax(reverse(inverse_fits[:,func,s,coarse]))
    banda[s,1] = mymin
    banda[s,2] = mymax
end

shearband = (min=banda[:,1], max=banda[:,2])
fits_file = joinpath(analysis_path,"coarsen10_fits.jld")
isfile(fits_file) && rm(fits_file)
save(fits_file,"fits", times_fits, "maxima", maxima,        "inverse", inverse_fits, "shearband", shearband)
