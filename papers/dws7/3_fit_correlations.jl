using DrWatson
@quickactivate "granular_media"

using Revise
using speckles
import YAML
using JLD2
using UnPack

using LsqFit
using HypothesisTests
using Statistics
using RollingFunctions

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

# Set filename]
coarse_file = joinpath(analysis_path,"coarsen$(x)x$(y).jld2") 
filename = String(split(coarse_file,".")[1])
temporal_file(root,t) = joinpath(analysis_path,"$(root)_$t.jld2")
temporal_file(filename,1)

## Fit exponentials
functions = [:stretched, :gaussian, :cosh, :half] 
N_points = 10

@. gauss_exp(x,p) =exp(-(x/p[1])^2)
@. stretch_exp(x,p) = exp(-(x/abs(p[1]))^abs(p[2]))
@. fit_cosh(x,p) = 2-cosh(x/p[1])
function halfheight(data)
        x = findfirst(x->x<data[1]/2., data)
        return  isnothing(x) ? [length(data)] : [x]
end
function get_fits(data)
    global N_points
    tt = N_points
    return (
    stretched = curve_fit(stretch_exp, 3tt:4tt, data[tt:2tt], [tt,1.]).param,
    gaussian  = curve_fit(gauss_exp,   1:tt, data[1:tt], [1.]).param,
    cosh      = curve_fit(fit_cosh,    1:tt, data[1:tt], [1.]).param,
    half      = halfheight(data)
    )
end

corr_file = joinpath(analysis_path,"$(filename)_correlations.jld2")
rm(joinpath(analysis_path,"$(filename)_fits.jld2"))
fits,_ = produce_or_load(corr_file, filename=joinpath(analysis_path,"$(filename)_fits.jld2")) do corr_file
    data = load(corr_file) |> tosymboldict |> dict2ntuple;
    @unpack correlations, speeds, times = data 
    ys = size(correlations[1,1],1)
    global functions
    fits = Array{Any,3}(undef, size(times,1), size(speeds,1),ys)
    inverse_fits = Array{Float64,4}(undef, size(functions,1), ys,size(times,1), size(speeds,1))
    for s in eachindex(speeds)
        for t in eachindex(times)
            mysample = mean(correlations[t,s], dims=2)[:,1,:]
            @info "Run fit $(speeds[s]) $(times[t]), matrix: $(size(mysample))"
            for y in 1:size(mysample)[1] # number of rows
                _fits = @views get_fits(mysample[y,:])
                ## run all fits
                fits[t,s,y] = _fits
                inverse_fits[:,y, t,s] = [ 1. / getfield(_fits,f)[1] for f in functions]
            end
        end
    end
    return @strdict fits inverse_fits times speeds
end

data = load(joinpath(analysis_path,"$(filename)_fits.jld2")) |> tosymboldict |> dict2ntuple;
@unpack fits, inverse_fits, times, speeds = data

## Get bandwidth for specific values
fit_params = conf["fit"]
func = fit_params["func"]
mobilmean = fit_params["mobilmean"]
coarse_time = fit_params["coarse_time"]
config = (@strdict func coarse_time mobilmean) |> dict2ntuple

data,_ = produce_or_load(config,prefix=filename) do config
    @unpack func, coarse_time, mobilmean = config
    banda = zeros(Int,8,3)
    for s in 1:8
        t = findfirst(x->x==coarse_time, times)
        f = findfirst(x->x==func, functions)
        local data = inverse_fits[f,:,t,s]
        deriv = diff(reverse(data))
        mymin = argmin(runmean(deriv,mobilmean))
        mymax = argmax(reverse(data))
        banda[s,1] = mymin
        banda[s,2] = mymax
    end
    return @strdict banda func=func coarse_time=coarse_time mobilmean=mobilmean
end

