# using .Speckles
using HDF5
using JLD
using Plots
using ColorSchemes
using Printf
using LaTeXStrings
using StatsPlots
using StatsBase

pyplot()
default(guidefontsize=18, tickfontsize=13, titlefontsize=15, grid=false, colorbar_tickfontsize=13, colorbar_titlefontsize=18, legendfontsize=12)
#%%md Load files

# include("make_matrices.jl")
full_matrix = "/run/media/cocconat/data/granulari/full_matrices"
analysis_path = "/run/media/cocconat/data/granulari/analysis/"
results_path = "final_plots"
!isdir(results_path) && mkdir(results_path)

temporal_file = joinpath(analysis_path,"coarsen10_time.h5")
corr_file = joinpath(analysis_path,"coarsen10_correlation.h5")
norm_file = joinpath(analysis_path,"norm.h5")
times = ["1", "2","5", "10", "20", "40"];

_norm = read(h5open(norm_file,"r")["norm"])
rawdata_file = h5open(joinpath(analysis_path,"coarsen10.h5"),"r")

corr_mats = read(h5open(corr_file))

all_matrix = read(rawdata_file)
speeds = sort(collect(keys(all_matrix)))

fits_file = joinpath(analysis_path,"coarsen10_fits.jld")
times_fits = load(fits_file,"fits")
maxima = load(fits_file,"maxima")
inverse_fits = load(fits_file,"inverse")
shearband = load(fits_file,"shearband")


#=================================================
                    Analysis
=================================================#


## Plot average pixel amplitude



(root,dirs,files) =  first(walkdir(full_matrix))
file = files[2]
_mat = read(h5open(joinpath(root,file),"r"))["matrix"]
heatmap(mean(shift_mat(_mat), dims=3)[:,:,1], c=:haline, xticks=(0:50:250, 0:2:10), yticks=(0:100:400, reverse(0:4:16)), clims=(50,256))
# flat = StatsBase.sum(_mat, dims=3)[:,:,1
length(0:50:250)
##

(root,dirs,files) =  first(walkdir(full_matrix))
plots = []
for (file,speed) in collect(zip(files, speeds))[[1,3,6,9]]
	_mat = read(h5open(joinpath(root,file),"r"))["matrix"]
    flat =StatsBase.mean(shift_mat(_mat), dims=3)[:,:,1]
    push!(plots, heatmap(flat, title="$speed ω (s⁻¹)", c=:haline,titlefontsize=14, xticks=(0:50:250, 0:2:10), yticks=(0:100:400, reverse(0:4:16)), clims=(50,256)))
	if length(plots) <6
		plot!(colorbar=false)
	end
end
plot(plots...)

layout = @layout [
 		grid(2,2) a{0.07w}
 ]

z = ones(256,1)
for x in 1:256
    z[x] = 256-x
end
p = heatmap(reverse(z), ymirror=true, yticks=(0:50:256, 0:50:256), c=:haline, cbar=false, xticks=:none, ylabel="<Pixel intensity>", box=:on,grid=:off)


plot!(plots[3], ylabel="                   z (mm)", xlabel="                            width (mm)")
mean_pixels = plot([plots...;p]..., layout=layout);
plot!(guidefontsize=18, tickfontsize=11)
savefig(mean_pixels,joinpath(results_path,"fig2.pdf"))
plot!()

##
# Verify that the image has an exponential distribution in light-intensity to verify the measure regard actual spekles. : Goodman Statistical Optics
speeds
l = @layout [
    b{0.5w} c{0.5w}
	d{0.5w} e{0.5w}
]
plot(histogram( all_matrix["0.0"][:], norm=true , alpha=1., bins=-0:5:200, title="0.0 ω(s⁻¹)"),
histogram( all_matrix["1.0"][:], norm=true, alpha=1., bins=0:2:200      , title="0.5 ω(s⁻¹)"),
histogram( all_matrix["0.5"][:], norm=true, alpha=1., bins=0:2:200      , title="1.0 ω(s⁻¹)",xlabel="                                Pixel Intensity", ylabel="                    Pixel density", ),
histogram( all_matrix["1.5"][:], norm=true, alpha=1., bins=0:2:200      , title="1.5 ω(s⁻¹)"),
layout=l, lw=3)

intensity = plot!( legend=false)
savefig(intensity, joinpath(results_path,"fig_intensity.pdf"))
# plot!(title="Light ")
plot!(guidefontsize=18, tickfontsize=13)

##

z = ones(1,9)
for x in 1:9
    z[x] = x
end
cs = cgrad(:roma, 1:9)[1/9:1/9:1]

p = plot()
speeds
speeds
for (n,speed) in enumerate(speeds[2:end])
	v = StatsBase.mean(corr_mats["1"][speed][30,:,1:300], dims=1)[:]
	p=plot!(1:300,v, label=false,  xlims=(1,30), c=cs[n], lw=3)
end
p = plot(p, xlabel="Time (s)", ylabel="g(z,t)", lw=3, xticks=(0:7.5:30, 0:0.5:2), guidefontsize=18, tickfontsize=13,titlefontsize=15, title="g(z,t); z = 8mm", grid=false)

heatmap!(z, c=cs[1:end], colorbar=false, title="ω (s⁻¹)", yticks=:none, titlefontsize=13,
    inset_subplots = bbox(0.65, 0.8, 0.3, 0.07, :bottom), subplot=2, axes=false, xticks=(1:2:9, speeds[[2:2:9;9]]), xrotation=-30
)
savefig(p, joinpath(results_path,"fig3.pdf"))
p

##

#
z = ones(1,42)
for x in 1:42
    z[x] = 42-x
end
PALETTE = 42
cs = cgrad(:roma, 1:PALETTE)[collect(1:PALETTE) ./PALETTE]
v = StatsBase.mean(corr_mats["1"]["0.05"][:,:,1:300], dims=2)[:,1,:]

p1 =plot(v', c=cs', xlims=(1,130), legend=false, title="Raw correlations", ylims=(0,1.6), yticks=(0:0.2:1, 0:0.2:1), xticks=(0:15:135, round.(Int,0:15/15:135/15) ) )

heatmap!(z, c=cs, colorbar=false, ylabel="", yticks=:none,
    inset_subplots = bbox(0.65, 0.65, 0.3, 0.1, :bottom), xticks=(1:10:41, 0:4:16), subplot=2, title= "z (mm)", titlefontsize=12)

v = StatsBase.mean(corr_mats["20"]["0.05"][:,:,1:300], dims=2)[:,1,:]
p2 = plot(v', c=cs', xlims=(1,130), legend=false, title="Time averaged", ylims=(0,1.5), yticks=(0:0.5:1, 0:0.5:1), ylabel="     g(z,t)",  xticks=(0:15:135, round.(Int,0:15/15:135/15) ), xlabel="Time (s)" )
decay = plot(p1,p2,layout=(2,1))

plot!(guidefontsize=18, tickfontsize=13, grid=false)
savefig(decay, joinpath(results_path,"fig4.pdf"))
plot!()


## Correlate stretch fit of raw data with gaussian on averaged
plots = []
PALETTE = 9
cs = cgrad(:roma, 1:PALETTE)[collect(1:PALETTE) ./PALETTE]
p1 = plot()

coarse = "1"
for (n,s) in enumerate(speeds[2:end])
    p1 = scatter!(p1,[x.stretched[1:1] for x in times_fits["1"][s]],    [x.gaussian for x in times_fits[coarse][s]], label=false, c=cs[n], title="", msc=cs[n], ms=6)
end
p = plot!(p1, ylabel="Gaussian Exp. - time avg", xlabel="Stretched Exp. - raw data")
#
savefig(p, joinpath(results_path,"figA1.pdf"))
p

##
PALETTE = 9
cs = cgrad(:roma, 1:PALETTE)[collect(1:PALETTE) ./PALETTE]
time_plots = []

z = ones(1,9)
for x in 1:9
    z[x] = x
end
# for t in times
t = "20"
g_p = plot(legend=false, title="Gaussian")
for (s,c) in zip(speeds[2:end],cs)
    data = times_fits[t][s]
    for (n,d) in enumerate(data)
        scatter!(g_p, [42-n],d.gaussian , c=c, ylabel="Gaussian Exp. τ (ms)", xticks=(0:10:40, 0:4:16),  msc=c, ms=6 )
    end
end
#     push!(time_plots,g_p)
# end
# plot(time_plots...)

# fit_plots = plot(time_plots[5], xlabel="z (mm)")
plot!(g_p, xlabel="z (mm)")
heatmap!(z[:,2:end], c=cs[2:end], colorbar=false, title="ω (s⁻¹)", yticks=:none, titlefontsize=13,
    inset_subplots = bbox(0.15, 0.80, 0.3, 0.07, :bottom), subplot=2, axes=false,
	 xticks=(2:2:9, speeds[2:2:9]), xrotation=-30)
savefig(g_p, joinpath(results_path,"figA2.pdf"))
plot!(g_p)
##
ss = parse.(Float64,speeds[2:end])
p = heatmap(speeds[2:end], 1:42, 100*inverse_fits[:,2, :, 5], c=:roma, xlabel="ω (s⁻¹)", ylabel="z (mm)", yticks=(0:10:40, reverse(0:4:16)), colorbar_title="ν (s⁻¹)", cbarfontsize=12, title=" ")
annotate!((8.75, 44, text(" x 10²")))
savefig(p, joinpath(results_path,"fig5.pdf"))
p

##
PALETTE = 9
cs = cgrad(:roma, 1:PALETTE)[collect(1:PALETTE) ./PALETTE]

p = plot()
z = ones(1,9)
for x in 1:9
    z[x] = x
end

coarse = 5
for p in 1:8
    plot!(reverse(1:42), 100*inverse_fits[:,2,p,coarse], c=cs[p], legend=false)
    scatter!(reverse(1:42), 100*inverse_fits[:,2,p,coarse], c=cs[p], legend=false, msc=cs[p], ms=6)
    # plot!(ylims=(0,0.20))
end
fit_all = plot!(ylabel="ν (s⁻¹)", xlabel="z (mm)", title=" ")
annotate!((1.2, 6.6, text(" x 10²")))

plot!(xticks=(0:10:40, 0:4:16))
heatmap!(z, c=cs[1:9], colorbar=false, title="ω (s⁻¹)", yticks=:none, titlefontsize=13,
    inset_subplots = bbox(0.65, 0.8, 0.3, 0.07, :bottom), subplot=2, axes=false, xticks=(1:2:10, speeds[[2:2:9;9]]), xrotation=-30,
)
savefig(p, joinpath(results_path,"fig6.pdf"))
p
##
cs = cgrad(:roma, 1:PALETTE)[collect(1:PALETTE) ./PALETTE]
taus = [inverse_fits[shearband.min[x,1],2,x,5] for x in 1:8]
taus_x = [42 - argmax(1 ./inverse_fits[:,2,x,5]) for x in 1:8]
taus = [maximum(1 ./inverse_fits[:,2,x,5]) for x in 1:8]
p = scatter(taus_x, taus, label="", c=cs, msc=cs, ylabel="maximum τ  (s)", ms=12, xlabel="z (mm)")
plot!(xticks=(0:10:40, 0:4:16), legend=:topleft)
heatmap!(z, c=cs[1:end], colorbar=false, title="ω (s⁻¹)", yticks=:none, titlefontsize=13,
    inset_subplots = bbox(0.65, 0.8, 0.3, 0.07, :bottom), subplot=2, axes=false, xticks=(1:2:10, speeds[[2:2:9;9]]))
savefig(p, joinpath(results_path,"fig7_maxtau.pdf"))

cs = cgrad(:roma, 1:PALETTE)[collect(1:PALETTE) ./PALETTE]
taus_x = [42 - argmax(inverse_fits[:,2,x,5]) for x in 1:8]
taus = 100 .*[maximum(inverse_fits[:,2,x,5]) for x in 1:8]
p = scatter(taus_x, taus, label="", c=cs, msc=cs, ylabel="maximum ν (s⁻¹)", ms=12, xlabel="z (mm)")
plot!(xticks=(0:10:20, 0:4:8), legend=:topleft, xlims=(0:8))
heatmap!(z, c=cs[1:end], colorbar=false, title="ω (s⁻¹)", yticks=:none, titlefontsize=13,
    inset_subplots = bbox(0.65, 0.8, 0.3, 0.07, :bottom), subplot=2, axes=false, xticks=(1:2:10, speeds[[2:2:9;9]]))
savefig(p, joinpath(results_path,"fig7_maxtauinv.pdf"))
plot!(title=" ")
annotate!((0.65, 6.45, text(" x 10²")))

# scatter(1 ./inverse_fits[:,2,:,1])



##
p = groupedbar(xticks=(1:8,string.(speeds)),[-shearband.max -shearband.min], bar_position=:stack,c=[:white cs[6] ], lc=:white, legend=false, ylabel= "z (mm)", xlabel="ω (s⁻¹)", guidefontsize=18, tickfontsize=13, grid=false, ylims=(-42,0),yticks=(reverse(-40:10:0), 0:4:16))
savefig(p, joinpath(results_path,"fig8_bars.pdf"))
##
gr()
layout = @layout [
			a{0.9w} _;_ c{0.0001h}
			]
p = plot(parse.(Float64,speeds[2:end]), shearband.max, seriestypes=[:line, :scatter], c=:red, msc=:red , ms=8, labels=["" "z maximum"],  ylabel="z maximum (mm)", legend=:topleft)
p = plot!(twinx(), parse.(Float64,speeds[2:end]), shearband.min, ylabel="z shear (mm)", label=["" "z shear"], seriestypes=[:line, :scatter], c=:black, msc=:black, ms=8,  legend=:bottomright,xlabel="ω (s⁻¹)")
p =plot(p,plot(frame=:none),layout=layout)
savefig(p, joinpath(results_path,"fig8.pdf"))

##
p = scatter(parse.(Float64,speeds[2:end]), maxima[2,:,5], label="Gaussian Exp. time avg", xlabel="ω (s⁻¹)", ylabel="z (mm)", c=:black, ms=8, shape=:square)
p = scatter!(parse.(Float64,speeds[2:end]), maxima[4,:,1], label="Half width Raw data", xlabel="ω (s⁻¹)", c=:red, msc=:red, ms=8)
plot!(yticks=(30:5:40, reverse(0:2:4)), legendfontsize=12, title=" ", )


savefig(p, joinpath(results_path,"fig9.pdf"))
p
##
p = scatter(maxima[2,:,5], maxima[4,:,1], ylabel="Half width Raw data", xlabel="Gaussian Exp. time avg", c=:black, msc=:black, ms=8, label="")

savefig(p, joinpath(results_path,"fig9_regression.pdf"))


##
plot()

z = ones(1,9)
for x in 1:9
    z[x] = x
end
for (s, c) in zip(1:8,cs)
    plot!(reverse(1:42),inverse_fits[:,func, s, coarse], c=c, lw=2)
   scatter!([banda[s,1]],[inverse_fits[43-banda[s,1], func, s,coarse]],ms=15, shape=:diamond, c=c)
   scatter!([banda[s,2]],[inverse_fits[43-banda[s,2], func, s,coarse]],ms=15, shape=:star5, c=c)
end
plot!(xticks=(0:10:40, 0:4:16), xlabel="z (mm)")

heatmap!(z, c=cs[2:end], colorbar=false, title="ω (s⁻¹)", yticks=:none, titlefontsize=13,
    inset_subplots = bbox(0.65, 0.8, 0.3, 0.07, :bottom), subplot=2, axes=false, xticks=(1:2:10, speeds[[2:2:9;9]]), xrotation=-30
)
plot!(legend=false)
p = plot!(legend=false)
savefig(p, joinpath(results_path,"figA3.pdf"))
p
##
