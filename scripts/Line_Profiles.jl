using Gradus
using Plots

d = ThinDisc(0.0, Inf)
x_broad = SVector(0.0, 100_000.0, deg2rad(19), 0.0)
x_int = SVector(0.0, 100_000.0, deg2rad(34), 0.0)
x_narr = SVector(0.0, 100_000.0, deg2rad(12), 0.0)
m = KerrMetric(1.0, 0.998)

# maximal integration radius for broad line
minrₑ_broad = 10
maxrₑ_broad = 150

# maximal integration radius for intermediate line
minrₑ_int = 2000
maxrₑ_int = 4000

# maximal integration radius for narrow line
minrₑ_narr = 4000
maxrₑ_narr = 7000

# emissivity function
ε(r) = r^(-3)

# g grid to do flux integration for broad line
gs_broad = range(0.8, 3.4, 500)
_, flux_broad = lineprofile(gs_broad, ε, m, x_broad, d, minrₑ = minrₑ_broad, maxrₑ = maxrₑ_broad, verbose = true)

# g grid to do flux integration for intermediate line
gs_int = range(0.8, 1.5, 500)
_, flux_int = lineprofile(gs_int, ε, m, x_int, d, minrₑ = minrₑ_int, maxrₑ = maxrₑ_int, verbose = true)

# g grid to do flux integration for narrow line
gs_narr = range(0.8, 2.1, 500)
_, flux_narr = lineprofile(gs_narr, ε, m, x_narr, d, minrₑ = minrₑ_narr, maxrₑ = maxrₑ_narr, verbose = true)

# plot flux as a function of energy for broad line
plot(gs_broad, flux_broad, legend=false, xlabel = "Energy", ylabel = "Flux")

# plot flux as a function of energy for intermediate line
plot!(gs_int, flux_int, legend=false, xlabel = "Energy", ylabel = "Flux")

# plot flux as a function of energy for narrow line
plot!(gs_narr, flux_narr, legend=false, xlabel = "Energy", ylabel = "Flux")
