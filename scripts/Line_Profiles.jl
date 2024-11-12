using Gradus
using Plots

d = ThinDisc(0.0, Inf)
x_broad = SVector(0.0, 100_000.0, deg2rad(19), 0.0)
x_int = SVector(0.0, 100_000.0, deg2rad(34), 0.0)
x_narr = SVector(0.0, 100_000.0, deg2rad(12), 0.0)
m = KerrMetric(1.0, 0.998)

# maximal integration radius for broad line and its normalisation
minrₑ_broad = 10
maxrₑ_broad = 150
K_broad = 3.2

# maximal integration radius for intermediate line and its normalisation
minrₑ_int = 2000
maxrₑ_int = 4000
K_int = 2.1

# maximal integration radius for narrow line and its normalisation
minrₑ_narr = 4000
maxrₑ_narr = 7000
K_narr = 1.5

# emissivity function
ε(r) = r^(-3)

# g grid to do flux integration for broad line
gs_broad = range(0.7, 1.1, 500)
_, flux_broad = lineprofile(gs_broad, ε, m, x_broad, d, minrₑ = minrₑ_broad, maxrₑ = maxrₑ_broad, verbose = true)
flux_broad = K_broad * flux_broad

# g grid to do flux integration for intermediate line
gs_int = range(0.7, 1.1, 500)
_, flux_int = lineprofile(gs_int, ε, m, x_int, d, minrₑ = minrₑ_int, maxrₑ = maxrₑ_int, verbose = true)
flux_int = K_int * flux_int

# g grid to do flux integration for narrow line
gs_narr = range(0.7, 1.1, 500)
_, flux_narr = lineprofile(gs_narr, ε, m, x_narr, d, minrₑ = minrₑ_narr, maxrₑ = maxrₑ_narr, verbose = true)
flux_narr = K_narr * flux_narr

# plot flux as a function of energy for broad line
plot(gs_broad, flux_broad, xlabel = "Energy", ylabel = "Flux", xrange = (0.75, 1.05), label = "Broad line")

# plot flux as a function of energy for intermediate line
plot!(gs_int, flux_int, xlabel = "Energy", ylabel = "Flux", label = "Intermediate line")

# plot flux as a function of energy for narrow line
plot!(gs_narr, flux_narr, xlabel = "Energy", ylabel = "Flux", label = "Narrow line")
