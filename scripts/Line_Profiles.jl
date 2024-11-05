using Gradus
using Plots

d = ThinDisc(0.0, 7000.0)
x = SVector(0.0, 100_000.0, deg2rad(12), 0.0)
m = KerrMetric(1.0, 0.998)

# maximal integration radius
minrₑ = 4000.0
maxrₑ = 7000.0

# emissivity function
ε(r) = r^(-3)

# g grid to do flux integration over
gs = range(0.8, 1.2, 500)
_, flux = lineprofile(gs, ε, m, x, d, minrₑ = minrₑ, maxrₑ = maxrₑ, verbose = true)

# plot flux as a function of energy
plot(gs, flux, legend=false, xlabel = "Energy", ylabel = "Flux")
