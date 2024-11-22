using Gradus, Plots

m = KerrMetric(1.0, 0.998)
x = SVector(0.0, 1000.0, deg2rad(65), 0.0)

# get the innermost stable circular orbit for this spacetime
r_isco = Gradus.isco(m)
# radius at which the disc has maximum thickness
r_tr = 5.0
# half thickness at that radius
h = 4.0
# outer radius of the disc
r_out = 15.0

# define the disc shape -- return a negative number 
# where the disc should not be intersected, else the cross 
# sectional height

d = ThickDisc() do ρ
    if ρ ≤ r_isco || ρ > r_out
        return -1.0
    else
        if ρ ≤ r_tr
            return h * (ρ - r_isco) / (r_tr - r_isco)
        else
            return h * (r_out - ρ) / (r_out - r_tr)
        end
    end
end

Gradus.inner_radius(d::typeof(d)) = r_isco

pl_int = interpolate_plunging_velocities(m)
redshift = interpolate_redshift(pl_int, x)
pf = redshift ∘ ConstPointFunctions.filter_intersected()

# and then render as usual
α, β, img = rendergeodesics(
    m,
    x,
    d,
    2000.0,
    αlims = (-25, 25), 
    βlims = (-15, 18),
    verbose = true,
    pf = pf
)

heatmap(α, β, img, aspect_ratio=1)

# specify an emissivity profile
ε(r) = r^(-3)

# g grid
gs = range(0.1, 1.5, 150)

# do flux integration for broad line
 _, flux_broad = lineprofile(gs, ε, m, x, d, verbose = true, method = BinningMethod())

 plot(gs, flux_broad)

 # it would be interesting to compare with the equivalent thin disc
 _, flux_thin_disc = lineprofile(gs, ε, m, x, ThinDisc(r_isco, r_out), verbose = true, method = BinningMethod())
 plot!(gs, flux_thin_disc)
