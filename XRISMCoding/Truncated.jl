using Gradus, Plots, Colors, SpectralFitting, CairoMakie, Makie

m = KerrMetric(1.0, 0.998)
x = SVector(0.0, 1000.0, deg2rad(50), 0.0)

# get the innermost stable circular orbit for this spacetime
r_isco = Gradus.isco(m)
# radius at which the disc has maximum thickness
r_tr = 15.0
# half thickness at that radius
h = 10.0
# outer radius of the disc
r_out = 50.0

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
    αlims = (-55, 55), 
    βlims = (-40, 40),
    verbose = true,
    pf = PointFunction((m, gp, t) -> gp.x[2] * cos(gp.x[3])) ∘ ConstPointFunctions.filter_intersected()
)

fig = Figure(resolution = (800, 600), backgroundcolor =:transparent)
ax = Axis(fig[1, 1], 
        title = "Hot Inner Region Truncated Accretion Disc",
        backgroundcolor =:transparent)

CairoMakie.heatmap!(α, β, img', colormap=:linear_kry_5_95_c72_n256)
CairoMakie.contour!(α, β, img', color=:black, levels = 10)
#Makie.save("./stationarytruncated80.png", fig)
fig

