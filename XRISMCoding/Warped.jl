using Gradus, Plots, Colors, CairoMakie, Makie

m = KerrMetric(1.0, 0.998)
x = SVector(0.0, 1000.0, deg2rad(60), 0.0)


# get the innermost stable circular orbit for this spacetime
r_isco = Gradus.isco(m)
r_tr = 30.0   # Radius at which the disc has maximum thickness
h = 10.0      # Base half-thickness at r_tr
r_out = 50.0  # Outer radius of the disc
theta_0 = π / 3  # Maximum warp angle
k = π / 6       # Warp frequency

# define the disc shape -- return a negative number 
# where the disc should not be intersected, else the cross 
# sectional height


d = ThickDisc() do ρ
    if ρ ≤ r_isco || ρ > r_out
        return -1.0  # No disc beyond valid range
    else
        # Warp angle (oscillates with radius)
        theta = theta_0 * sin(k * (ρ - r_isco))

        # Define thickness function
        if ρ ≤ r_tr
            f_r = (ρ - r_isco) / (r_tr - r_isco)
        else
            f_r = (r_out - ρ) / (r_out - r_tr)
        end
        
        # Apply warping via cosine adjustment
        return h * max(0, f_r) * cos(theta)
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
    βlims = (-35, 35),
    verbose = true,
    pf = PointFunction((m, gp, t) -> gp.x[2] * cos(gp.x[3])) ∘ ConstPointFunctions.filter_intersected()
)

fig = Figure(resolution = (800, 600), backgroundcolor =:transparent)
ax = Axis(fig[1, 1], 
        title = "Warped Accretion Disc",
        backgroundcolor =:transparent)

hm = CairoMakie.heatmap!(α, β, img', colormap=:linear_kry_5_95_c72_n256)
CairoMakie.contour!(α, β, img', color=:black, levels = 10)
Colorbar(fig[:, end+1], hm)
#Makie.save("./warped801.png", fig)
fig
