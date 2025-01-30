using Plots, Gradus, SpectralFitting

# universal values (for the entire code simulation)
m = KerrMetric()
r_isco = Gradus.isco(m) # inner radius of rhe disc
r_tr = 20 # radius at which the disc has maximum thickness
r_out = 50.0 # outer radius of the disc

# code for a super-Eddington accretion flow
h_1 = 30.0 # half thickness at that radius

d_1 = ThickDisc() do ρ_1
    if ρ_1 ≤ r_isco || ρ_1 > r_out
        return -1.0
    else
        if ρ_1 ≤ r_tr
            return h_1 * (ρ_1 - r_isco) / (r_tr - r_isco)
        else
            return h_1 * (r_out - ρ_1) / (r_out - r_tr)
        end
    end
end

Gradus.inner_radius(d::typeof(d_1)) = r_isco
model1 = LampPostModel(h = 30.0)
profile1 = emissivity_profile(m, d_1, model1; n_samples = 500)


# code for a less thick accretion disc
h_2 = 20.0

d_2 = ThickDisc() do ρ_2
    if ρ_2 ≤ r_isco || ρ_2 > r_out
        return -1.0
    else
        if ρ_2 ≤ r_tr
            return h_2 * (ρ_2 - r_isco) / (r_tr - r_isco)
        else
            return h_2 * (r_out - ρ_2) / (r_out - r_tr)
        end
    end
end

Gradus.inner_radius(d::typeof(d_2)) = r_isco
model2 = LampPostModel(h = 30.0)
profile2 = emissivity_profile(m, d_2, model2; n_samples = 500)

# code for an lesser accretion disc
h_3 = 10.0
d_3 = ThickDisc() do ρ_3
    if ρ_3 ≤ r_isco || ρ_3 > r_out
        return -1.0
    else
        if ρ_3 ≤ r_tr
            return h_3 * (ρ_3 - r_isco) / (r_tr - r_isco)
        else
            return h_3 * (r_out - ρ_3) / (r_out - r_tr)
        end
    end
end

Gradus.inner_radius(d::typeof(d_3)) = r_isco
model3 = LampPostModel(h = 30.0)
profile3 = emissivity_profile(m, d_3, model2; n_samples = 500)

# code for a thin accretion disc
h_4 = 5.0

d_4 = ThickDisc() do ρ_4
    if ρ_4 ≤ r_isco || ρ_4 > r_out
        return -1.0
    else
        if ρ_4 ≤ r_tr
            return h_4 * (ρ_4 - r_isco) / (r_tr - r_isco)
        else
            return h_4 * (r_out - ρ_4) / (r_out - r_tr)
        end
    end
end

Gradus.inner_radius(d::typeof(d_4)) = r_isco
model4 = LampPostModel(h = 30.0)
profile4 = emissivity_profile(m, d_4, model2; n_samples = 500)

# visualisation of the emissivity profiles for each thickness as a function of disc radius
Plots.plot(profile1, xlabel="Disc Radius (r)", ylabel="Emissivity", title="Emissivity Profile", label = "Height = 30.0m", linecolor = RGB(93/255, 24/255, 33/255))
Plots.plot!(profile2, xlabel="Disc Radius (r)", ylabel="Emissivity", title="Emissivity Profile", label = "Height = 20.0m", linecolor = RGB(72/255, 25/255, 46/255))
Plots.plot!(profile3, xlabel="Disc Radius (r)", ylabel="Emissivity", title="Emissivity Profile", label = "Height = 10.0m", linecolor = RGB(199/255, 78/255, 81/255))
Plots.plot!(profile4, xlabel="Disc Radius (r)", ylabel="Emissivity", title="Emissivity Profile", label = "Height = 5.0m", linecolor = RGB(226/255, 179/255, 194/255))
