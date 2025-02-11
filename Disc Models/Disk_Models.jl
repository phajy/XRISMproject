using SpectralFitting, Gradus, CairoMakie, Turing, Makie, Plots, Colors

## Energy Range for Fe Kα Region
energy_range = range(5.5, 7.2, length = 1000)

## Shakura Sunyaev Thin Disc Model
shakura_sunyaev_model = (PowerLaw(K = 8e-4, a = 1.7)  # Comptonization power law
        + BlackBody(K = 1.2, kT = 0.02)  # Multi-temperature blackbody spectrum
        + BlackBody(K = 0.8, kT = 0.01)  # Outer disk component to represent cooler regions
        + GaussianLine(K = 2.5e-4, μ = 6.4, σ = 0.3))  # Broad Fe Kα due to relativistic effects
shakura_sunyaev_spectrum = invokemodel(energy_range, shakura_sunyaev_model)

## ADAF Truncated Model
truncated_model = (PowerLaw(K = 1e-3, a = 1.7)  # Stronger coronal emission
        + BlackBody(K = 0.8, kT = 0.01)  # Cooler disk due to truncation reducing the high-energy part of the blackbody spectrum
        + PowerLaw(K = 5e-4, a = 2.4)  # Additional high-energy component from hot inner region contribition 
        + GaussianLine(K = 1.8e-4, μ = 6.4, σ = 0.1))  # Less prominent or absent broad Fe Kα due to missing emission region closest to the black hole
truncated_spectrum = invokemodel(energy_range, truncated_model)

## Super-Eddington Model
super_eddington_model = (PowerLaw(K = 1.2e-3, a = 2.1)  # Steeper due to high-energy processes
        + BlackBody(K = 6.0, kT = 0.2)  # Hot thick disk component
        + BlackBody(K = 3.0, kT = 0.1)  # Additional thermal emission
        + GaussianLine(K = 3e-4, μ = 6.3, σ = 0.4))  # Extremely broad Fe Kα due to velocity effects
super_eddington_spectrum = invokemodel(energy_range, super_eddington_model)

## Warped Model
warped_model = (PowerLaw(K = 9e-4, a = 1.8)  # Moderate Comptonization
        + BlackBody(K = 1.0, kT = 0.04)  # Varied temperature due to warping
        + BlackBody(K = 0.7, kT = 0.02)  # Secondary temperature variation
        + GaussianLine(K = 2.2e-4, μ = 6.4, σ = 0.2)  # Complex, asymmetric Fe Kα
        + GaussianLine(K = 1.2e-4, μ = 6.5, σ = 0.1))  # Additional component due to warping
warped_spectrum = invokemodel(energy_range, warped_model)

energy = range(5.5, 7.2, length = 999)

## Plot Models
fig = Figure(resolution = (800, 600))
ax = Axis(fig[1, 1], 
        xlabel = "Energy (keV)",
        ylabel = "Flux (counts/m²/s/keV)",
        title = "Comparison of Different Accretion Disk Models for the AGN in NGC 4151")

line1 = lines!(ax, energy, shakura_sunyaev_spectrum, color =:deepskyblue)
line2 = lines!(ax, energy, truncated_spectrum, color =:tomato)
line3 = lines!(ax, energy, super_eddington_spectrum, color =:goldenrod1)
line4 = lines!(ax, energy, warped_spectrum, color =:chartreuse3)

legend = Legend(fig[1, 2], 
   [line1, line2, line3, line4],
   ["Shakura-Sunyaev Model", "Truncated Model", "Super-Eddington Model", 
    "Warped Model"], position =:rt)

 Makie.save("./disk_models.png", fig)
fig




#### EVERYTHING BELOW IS A WIP




## Simulating XRISM data
DATADIR = "data"
spec_path = joinpath(DATADIR, "xa_merged_p0px1000_Hp.pi")
response_path = joinpath(DATADIR, "xa_merged_p0px1000_HpS.rmf")
ancillary_path = joinpath(DATADIR, "rsl_standard_GVclosed.arf")
data = OGIPDataset(spec_path, response=response_path, ancillary=ancillary_path)

normalize!(data)
mask_energies!(data, 6.45, 6.65)
plot(data)

#model.K_2.value = 1.7
sim_data1 = simulate(model, data, exposure_time=1, seed=5)
