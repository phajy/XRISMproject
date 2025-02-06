using SpectralFitting, Gradus, CairoMakie, Turing, Makie, Plots, Colors

##Simulating XRISM data
simulated_energy = range(0.1, stop = 1.6, length = 1000)
simulated_flux = rand(1000) .* exp.(-simulated_energy ./ 2) 

## Shakura Sunyaev Thin Disc Model
shakura_sunyaev_model = (PowerLaw(K = 1e-3, a = 1.7) 
        + BlackBody(K = 1.0, kT = 0.1))
shakura_sunyaev_spectrum = invokemodel(simulated_energy, shakura_sunyaev_model)

## Truncated Model
truncated_model = (PowerLaw(K = 1e-3, a = 1.7) 
        + BlackBody(K = 1.0, kT = 0.02) 
        + PowerLaw(K = 5e-4, a = 2.5))
truncated_spectrum = invokemodel(simulated_energy, truncated_model)

## Super-Eddington Model
super_eddington_model = (PowerLaw(K = 1e-3, a = 2.1) 
        + BlackBody(K = 5.0, kT = 0.15)
        + BlackBody(K = 2.0, kT = 0.1))
super_eddington_spectrum = invokemodel(simulated_energy, super_eddington_model)

## Warped Model
warped_model = (PowerLaw(K = 1e-3, a = 1.7) 
        + BlackBody(K = 1.0, kT = 0.05) 
        + BlackBody(K = 0.8, kT = 0.03)
        + GaussianLine(K = 5e-4, μ = 6.4, σ = 0.2))
warped_spectrum = invokemodel(simulated_energy, warped_model)

energy = range(0.1, stop = 1.6, length = 999)

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

 Makie.save("./different_disk_models.png", fig)
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










## Comparing XSPEC and SF models (WIP)
xs_shakura_sunyaev_model = (XS_PowerLaw(K = 1e-3, a = 1.9) + XS_BlackBody(K = 1.0, kT = 0.05))
xs_shakura_sunyaev_spectrum = invokemodel(simulated_energy, xs_shakura_sunyaev_model)

fig2 = Figure(resolution = (800, 600))
ax2 = Axis(fig2[1, 1], 
        xlabel = "Energy (keV)",
        ylabel = "Flux (counts/m²/s/keV)",
        title = "Comparison of XSPEC and SpectralFitting Models for the Shakura-Sunyaev Model")


line1 = lines!(ax2, energy, shakura_sunyaev_spectrum, color =:deepskyblue)
line2 = lines!(ax2, energy, xs_shakura_sunyaev_spectrum, color =:goldenrod1)

legend = Legend(fig2[1, 2], 
   [line1, line2], 
   ["SpectralFitting", "XSPEC"], position =:rt)
fig2


## Compare Lamp Post Models

function LampPostModel(h)
    return (
        PowerLaw(K = 1e-3, a = 1.8) 
        + BlackBody(K = 1.0, kT = 0.05) 
        + GaussianLine(K = 0.3, μ = 6.4, σ = 0.2))
end

max_lamppost_model = Gradus.LampPostModel(h = 20.0, θ = 50.0, ϕ = 0.0)
min_lamppost_model = Gradus.LampPostModel(h = 3.0, θ = 50.0, ϕ = 0.0)

max_lamppost_spectrum = (simulated_energy, max_lamppost_model)[2]
min_lamppost_spectrum = (simulated_energy, min_lamppost_model)[2]

fig1 = Figure(resolution = (800, 600))
ax1 = Axis(fig1[1, 1], 
        xlabel = "Energy (keV)",
        ylabel = "Flux (counts/m²/s/keV)",
        title = "Comparison of Models for the Simulated XRISM Spectra of the AGN in NGC 4151")

line1 = lines!(ax1, energy, shakura_sunyaev_spectrum, color =:deepskyblue)
line5 = lines!(ax1, energy, max_lamppost_spectrum, color =:dodgerblue4)
line6 = lines!(ax1, energy, min_lamppost_spectrum, color =:palevioletred1)

legend = Legend(fig1[1, 2], 
    [line1, line5, line6], 
    ["Shakura-Sunyaev Model", "Maximum Height Lamp-Post Model", "Minimum Height Lamp-Post Model"])
fig1


## Line Profiles
function calculate_line_profile(m, x, d, bins, plane)
    _, f = lineprofile(
        m, 
        x, 
        d, 
        method = BinningMethod(), 
        # no false images
        callback = domain_upper_hemisphere(),
        verbose = true,
        bins = bins,
        plane = plane,
    )
    return f
end

m = KerrMetric(a = 0.998)
d = 


## Emissivity Profiles (WIP)

d = ShakuraSunyaev(m, )
shakura_sunyaev_profile = emissivity_profile(m, d, shakura_sunyaev_spectrum; n_samples = 500)
