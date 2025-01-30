using SpectralFitting
using Printf
using Plots
using Turing
using CairoMakie
using PairPlots
using PythonCall

data_path = "data"

# Function to prepare spectral data
function prepare!(ds, min, max)
    regroup!(ds)
    drop_bad_channels!(ds)
    mask_energies!(ds, min, max)
    subtract_background!(ds)
    normalize!(ds)
    ds
end
function prepare_xmm(spec, rmf, arf)
    ds = XmmData(XmmEPIC(), spec, response = rmf, ancillary = arf)
    prepare!(ds, 3.0, 10.0)
end
function prepare_nustar(data_path, obsid, fpm)
    ds = NuStarData(joinpath(data_path, "nu$(obsid)$(fpm)01_sr_min20.pha"))
    prepare!(ds, 3.0, 50.0)
end

# Process XMM data
xmm_data = []
for i = 2:4
    spec = joinpath(data_path, "xa_merged_p0px1000_Hp.fits")
    rmf = joinpath(data_path, "xa_merged_p0px1000_HpS.rmf")
    arf = joinpath(data_path, "rsl_standard_GVclosed.arf")
    push!(xmm_data, prepare_xmm(spec, rmf, arf))
end

# Model and Simulation
begin
    model =
        PowerLaw(a = FitParam(0.3)) +
        GaussianLine(K = FitParam(0.02), μ = FitParam(6.4), σ = FitParam(0.05))
    model.K_2.value = 1.7

    sim_data1 = simulate(model, xmm_data[1], exposure_time = 1, seed = 5)
    model.K_1.value = 0.03
    model.K_2.value = 3.0
    sim_data2 = simulate(model, xmm_data[2], exposure_time = 2, seed = 7)

    model.K_1.value = 0.05
    model.K_2.value = 2.3
    sim_data3 = simulate(model, xmm_data[3], exposure_time = 1, seed = 9)
end
