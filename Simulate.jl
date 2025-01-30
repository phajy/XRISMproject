using SpectralFitting
using Printf
using Plots
using Turing
using CairoMakie

# Process XRISM data
DATADIR = "data"
spec_path = joinpath(DATADIR, "xa_merged_p0px1000_Hp.pi")
response_path = joinpath(DATADIR, "xa_merged_p0px1000_HpS.rmf")
ancillary_path = joinpath(DATADIR, "rsl_standard_GVclosed.arf")
data = OGIPDataset(spec_path, response=response_path, ancillary=ancillary_path)

normalize!(data)
mask_energies!(data, 6.45, 6.65)
plot(data)

# Model and Simulation
model =
    PowerLaw(a=FitParam(0.3)) +
    GaussianLine(K=FitParam(0.02), μ=FitParam(6.4), σ=FitParam(0.05))
model.K_2.value = 1.7

sim_data1 = simulate(model, data, exposure_time=1, seed=5)
