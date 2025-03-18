using Gradus
using SpectralFitting
using XSPECModels
using Plots
using Printf
using Statistics

DATADIR = "data"
spec_path = joinpath(DATADIR, "xa_merged_p0px1000_Hp.pi")
response_path = joinpath(DATADIR, "xa_merged_p0px1000_HpS.rmf")
ancillary_path = joinpath(DATADIR, "rsl_standard_GVclosed.arf")
data = OGIPDataset(spec_path, response = response_path, ancillary = ancillary_path)

normalize!(data)
mask_energies!(data, 1.0, 8.0)

# create a new lamppost model called MyLPModel as an example
struct MyLPModel{T} <: AbstractSpectralModel{T,Additive}
    K::T
    a::T
    height::T
    theta::T
    E_line::T
end

function MyLPModel(; K = 1.0, a = 0.998, height = 10.0, theta = 30.0, E_line = 6.4)
    MyLPModel(K, a, height, theta, E_line)
end

function SpectralFitting.invoke!(output, domain, model::MyLPModel)
    m = KerrMetric(one(model.a), model.a)

    x = SVector(0.0, 1000.0, deg2rad(model.theta), 0.0)
    d = ThinDisc(0.0, Inf)

    profile = emissivity_profile(m, d, LampPostModel(;h = model.height))

    _, f = @time lineprofile(m, x, d, profile; bins = domain ./ model.E_line)
    @. f = model.K * f
    f[1] = 0.0
    @views output .= f[1:end-1]
end

model =
    # XS_Laor(K = FitParam(1.0e-2, frozen = true), lineE = FitParam(6.4, lower_limit = 6.0, upper_limit = 7.0), a = FitParam(0.998, frozen = true), θ = FitParam(30.0)) +
    MyLPModel(K = FitParam(8.0e-3, frozen = true), a = FitParam(0.998, frozen = true), height = FitParam(10.0, frozen = true), theta = FitParam(30.0, frozen = true), E_line = FitParam(6.4, frozen = true)) +
    PowerLaw(a = FitParam(2.29), K = FitParam(1.0e-3))

domain = collect(range(1.0, 8.0, 300))
flux = invokemodel(domain, model)
plot(domain[1:end-1], flux)

# note that fitting does not work because there is an issue with automatic differentiation
# prob = FittingProblem(model => data)
# result = SpectralFitting.fit!(prob, LevenbergMarquadt())
# model

# plot(data)
# plot!(result)

sims = @time simulate(model, data.data.response, data.data.ancillary; exposure_time = 1e6, seed = 42)
plot(sims, xlims = (2.0, 8.0))

# this can also be rebinned by hand if desired (this will be automated in future)
x = sims.output_domain[1:end-1]
y = sims.data
y_err = sims.variance

function rebin(x, y, y_err, n)
    len = length(x)
    new_len = div(len, n)
    new_x = [mean(x[i*n-n+1:i*n]) for i in 1:new_len]
    new_y = [mean(y[i*n-n+1:i*n]) for i in 1:new_len]
    new_y_err = [sqrt(sum(y_err[i*n-n+1:i*n].^2)) for i in 1:new_len]
    return new_x, new_y, new_y_err
end

n = 8   # example rebinning factor
x_rebinned, y_rebinned, y_err_rebinned = rebin(x, y, y_err, n)
plot(x_rebinned, y_rebinned, xlims = (3.0, 7.0), legend = false, xlabel="Energy (keV)", ylabel="Photons/cm^2/s/keV")
