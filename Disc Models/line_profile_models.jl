using Gradus
using SpectralFitting
using Plots

# create a new lamppost model called MyLPModel as an example
struct MyLPModel{T} <: AbstractSpectralModel{T,Additive}
    K::T
    a::T
    thickness::T
    height::T
    theta::T
end

function MyLPModel(; K = 1.0, a = 0.998, thickness = 5.0, height = 10.0, theta = 30.0)
    MyLPModel(K, a, thickness, height, theta)
end

function SpectralFitting.invoke!(output, domain, model::MyLPModel)
    m = KerrMetric(one(model.a), model.a)

    x = SVector(0.0, 1000.0, deg2rad(model.theta), 0.0)
    # d = ThinDisc(0.0, Inf)

    r_isco = Gradus.isco(m)
    # the following can easily be made model parameters if desired
    r_tr = 10.0
    r_out = 75.0
    d = ThickDisc() do ρ
        if ρ ≤ r_isco || ρ > r_out
            return -1.0
        else
            if ρ ≤ r_tr
                return model.thickness * (ρ - r_isco) / (r_tr - r_isco)
            else
                return model.thickness * (r_out - ρ) / (r_out - r_tr)
            end
        end
    end

    profile = emissivity_profile(m, d, LampPostModel(;h = model.height))

    _, f = @time lineprofile(m, x, d, profile; bins = domain)
    @. f = model.K * f
    @views output .= f[1:end-1]
end

energy_range = range(0.0, 1.5, length = 500)

# thick disc
example_model = MyLPModel(K = 1.0, a = 0.998, thickness = 10.0, height = 5.0, theta = 30.0)
example_spectrum = invokemodel(energy_range, example_model)
plot(6.4*energy_range[1:end-1], example_spectrum, xrange=[3, 7.5], xlabel = "Energy (keV)", ylabel = "Line flux (arbitrary units)", label = "Thick disc")

# thin disc (not razor thin; can replace with truly thin disc by changing model)
example_model_2 = MyLPModel(K = 1.0, a = 0.998, thickness = 0.5, height = 5.0, theta = 30.0)
example_spectrum_2 = invokemodel(energy_range, example_model_2)
plot!(6.4*energy_range[1:end-1], example_spectrum_2, label = "Thin disc")
