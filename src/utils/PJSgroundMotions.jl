
# __precompile__()
#
# module PJSgroundMotions
#
# using Distributions
# import Distributions.pdf, Distributions.cdf, Distributions.ccdf
#
# export
#     PJSgroundMotion,
#     get_spectrum,
#     pdf,
#     cdf,
#     ccdf


"""
    PJSgroundMotion{T<:Real}

Custom type representing the output of a ground-motion model.
The type contains five fields:
- `IM` is the value of the intensity measure
- `lnIM` is the logarithmic intensity measure value
- `τ` is the between-event standard deviation
- `ϕ` is the within-event standard deviation
- `σ` is the total standard deviation

# Examples
```julia-repl
    gm = PJSgroundMotion(0.1, log(0.1), 0.3, 0.4, 0.5)
```
"""
struct PJSgroundMotion{T<:Real}
    IM::T
    lnIM::T
    τ::T
    ϕ::T
    σ::T
end

PJSgroundMotion() = PJSgroundMotion(NaN, NaN, NaN, NaN, NaN)

"""
    get_spectrum(v::Vector{PJSgroundMotion})

Extract a spectrum from a vector of ground motion instances corresponding to spectral ordinates.
Returns a vector of length equal to the input vector
"""
function get_spectrum(v::Vector{PJSgroundMotion})
    Sai = zeros(length(v),1)
    for i = 1:length(v)
        @inbounds Sai[i] = v[i].IM
    end
    return Sai
end


"""
    distribution(gm::PJSgroundMotion)

Build a normal distribution from the ground-motion instance
"""
function distribution(gm::PJSgroundMotion)
    return Normal(gm.lnIM, gm.σ)
end


"""
    pdf(gm::PJSgroundMotion{T}, lnim::T) where T<:Real

Compute the pdf for the logarithmic IM value
"""
function pdf(gm::PJSgroundMotion{T}, lnim::T) where T<:Real
    d = distribution(gm)
    return pdf(d, lnim)
end


"""
    cdf(gm::PJSgroundMotion{T}, lnim::T) where T<:Real

Compute the cdf for the logarithmic IM value
"""
function cdf(gm::PJSgroundMotion{T}, lnim::T) where T<:Real
    d = distribution(gm)
    return cdf(d, lnim)
end


"""
    ccdf(gm::PJSgroundMotion{T}, lnim::T) where T<:Real

Compute the ccdf for the logarithmic IM value
"""
function ccdf(gm::PJSgroundMotion{T}, lnim::T) where T<:Real
    d = distribution(gm)
    return ccdf(d, lnim)
end

# end
