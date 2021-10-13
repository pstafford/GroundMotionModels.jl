
include("PJSgroundMotions.jl")


function PJSkbc2016(T::U, M::U, Rjb::U, Vs30::U, region::String="Regional") where U<:Real

    if (region != "Regional") && (region != "Others") && (region != "Turkey") && (region != "Italy")
        return PJSgroundMotion(NaN, NaN, NaN, NaN, NaN)
    end

    # Coefficient values
    Ti = [ -1.0, 0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0 ]
    e1 = [ 0.773, 2.982, 3.002, 3.064, 3.128, 3.223, 3.304, 3.757, 3.877, 3.578, 3.482, 3.34, 3.22, 2.998, 2.88, 2.312, 1.684, 1.057, 0.755 ]
    b1 = [ 0.483, -0.363, -0.366, -0.368, -0.378, -0.414, -0.478, -0.666, -0.404, -0.217, 0.107, 0.243, 0.392, 0.667, 0.837, 1.127, 1.079, 1.474, 1.775 ]
    b2 = [ -0.101, -0.195, -0.193, -0.192, -0.183, -0.168, -0.165, -0.232, -0.226, -0.231, -0.226, -0.233, -0.191, -0.169, -0.176, -0.127, -0.159, -0.039, 0.035 ]
    b3 = [ -0.021, -0.406, -0.412, -0.425, -0.440, -0.487, -0.497, -0.341, -0.214, -0.122, -0.042, 0.01, -0.236, -0.178, -0.114, -0.094, -0.222, 0.052, 0.302 ]
    c1 = [ -1.198, -1.231, -1.236, -1.251, -1.267, -1.299, -1.321, -1.342, -1.212, -1.048, -0.966, -0.947, -0.946, -0.972, -0.990, -0.948, -0.911, -0.855, -0.852 ]
    c2 = [ 0.229, 0.272, 0.272, 0.273, 0.278, 0.291, 0.301, 0.295, 0.243, 0.207, 0.159, 0.142, 0.163, 0.144, 0.128, 0.139, 0.162, 0.16, 0.143 ]
    c3 = [ -0.00008, -0.00395, -0.00385, -0.00375, -0.00371, -0.00377, -0.00388, -0.00522, -0.00693, -0.00792, -0.00701, -0.00539, -0.00497, -0.00197, -0.00094, 0.0, 0.0, 0.0, 0.0 ]
    h = [ 5.845, 6.39, 6.425, 6.336, 6.108, 6.096, 6.086, 7.658, 7.468, 6.03, 5.123, 4.75, 4.58, 4.685, 5.392, 4.553, 4.309, 4.365, 4.99 ]
    τi = [ 0.349, 0.35, 0.347, 0.351, 0.348, 0.35, 0.352, 0.375, 0.362, 0.364, 0.357, 0.366, 0.382, 0.382, 0.369, 0.365, 0.36, 0.433, 0.429 ]
    ϕ0i = [ 0.496, 0.451, 0.452, 0.454, 0.461, 0.463, 0.469, 0.459, 0.463, 0.472, 0.503, 0.54, 0.528, 0.541, 0.523, 0.534, 0.553, 0.519, 0.507 ]
    σi = [ 0.683, 0.657, 0.657, 0.661, 0.672, 0.681, 0.704, 0.702, 0.693, 0.7, 0.725, 0.737, 0.767, 0.769, 0.786, 0.806, 0.793, 0.774, 0.683 ]
    Δc3_IT = [ -0.00189, -0.00326, -0.00334, -0.00343, -0.00356, -0.00372, -0.00374, -0.00330, -0.00371, -0.00402, -0.00391, -0.00366, -0.00343, -0.00229, -0.00226, 0.0, 0.0, 0.0, 0.0 ]
    Δc3_Others = [ 0.00142, 0.00326, 0.00341, 0.00349, 0.00364, 0.00371, 0.00378, 0.00347, 0.00338, 0.00348, 0.00308, 0.00296, 0.00234, 0.00175, 0.00186, 0.0, 0.0, 0.0, 0.0 ]
    Δc3_TR = [ 0.0005, 0, -0.00007, -0.00006, -0.00008, 0.00001, -0.00005, -0.00016, 0.00033, 0.00054, 0.00083, 0.0007, 0.00109, 0.00052, 0.00039, 0.0, 0.0, 0.0, 0.0 ]
    # SE_∆c3_IT = [ 0.00077, 0.00079, 0.0008, 0.0008, 0.00081, 0.00082, 0.00083, 0.00084, 0.00084, 0.00084, 0.00085, 0.00089, 0.00089, 0.00088, 0.00088, 0.0, 0.0, 0.0, 0.0 ]
    # SE_∆c3_Others = [ 0.00074, 0.00076, 0.00076, 0.00076, 0.00077, 0.00078, 0.00079, 0.00079, 0.0008, 0.0008, 0.00081, 0.00085, 0.00085, 0.00084, 0.00086, 0.0, 0.0, 0.0, 0.0 ]
    # SE_∆c3_TR = [ 0.00035, 0.00034, 0.00034, 0.00034, 0.00034, 0.00034, 0.00035, 0.00035, 0.00035, 0.00035, 0.00036, 0.00038, 0.00038, 0.00039, 0.00039, 0.0, 0.0, 0.0, 0.0 ]
    g1 = [ 2.188, 1.407, 1.399, 1.382, 1.312, 1.244, 1.163, 0.962, 1.066, 1.207, 1.462, 1.779, 2.236, 2.931, 3.348, 3.395, 3.337, 2.964, 2.707 ]
    g2 = [ -0.364, -0.234, -0.233, -0.230, -0.218, -0.207, -0.194, -0.160, -0.177, -0.200, -0.243, -0.296, -0.373, -0.488, -0.558, -0.566, -0.556, -0.493, -0.451 ]
    ϕS2S = [ 0.314, 0.33, 0.33, 0.332, 0.336, 0.342, 0.35, 0.393, 0.399, 0.359, 0.331, 0.318, 0.342, 0.386, 0.424, 0.446, 0.463, 0.415, 0.397 ]
    Δg1_IT = [ -0.296, -0.360, -0.351, -0.379, -0.376, -0.409, -0.439, -0.344, -0.072, -0.094, -0.109, -0.206, -0.294, -0.329, -0.603, -0.720, -0.952, -0.394, -0.341 ]
    Δg1_Others = [ -1.082, -0.678, -0.663, -0.655, -0.652, -0.606, -0.562, -0.390, -0.341, -0.403, -0.660, -0.777, -1.113, -1.501, -1.581, -1.309, -1.009, -0.831, -0.791 ]
    Δg1_TR = [ 1.378, 1.038, 1.013, 1.034, 1.028, 1.014, 1.001, 0.734, 0.413, 0.497, 0.769, 0.983, 1.406, 1.831, 2.184, 2.028, 1.962, 1.226, 1.021 ]
    Δg2_IT = [ 0.051, 0.063, 0.062, 0.067, 0.066, 0.072, 0.078, 0.061, 0.013, 0.017, 0.019, 0.036, 0.05, 0.057, 0.103, 0.123, 0.162, 0.069, 0.06 ]
    Δg2_Others = [ 0.186, 0.119, 0.116, 0.115, 0.115, 0.107, 0.099, 0.07, 0.062, 0.072, 0.115, 0.135, 0.191, 0.258, 0.271, 0.223, 0.172, 0.145, 0.138 ]
    Δg2_TR = [ -0.236, -0.182, -0.178, -0.182, -0.182, -0.179, -0.177, -0.131, -0.075, -0.089, -0.134, -0.171, -0.242, -0.315, -0.374, -0.346, -0.334, -0.214, -0.178 ]
    # SE_∆g1_IT = [ 0.286, 0.258, 0.256, 0.253, 0.252, 0.255, 0.259, 0.261, 0.212, 0.224, 0.264, 0.261, 0.308, 0.351, 0.392, 0.429, 0.449, 0.354, 0.332 ]
    # SE_∆g1_Others = [ 0.23, 0.212, 0.21, 0.208, 0.207, 0.21, 0.213, 0.22, 0.184, 0.192, 0.217, 0.213, 0.247, 0.281, 0.321, 0.384, 0.419, 0.39, 0.386 ]
    # SE_∆g1_TR = [ 0.349, 0.314, 0.31, 0.307, 0.304, 0.307, 0.311, 0.304, 0.233, 0.256, 0.316, 0.318, 0.375, 0.43, 0.479, 0.505, 0.52, 0.417, 0.387 ]
    # SE_∆g2_IT = [ 0.049, 0.045, 0.045, 0.045, 0.044, 0.045, 0.046, 0.047, 0.039, 0.04, 0.046, 0.045, 0.053, 0.06, 0.067, 0.073, 0.076, 0.062, 0.058 ]
    # SE_∆g2_Others = [ 0.039, 0.037, 0.037, 0.037, 0.037, 0.037, 0.038, 0.039, 0.034, 0.034, 0.038, 0.037, 0.043, 0.048, 0.055, 0.065, 0.071, 0.068, 0.067 ]
    # SE_∆g2_TR = [ 0.06, 0.055, 0.054, 0.054, 0.054, 0.054, 0.055, 0.054, 0.043, 0.046, 0.055, 0.055, 0.065, 0.074, 0.082, 0.086, 0.089, 0.073, 0.068 ]

    # period independent coefficients
    Mh = 6.75
    Mref = 5.5
    Rref = 1.0

    if T < 0.0 # PGV requested
        tidHI = tidLO = 1
    elseif T > 4.0
        return PJSgroundMotion(NaN, NaN, NaN, NaN, NaN)
    else
        tidHI = findfirst(Ti .>= T)
        tidLO = findlast(Ti .<= T)
    end

    if tidLO == tidHI
        # we have a known period
        tid = tidHI

        if region == "Italy"
            Δc3 = Δc3_IT[tid]
            Δg1 = Δg1_IT[tid]
            Δg2 = Δg2_IT[tid]
        elseif region == "Turkey"
            Δc3 = Δc3_TR[tid]
            Δg1 = Δg1_TR[tid]
            Δg2 = Δg2_TR[tid]
        elseif region == "Others"
            Δc3 = Δc3_Others[tid]
            Δg1 = Δg1_Others[tid]
            Δg2 = Δg2_Others[tid]
        else
            Δc3 = 0.0
            Δg1 = 0.0
            Δg2 = 0.0
        end

        # magnitude scaling
        if M < Mh
            FM = b1[tid]*(M - Mh) + b2[tid]*(M - Mh)^2
        else
            FM = b3[tid]*(M - Mh)
        end

        # distance scaling
        FD = (c1[tid] + c2[tid]*(M - Mref)) * log( sqrt(Rjb^2 + h[tid]^2)/Rref ) + (c3[tid]+Δc3)*(sqrt( Rjb^2 + h[tid]^2 ) - Rref)

        # site scaling
        FV = (g1[tid] + Δg1) + (g2[tid] + Δg2)*log(Vs30)

        # median prediction
        lnSa = e1[tid] + FM + FD + FV - log(9.80665)
        Sa = exp(lnSa)

        τ = τi[tid]
        ϕ = sqrt(ϕ0i[tid]^2 + ϕS2S[tid]^2)
        σ = σi[tid]
    else
        # need to interpolate
        @inbounds TiLO = Ti[tidLO]
        @inbounds TiHI = Ti[tidHI]

        if T > 0.0 && T < 0.01
            # assume that PGA coefficients are really at 0.001s
            Tfrac = log(max(T,0.001)/0.001)/log(TiHI/0.001)
        else
            Tfrac = log(T/TiLO)/log(TiHI/TiLO)
        end

        kbcLO = PJSkbc2016(TiLO, M, Rjb, Vs30, region)
        kbcHI = PJSkbc2016(TiHI, M, Rjb, Vs30, region)

        lnSa = kbcLO.lnIM + Tfrac*(kbcHI.lnIM - kbcLO.lnIM)
        Sa = exp(lnSa)

        τ = kbcLO.τ + Tfrac*(kbcHI.τ - kbcLO.τ)
        ϕ = kbcLO.ϕ + Tfrac*(kbcHI.ϕ - kbcLO.ϕ)
        σ = kbcLO.σ + Tfrac*(kbcHI.σ - kbcLO.σ)
    end

    gm = PJSgroundMotion(Sa, lnSa, τ, ϕ, σ)
    return gm
end


function PJSkbc2016(Ti::Vector{U}, M::U, Rjb::U, Vs30::U, region::String="Regional") where U<:Real
    predictions = Array{PJSgroundMotion,1}()
    for T in Ti
        push!(predictions, PJSkbc2016(T, M, Rjb, Vs30, region))
    end
    return predictions
end


function PJSkbc2016(T::U, rup::Rupture{U}, site::Site{U}, region::String="Regional") where U<:Real
    M = rup.magnitude
    Rjb = joyner_boore_distance(rup, site)
    Vs30 = vs30(site)

    return PJSkbc2016(T, M, Rjb, Vs30, region)
end

function PJSkbc2016(Ti::Vector{U}, rup::Rupture{U}, site::Site{U}, region::String="Regional") where U<:Real
    predictions = Array{PJSgroundMotion,1}()
    for T in Ti
        push!(predictions, PJSkbc2016(T, rup, site, region))
    end
    return predictions
end

# PJSkbc2016(0.01, 6.0, 10.0, 450.0, "Regional").IM