

function PJSidriss2014( T::U, M::U, Rrup::U, Vs30::U, Frv::Int64 ) where U<:Real
    # coefficients vary above/below critical magnitude Mc
    Mc = 6.75

    # model limits the scaling of Vs30, and is only applicable for Vs30>450 m/s
    if Vs30 < 450.0
        print("PJS Error: Idriss (2014) model is only applicable for Vs30 >= 450 m/s")
        return PJSgroundMotion()
    end

    Vs30 = ( Vs30 > 1200.0 ) ? 1200.0 : Vs30

    # coefficients
    Ti = [ 0.0, 0.01, 0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 7.5, 10.0 ]
    α1_lo = [ 7.0887, 7.0887, 7.1157, 7.2087, 6.2638, 5.9051, 7.5791, 8.019, 9.2812, 9.5804, 9.8912, 9.5342, 9.2142, 8.3517, 7.0453, 5.1307, 3.361, 0.1784, -2.4301, -4.357, -7.8275, -9.2857 ]
    α2_lo = [ 0.2058, 0.2058, 0.2058, 0.2058, 0.0625, 0.1128, 0.0848, 0.1713, 0.1041, 0.0875, 0.0003, 0.0027, 0.0399, 0.0689, 0.16, 0.2429, 0.3966, 0.756, 0.9283, 1.1209, 1.4016, 1.5574 ]
    β1_lo = [ 2.9935, 2.9935, 2.9935, 2.9935, 2.8664, 2.9406, 3.019, 2.7871, 2.8611, 2.8289, 2.8423, 2.83, 2.856, 2.7544, 2.7339, 2.68, 2.6837, 2.6907, 2.5782, 2.5468, 2.4478, 2.3922 ]
    β2_lo = [ -0.2287, -0.2287, -0.2287, -0.2287, -0.2418, -0.2513, -0.2516, -0.2236, -0.2229, -0.22, -0.2284, -0.2318, -0.2337, -0.2392, -0.2398, -0.2417, -0.245, -0.2389, -0.2514, -0.2541, -0.2593, -0.2586 ]
    α1_hi = [ 9.0138, 9.0138, 9.0408, 9.1338, 7.9837, 7.756, 9.4252, 9.6242, 11.13, 11.3629, 11.7818, 11.6097, 11.4484, 10.9065, 9.8565, 8.3363, 6.8656, 4.1178, 1.8102, 0.0977, -3.0563, -4.4387 ]
    α2_hi = [ -0.0794, -0.0794, -0.0794, -0.0794, -0.1923, -0.1614, -0.1887, -0.0665, -0.1698, -0.1766, -0.2798, -0.3048, -0.2911, -0.3097, -0.2565, -0.232, -0.1226, 0.1724, 0.3001, 0.4609, 0.6948, 0.8393 ]
    β1_hi = [ 2.9935, 2.9935, 2.9935, 2.9935, 2.7995, 2.8143, 2.8131, 2.4091, 2.4938, 2.3773, 2.3772, 2.3413, 2.3477, 2.2042, 2.1493, 2.0408, 2.0013, 1.9408, 1.7763, 1.703, 1.5212, 1.4195 ]
    β2_hi = [ -0.2287, -0.2287, -0.2287, -0.2287, -0.2319, -0.2326, -0.2211, -0.1676, -0.1685, -0.1531, -0.1595, -0.1594, -0.1584, -0.1577, -0.1532, -0.147, -0.1439, -0.1278, -0.1326, -0.1291, -0.122, -0.1145 ]
    α3 = [ 0.0589, 0.0589, 0.0589, 0.0589, 0.0417, 0.0527, 0.0442, 0.0329, 0.0188, 0.0095, -0.0039, -0.0133, -0.0224, -0.0267, -0.0198, -0.0367, -0.0291, -0.0214, -0.024, -0.0202, -0.0219, -0.0035 ]
    ξ = [ -0.854, -0.854, -0.854, -0.854, -0.631, -0.591, -0.757, -0.911, -0.998, -1.042, -1.03, -1.019, -1.023, -1.056, -1.009, -0.898, -0.851, -0.761, -0.675, -0.629, -0.531, -0.586 ]
    γ = [ -0.0027, -0.0027, -0.0027, -0.0027, -0.0061, -0.0056, -0.0042, -0.0046, -0.003, -0.0028, -0.0029, -0.0028, -0.0021, -0.0029, -0.0032, -0.0033, -0.0032, -0.0031, -0.0051, -0.0059, -0.0057, -0.0061 ]
    φ = [ 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.06, 0.04, 0.02, 0.02, 0.0, 0.0, 0.0, 0.0 ]

    # check bounds on the allowable value of T (note that coefficients for PGA=0.01s)
    T = ( T < 0.01 ) ? 0.01 : ( T > 10.0 ) ? 10.0 : T
    tidHI = findfirst(Ti .>= T)
    tidLO = findlast(Ti .<= T)

    if tidLO == tidHI
        # we have a known period
        tid = tidHI

        if M <= Mc
            lnSa = α1_lo[tid] + α2_lo[tid]*M + α3[tid]*(8.5 - M)^2 - (β1_lo[tid] + β2_lo[tid]*M) * log( Rrup + 10.0 ) + ξ[tid]*log(Vs30) + γ[tid]*Rrup + φ[tid]*Frv
        else
            lnSa = α1_hi[tid] + α2_hi[tid]*M + α3[tid]*(8.5 - M)^2 - (β1_hi[tid] + β2_hi[tid]*M) * log( Rrup + 10.0 ) + ξ[tid]*log(Vs30) + γ[tid]*Rrup + φ[tid]*Frv
        end

        Sa = exp( lnSa )

        τ = NaN
        ϕ = NaN
        σ = 1.18 + 0.035*log(max(min(T, 3.0), 0.05)) - 0.06*max(min(M, 7.5), 5.0)

    else
        # need to interpolate
        @inbounds TiLO = Ti[tidLO]
        @inbounds TiHI = Ti[tidHI]
        Tfrac = log(T/TiLO)/log(TiHI/TiLO)

        iLO = PJSidriss2014(TiLO, M, Rrup, Vs30, Frv)
        iHI = PJSidriss2014(TiHI, M, Rrup, Vs30, Frv)

        lnSa = iLO.lnIM + Tfrac*(iHI.lnIM-iLO.lnIM)
        Sa = exp( lnSa )

        τ = NaN
        ϕ = NaN
        σ = 1.18 + 0.035*log(max(min(T, 3.0), 0.05)) - 0.06*max(min(M, 7.5), 5.0)
    end

    gm = PJSgroundMotion(Sa, lnSa, τ, ϕ, σ)
    return gm
end

# @time PJSidriss2014(0.01, 6.0, 10.0, 760.0, 0)

function PJSidriss2014( Ti::Vector{U}, M::U, Rrup::U, Vs30::U, Frv::Int64) where U<:Real
    predictions = Array{PJSgroundMotion,1}()
    for T in Ti
        push!(predictions, PJSidriss2014(T, M, Rrup, Vs30, Frv))
    end
    return predictions
end


function PJSidriss2014(T::U, rup::Rupture{U}, site::Site{U}) where U<:Real
	Rrup = rupture_distance(rup, site)
	Vs30 = vs30(site)
	Fss, Fnm, Frv, Fuk = mechanism(rup)

	gm = PJSidriss2014(T, rup.magnitude, Rrup, Vs30, Frv )
	return gm
end


function PJSidriss2014(Ti::Vector{U}, rup::Rupture{U}, site::Site{U}) where U<:Real
	predictions = Array{PJSgroundMotion,1}()
	for T in Ti
		push!(predictions, PJSidriss2014(T, rup, site))
	end
	return predictions
end
