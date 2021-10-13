
include("PJSgroundMotions.jl")


function PJSsadigh1997(T::U, M::U, R::U, FaultType::Int64, SiteType::Int64) where U<:Real

    if T > 4.0
        return PJSgroundMotion(NaN, NaN, NaN, NaN, NaN)
    end

    # define the common period vector
    Ti = [ 0.0, 0.07, 0.10, 0.20, 0.30, 0.40, 0.50, 0.75, 1.00, 1.50, 2.00, 3.00, 4.00 ]

	mechFactor = 1.0
	if SiteType == 0
		# rock site
		if M <= 6.5
			c1 = [ -0.624, 0.110, 0.275, 0.153, -0.057, -0.298, -0.588, -1.208, -1.705, -2.407, -2.945, -3.700, -4.230 ]
	        c2 = ones(13)
	        c3 = [ 0.000, 0.006, 0.006, -0.004, -0.017, -0.028, -0.040, -0.050, -0.055, -0.065, -0.070, -0.080, -0.100 ]
	        c4 = [ -2.100, -2.128, -2.148, -2.080, -2.028, -1.990, -1.945, -1.865, -1.800, -1.725, -1.670, -1.610, -1.570 ]
	        c5 = 1.29649 * ones(13)
	        c6 = 0.250 * ones(13)
	        c7 = [ 0.0, -0.082, -0.041, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ]
	        ca = [ 1.39, 1.40, 1.41, 1.43, 1.45, 1.48, 1.50, 1.52, 1.53, 1.53, 1.53, 1.53, 1.53 ]
	        cb = -0.14 * ones(13)
	        cc = [ 0.38, 0.39, 0.40, 0.42, 0.44, 0.47, 0.49, 0.51, 0.52, 0.52, 0.52, 0.52, 0.52 ]
		else
			c1 = [ -1.274, -0.540, -0.375, -0.497, -0.707, -0.948, -1.238, -1.858, -2.355, -3.057, -3.595, -4.350, -4.880 ]
	        c2 = 1.1 * ones(13)
	        c3 = [ 0.000, 0.006, 0.006, -0.004, -0.017, -0.028, -0.040, -0.050, -0.055, -0.065, -0.070, -0.080, -0.100 ]
	        c4 = [ -2.100, -2.128, -2.148, -2.080, -2.028, -1.990, -1.945, -1.865, -1.800, -1.725, -1.670, -1.610, -1.570 ]
	        c5 = -0.48451 * ones(13)
	        c6 = 0.524 * ones(13)
	        c7 = [ 0.0, -0.082, -0.041, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ]
	        ca = [ 1.39, 1.40, 1.41, 1.43, 1.45, 1.48, 1.50, 1.52, 1.53, 1.53, 1.53, 1.53, 1.53 ]
	        cb = -0.14 * ones(13)
	        cc = [ 0.38, 0.39, 0.40, 0.42, 0.44, 0.47, 0.49, 0.51, 0.52, 0.52, 0.52, 0.52, 0.52 ]
		end
		if FaultType == 1
			mechFactor = 1.2
		end
	else
		# soil site
		if M <= 6.5
			if FaultType == 0
	            c1 = -2.17 * ones(13)
	            c6 = [ 0.0, 0.4572, 0.6395, 0.9187, 0.9547, 0.9251, 0.8494, 0.7010, 0.5665, 0.3235, 0.1001, -0.2801, -0.6274 ]
	        else
	            c1 = -1.92 * ones(13)
	            c6 = [ 0.0, 0.4572, 0.6395, 0.9187, 0.9547, 0.9005, 0.8285, 0.6802, 0.5075, 0.2215, -0.0526, -0.4905, -0.8907 ]
	        end
	        c2 = ones(13)
	        c3 = 1.70 * ones(13)
	        c4 = 2.1863 * ones(13)
	        c5 = 0.32 * ones(13)
	        c7 = [ 0.0, 0.005, 0.005, -0.004, -0.014, -0.024, -0.033, -0.051, -0.065, -0.090, -0.108, -0.139, -0.160 ]
	        ca = [ 1.52, 1.54, 1.54, 1.565, 1.58, 1.595, 1.61, 1.635, 1.66, 1.69, 1.70, 1.71, 1.71 ]
	        cb = -0.16 * ones(13)
	        cc = zeros(13)
		else
			if FaultType == 0
	            c1 = -2.17 * ones(13)
	            c6 = [ 0.0, 0.4572, 0.6395, 0.9187, 0.9547, 0.9251, 0.8494, 0.7010, 0.5665, 0.3235, 0.1001, -0.2801, -0.6274 ]
	        else
	            c1 = -1.92 * ones(13)
	            c6 = [ 0.0, 0.4572, 0.6395, 0.9187, 0.9547, 0.9005, 0.8285, 0.6802, 0.5075, 0.2215, -0.0526, -0.4905, -0.8907 ]
	        end
	        c2 = ones(13)
	        c3 = 1.70 * ones(13)
	        c4 = 0.3825 * ones(13)
	        c5 = 0.5882 * ones(13)
	        c7 = [ 0.0, 0.005, 0.005, -0.004, -0.014, -0.024, -0.033, -0.051, -0.065, -0.090, -0.108, -0.139, -0.160 ]
	        ca = [ 1.52, 1.54, 1.54, 1.565, 1.58, 1.595, 1.61, 1.635, 1.66, 1.69, 1.70, 1.71, 1.71 ]
	        cb = -0.16 * ones(13)
	        cc = zeros(13)
		end
	end


    T = (T > 4.0) ? 4.0 : T
    tidHI = findfirst(Ti .>= T)
    tidLO = findlast(Ti .<= T)

    if tidHI == tidLO
        # we have a known period
        tid = tidHI

		if SiteType == 0
	        lnSa = c1[tid] + c2[tid]*M + c3[tid]*((8.5 - M)^2.5) + c4[tid]*log( R + exp( c5[tid] + c6[tid]*M )) + c7[tid]*log( R + 2.0 )
	    else
	        lnSa = c1[tid] + c2[tid]*M - c3[tid]*log( R + c4[tid]*exp( c5[tid]*M ) ) + c6[tid] + c7[tid]*((8.5 - M)^2.5)
	    end
		lnSa += log(mechFactor)
		Sa = exp(lnSa)

		τ = NaN
		ϕ = NaN
		if SiteType == 0
			if M >= 7.21
				σ = cc[tid]
			else
				σ = ca[tid] + cb[tid]*M
			end
	    else
			σ = ca[tid] + cb[tid]*min(M, 7.0)
	    end

    else
        # we need to interpolate
        @inbounds TiLO = Ti[tidLO]
        @inbounds TiHI = Ti[tidHI]

        if T < 0.07
            # assume that PGA coefficients are really at 0.001s
            Tfrac = log(max(T,0.001)/0.001)/log(TiHI/0.001)
        else
            Tfrac = log(T/TiLO)/log(TiHI/TiLO)
        end

        sLO = PJSsadigh1997( TiLO, M, R, FaultType, SiteType )
        sHI = PJSsadigh1997( TiHI, M, R, FaultType, SiteType )

        lnSa = sLO.lnIM + Tfrac*(sHI.lnIM-sLO.lnIM)
        Sa = exp(lnSa)

        τ = sLO.τ + Tfrac*(sHI.τ - sLO.τ)
		ϕ = sLO.ϕ + Tfrac*(sHI.ϕ - sLO.ϕ)
		σ = sLO.σ + Tfrac*(sHI.σ - sLO.σ)
    end

    gm = PJSgroundMotion(Sa, lnSa, τ, ϕ, σ)
	return gm
end

@time s = PJSsadigh1997(0.1, 7.0, 10.0, 0, 0)


function PJSsadigh1997(Ti::Vector{U}, M::U, R::U, FaultType::Int64, SiteType::Int64) where U<:Real
	predictions = Array{PJSgroundMotion,1}()
	for T in Ti
		push!(predictions, PJSsadigh1997(T, M, R, FaultType, SiteType))
	end
	return predictions
end


function PJSsadigh1997(T::U, rup::Rupture{U}, site::Site{U}) where U<:Real
	R = rupture_distance(rup, site)
	# model uses generic site conditions. Rock sites are described as those having bedrock within about a metre of the surface
	# paper also states that deep soil sites are classed as those with greater than 20m of soil over bedrock (here bedrock is not hard rock)
	if depth_to_velocity_horizon(site, 750.0) > 20.0
		SiteType = 1
	else
		SiteType = 0
	end
	Fss, Fnm, Frv, Fuk = mechanism(rup)
	if Frv == 1
		FaultType = 1
	else
		FaultType = 0
	end
	gm = PJSsadigh1997(T, rup.magnitude, R, FaultType, SiteType)
	return gm
end


function PJSsadigh1997(Ti::Vector{U}, rup::Rupture{U}, site::Site{U}) where U<:Real
	predictions = Array{PJSgroundMotion,1}()
	for T in Ti
		push!(predictions, PJSsadigh1997(T, rup, site))
	end
	return predictions
end
