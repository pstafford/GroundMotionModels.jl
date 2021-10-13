
include("PJSgroundMotions.jl")


function PJScb2014_pga(M::T, Rrup::T, Rjb::T, Rx::T, Ztor::T, Zhyp::T, W::T, Dip::T, Fnm::Int64, Frv::Int64, Vs30::T, Z2p5::T) where T<:Real
	# define the coefficients for the PGA case
	c0 = -4.416
	c1 = 0.984
	c2 = 0.537
	c3 = -1.499
	c4 = -0.496
	c5 = -2.773
	c6 = 0.248
	c7 = 6.768
	c8 = 0.0
	c9 = -0.212
	c10 = 0.72
	c11 = 1.09
	c12 = 2.186
	c13 = 1.42
	c14 = -0.0064
	c15 = -0.202
	c16 = 0.393
	c17 = 0.0977
	c18 = 0.0333
	c19 = 0.00757
	c20 = -0.0055
	Δc20_JI = -0.0035
	Δc20_CH = 0.0036
	k1 = 865.0
	k2 = -1.186
	k3 = 1.839
	a2 = 0.167
	h1 = 0.241
	h2 = 1.474
	h3 = -0.715
	h5 = -0.337
	h6 = -0.27
	τ1 = 0.409
	τ2 = 0.322
	ϕ1 = 0.734
	ϕ2 = 0.492
	ϕlnAF = 0.3
	σ1 = 0.84
	σ2 = 0.588
	ρ = 1.0

	cc = 1.88
	n = 1.18
	h4 = 1.0

	# magnitude term
	if M <= 4.5
		fmag = c0 + c1*M
	elseif M <= 5.5
		fmag = c0 + c1*M + c2*( M - 4.5 )
	elseif M <= 6.5
		fmag = c0 + c1*M + c2*( M - 4.5 ) + c3*( M - 5.5 )
	else
		fmag = c0 + c1*M + c2*( M - 4.5 ) + c3*( M - 5.5 ) + c4*( M - 6.5 )
	end

	# geometric attenuation term
	fdis = ( c5 + c6*M ) * log( sqrt( Rrup^2 + c7^2 ) )

	# style-of-faulting term
	fflt_F = c8*Frv + c9*Fnm
	if fflt_F == 0
		if M <= 4.5
			fflt_M = 0.0
		elseif M <= 5.5
			fflt_M = M - 4.5
		else
			fflt_M = 1.0
		end
		fflt = fflt_F * fflt_M
	else
		fflt = 0.0
	end

	# hanging-wall term
	if isnan(Rx)
		fhng_Rx = 0.0
		fhng = 0.0
	else
		if Rx < 0
			fhng_Rx = 0.0
			fhng = 0.0
		else
			R1 = W * cosd( Dip )
			if Rx < R1
				f1_Rx = h1 + h2*(Rx/R1) + h3*(Rx/R1)^2
				fhng_Rx = f1_Rx
			else
				R2 = 62.0*M - 350.0
				f2_Rx = h4 + h5*((Rx-R1)/(R2-R1)) + h6*((Rx-R1)/(R2-R1))^2
				fhng_Rx = max( f2_Rx, 0.0 )
			end

			if Rrup == 0.0
				fhng_Rrup = 1.0
			else
				fhng_Rrup = ( Rrup - Rjb ) / Rrup
			end

			if M <= 5.5
				fhng_M = 0.0
			elseif M <= 6.5
				fhng_M = ( M - 5.5 )*( 1.0 + a2*( M - 6.5 ) )
			else
				fhng_M = 1.0 + a2*( M - 6.5 )
			end

			if Ztor <= 16.66
				fhng_Z = 1.0 - 0.06*Ztor
			else
				fhng_Z = 0.0
			end

			fhng_dip = ( 90.0 - Dip ) / 45.0

			fhng = c10 * fhng_Rx * fhng_Rrup * fhng_M * fhng_Z * fhng_dip
		end
	end


	# shallow site response
	# note that this ignores any regional difference with Japan
	if Vs30 > k1
		fsite_G = ( c11 + k2*n ) * log( Vs30 / k1 )
	else
		# must compute pga1100
		pga1100 = PJScb2014_pga( M, Rrup, Rjb, Rx, Ztor, Zhyp, W, Dip, Fnm, Frv, 1100.0, Z2p5 )
		fsite_G = c11 * log( Vs30 / k1 ) + k2*( log( pga1100.IM + cc*(Vs30/k1)^n ) - log( pga1100.IM + cc ) )
	end
	fsite = fsite_G

	# basin response term
	# ignoring Japanese regionalisation
	if Z2p5 <= 1.0
		fsed = c14 * ( Z2p5 - 1.0 )
	elseif Z2p5 <= 3.0
		fsed = 0.0
	else
		fsed = c16*k3*exp(-0.75)*( 1.0 - exp( -0.25*(Z2p5 - 3.0) ) )
	end

	# hypocentral depth term
	if Zhyp <= 7.0
		fhyp_H = 0.0
	elseif Zhyp <= 20.0
		fhyp_H = Zhyp - 7.0
	else
		fhyp_H = 13.0
	end
	if M <= 5.5
		fhyp_M = c17
	elseif M <= 6.5
		fhyp_M = c17 + (c18-c17)*(M - 5.5)
	else
		fhyp_M = c18
	end
	fhyp = fhyp_H * fhyp_M

	# fault dip term
	if M <= 4.5
		fdip = c19*Dip
	elseif M <= 5.5
		fdip = c19*(5.5-M)*Dip
	else
		fdip = 0.0
	end

	# anelastic attenuation term
	# ignoring regionalisation
	if Rrup > 80.0
		fatn = c20 * ( Rrup - 80.0 )
	else
		fatn = 0
	end

	lnPGA = fmag + fdis + fflt + fhng + fsite + fsed + fhyp + fdip + fatn
	pga = exp( lnPGA )

    # aleatory variability
	if M <= 4.5
	    τlnY = τ1
	    ϕlnY = ϕ1
	elseif M <= 5.5
	    τlnY = τ2 + ( τ1 - τ2 )*( 5.5 - M )
	    ϕlnY = ϕ2 + ( ϕ1 - ϕ2 )*( 5.5 - M )
	else
	    τlnY = τ2
	    ϕlnY = ϕ2
	end

	if Vs30 >= k1
	    α = 0.0
	else
	    α = k2*pga1100.IM * ( 1.0/( pga1100.IM + cc*(Vs30/k1)^n ) - 1.0/( pga1100.IM + cc ) )
	end

	if α == 0.0
		τ = τlnY
		ϕ = ϕlnY
	else
		τ = sqrt( τlnY^2 + α^2*pga1100.τ^2 + 2*α*ρ*τlnY*pga1100.τ )
		ϕlnYb = sqrt( ϕlnY^2 - ϕlnAF^2 )
		ϕ = sqrt( ϕlnYb^2 + ϕlnAF^2 + α^2*pga1100.ϕ^2 + 2*α*ρ*ϕlnYb*pga1100.ϕ )
	end

	σ = sqrt( τ^2 + ϕ^2 )

	gm = PJSgroundMotion(pga, lnPGA, τ, ϕ, σ)
	return gm
end

# M = 6.0
# Rrup = 10.0
# Rjb = 10.0
# Rx = 10.0
# Ztor = 0.0
# Zhyp = 10.0
# W = 12.0
# Dip = 90.0
# Fnm = 0
# Frv = 0
# Vs30 = 760.0
# Z2p5 = 1.5
#
# gg = PJSgmpe(0.1, log(0.1), 0.3, 0.4, 0.5)
# gg.Sa
# gg.lnSa
#
# @time cb_pga = PJScb2014_pga(M, Rrup, Rjb, Rx, Ztor, Zhyp, W, Dip, Fnm, Frv, Vs30, Z2p5)


function PJScb2014( T::U, M::U, Rrup::U, Rjb::U, Rx::U, Ztor::U, Zhyp::U, W::U, Dip::U, Fnm::Int64, Frv::Int64, Vs30::U, Z2p5::U ) where U<:Real
	# coefficients
	Ti = [ 0.01, 0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 7.5, 10.0 ]
	c0 = [ -4.365, -4.348, -4.024, -3.479, -3.293, -3.666, -4.866, -5.411, -5.962, -6.403, -7.566, -8.379, -9.841, -11.011, -12.469, -12.969, -13.306, -14.02, -14.558, -15.509, -15.975 ]
	c1 = [ 0.977, 0.976, 0.931, 0.887, 0.902, 0.993, 1.267, 1.366, 1.458, 1.528, 1.739, 1.872, 2.021, 2.18, 2.27, 2.271, 2.15, 2.132, 2.116, 2.223, 2.132 ]
	c2 = [ 0.533, 0.549, 0.628, 0.674, 0.726, 0.698, 0.51, 0.447, 0.274, 0.193, -0.02, -0.121, -0.042, -0.069, 0.047, 0.149, 0.368, 0.726, 1.027, 0.169, 0.367 ]
	c3 = [ -1.485, -1.488, -1.494, -1.388, -1.469, -1.572, -1.669, -1.75, -1.711, -1.77, -1.594, -1.577, -1.757, -1.707, -1.621, -1.512, -1.315, -1.506, -1.721, -0.756, -0.8 ]
	c4 = [ -0.499, -0.501, -0.517, -0.615, -0.596, -0.536, -0.49, -0.451, -0.404, -0.321, -0.426, -0.44, -0.443, -0.527, -0.63, -0.768, -0.89, -0.885, -0.878, -1.077, -1.282 ]
	c5 = [ -2.773, -2.772, -2.782, -2.791, -2.745, -2.633, -2.458, -2.421, -2.392, -2.376, -2.303, -2.296, -2.232, -2.158, -2.063, -2.104, -2.051, -1.986, -2.021, -2.179, -2.244 ]
	c6 = [ 0.248, 0.247, 0.246, 0.24, 0.227, 0.21, 0.183, 0.182, 0.189, 0.195, 0.185, 0.186, 0.186, 0.169, 0.158, 0.158, 0.148, 0.135, 0.14, 0.178, 0.194 ]
	c7 = [ 6.753, 6.502, 6.291, 6.317, 6.861, 7.294, 8.031, 8.385, 7.534, 6.99, 7.012, 6.902, 5.522, 5.65, 5.795, 6.632, 6.759, 7.978, 8.538, 8.468, 6.564 ]
	c8 = [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ]
	c9 = [ -0.214, -0.208, -0.213, -0.244, -0.266, -0.229, -0.211, -0.163, -0.15, -0.131, -0.159, -0.153, -0.09, -0.105, -0.058, -0.028, 0.0, 0.0, 0.0, 0.0, 0.0 ]
	c10 = [ 0.72, 0.73, 0.759, 0.826, 0.815, 0.831, 0.749, 0.764, 0.716, 0.737, 0.738, 0.718, 0.795, 0.556, 0.48, 0.401, 0.206, 0.105, 0.0, 0.0, 0.0 ]
	c11 = [ 1.094, 1.149, 1.29, 1.449, 1.535, 1.615, 1.877, 2.069, 2.205, 2.306, 2.398, 2.355, 1.995, 1.447, 0.33, -0.514, -0.848, -0.793, -0.748, -0.664, -0.576 ]
	c12 = [ 2.191, 2.189, 2.164, 2.138, 2.446, 2.969, 3.544, 3.707, 3.343, 3.334, 3.544, 3.016, 2.616, 2.47, 2.108, 1.327, 0.601, 0.568, 0.356, 0.075, -0.027 ]
	c13 = [ 1.416, 1.453, 1.476, 1.549, 1.772, 1.916, 2.161, 2.465, 2.766, 3.011, 3.203, 3.333, 3.054, 2.562, 1.453, 0.657, 0.367, 0.306, 0.268, 0.374, 0.297 ]
	c14 = [ -0.007, -0.0167, -0.0422, -0.0663, -0.0794, -0.0294, 0.0642, 0.0968, 0.1441, 0.1597, 0.141, 0.1474, 0.1764, 0.2593, 0.2881, 0.3112, 0.3478, 0.3747, 0.3382, 0.3754, 0.3506 ]
	c15 = [ -0.207, -0.199, -0.202, -0.339, -0.404, -0.416, -0.407, -0.311, -0.172, -0.084, 0.085, 0.233, 0.411, 0.479, 0.566, 0.562, 0.534, 0.522, 0.477, 0.321, 0.174 ]
	c16 = [ 0.39, 0.387, 0.378, 0.295, 0.322, 0.384, 0.417, 0.404, 0.466, 0.528, 0.54, 0.638, 0.776, 0.771, 0.748, 0.763, 0.686, 0.691, 0.67, 0.757, 0.621 ]
	c17 = [ 0.0981, 0.1009, 0.1095, 0.1226, 0.1165, 0.0998, 0.076, 0.0571, 0.0437, 0.0323, 0.0209, 0.0092, -0.0082, -0.0131, -0.0187, -0.0258, -0.0311, -0.0413, -0.0281, -0.0205, 0.0009 ]
	c18 = [ 0.0334, 0.0327, 0.0331, 0.027, 0.0288, 0.0325, 0.0388, 0.0437, 0.0463, 0.0508, 0.0432, 0.0405, 0.042, 0.0426, 0.038, 0.0252, 0.0236, 0.0102, 0.0034, 0.005, 0.0099 ]
	c19 = [ 0.00755, 0.00759, 0.0079, 0.00803, 0.00811, 0.00744, 0.00716, 0.00688, 0.00556, 0.00458, 0.00401, 0.00388, 0.0042, 0.00409, 0.00424, 0.00448, 0.00345, 0.00603, 0.00805, 0.0028, 0.00458 ]
	c20 = [ -0.0055, -0.0055, -0.0057, -0.0063, -0.007, -0.0073, -0.0069, -0.006, -0.0055, -0.0049, -0.0037, -0.0027, -0.0016, -0.0006, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ]
	Δc20_JI = [ -0.0035, -0.0035, -0.0034, -0.0037, -0.0037, -0.0034, -0.003, -0.0031, -0.0033, -0.0035, -0.0034, -0.0034, -0.0032, -0.003, -0.0019, -0.0005, 0.0, 0.0, 0.0, 0.0, 0.0 ]
	Δc20_CH = [ 0.0036, 0.0036, 0.0037, 0.004, 0.0039, 0.0042, 0.0042, 0.0041, 0.0036, 0.0031, 0.0028, 0.0025, 0.0016, 0.0006, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ]
	k1 = [ 865.0, 865.0, 908.0, 1054.0, 1086.0, 1032.0, 878.0, 748.0, 654.0, 587.0, 503.0, 457.0, 410.0, 400.0, 400.0, 400.0, 400.0, 400.0, 400.0, 400.0, 400.0 ]
	k2 = [ -1.186, -1.219, -1.273, -1.346, -1.471, -1.624, -1.931, -2.188, -2.381, -2.518, -2.657, -2.669, -2.401, -1.955, -1.025, -0.299, 0.0, 0.0, 0.0, 0.0, 0.0 ]
	k3 = [ 1.839, 1.84, 1.841, 1.843, 1.845, 1.847, 1.852, 1.856, 1.861, 1.865, 1.874, 1.883, 1.906, 1.929, 1.974, 2.019, 2.11, 2.2, 2.291, 2.517, 2.744 ]
	a2 = [ 0.168, 0.166, 0.167, 0.173, 0.198, 0.174, 0.198, 0.204, 0.185, 0.164, 0.16, 0.184, 0.216, 0.596, 0.596, 0.596, 0.596, 0.596, 0.596, 0.596, 0.596 ]
	h1 = [ 0.242, 0.244, 0.246, 0.251, 0.26, 0.259, 0.254, 0.237, 0.206, 0.21, 0.226, 0.217, 0.154, 0.117, 0.117, 0.117, 0.117, 0.117, 0.117, 0.117, 0.117 ]
	h2 = [ 1.471, 1.467, 1.467, 1.449, 1.435, 1.449, 1.461, 1.484, 1.581, 1.586, 1.544, 1.554, 1.626, 1.616, 1.616, 1.616, 1.616, 1.616, 1.616, 1.616, 1.616 ]
	h3 = [ -0.714, -0.711, -0.713, -0.701, -0.695, -0.708, -0.715, -0.721, -0.787, -0.795, -0.77, -0.77, -0.78, -0.733, -0.733, -0.733, -0.733, -0.733, -0.733, -0.733, -0.733 ]
	h5 = [ -0.336, -0.339, -0.338, -0.338, -0.347, -0.391, -0.449, -0.393, -0.339, -0.447, -0.525, -0.407, -0.371, -0.128, -0.128, -0.128, -0.128, -0.128, -0.128, -0.128, -0.128 ]
	h6 = [ -0.27, -0.263, -0.259, -0.263, -0.219, -0.201, -0.099, -0.198, -0.21, -0.121, -0.086, -0.281, -0.285, -0.756, -0.756, -0.756, -0.756, -0.756, -0.756, -0.756, -0.756 ]
	τ1 = [ 0.404, 0.417, 0.446, 0.508, 0.504, 0.445, 0.382, 0.339, 0.34, 0.34, 0.356, 0.379, 0.43, 0.47, 0.497, 0.499, 0.5, 0.543, 0.534, 0.523, 0.466 ]
	τ2 = [ 0.325, 0.326, 0.344, 0.377, 0.418, 0.426, 0.387, 0.338, 0.316, 0.3, 0.264, 0.263, 0.326, 0.353, 0.399, 0.4, 0.417, 0.393, 0.421, 0.438, 0.438 ]
	ϕ1 = [ 0.734, 0.738, 0.747, 0.777, 0.782, 0.769, 0.769, 0.761, 0.744, 0.727, 0.69, 0.663, 0.606, 0.579, 0.541, 0.529, 0.527, 0.521, 0.502, 0.457, 0.441 ]
	ϕ2 = [ 0.492, 0.496, 0.503, 0.52, 0.535, 0.543, 0.543, 0.552, 0.545, 0.568, 0.593, 0.611, 0.633, 0.628, 0.603, 0.588, 0.578, 0.559, 0.551, 0.546, 0.543 ]
	ϕlnAF = [ 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3 ]
	σ1 = [ 0.838, 0.848, 0.87, 0.928, 0.93, 0.888, 0.859, 0.833, 0.818, 0.803, 0.776, 0.764, 0.743, 0.746, 0.735, 0.727, 0.726, 0.753, 0.733, 0.695, 0.642 ]
	σ2 = [ 0.59, 0.594, 0.609, 0.642, 0.679, 0.69, 0.667, 0.647, 0.63, 0.642, 0.649, 0.665, 0.712, 0.72, 0.723, 0.711, 0.713, 0.683, 0.693, 0.7, 0.698 ]
	ρ = [ 1.0, 0.998, 0.986, 0.938, 0.887, 0.87, 0.876, 0.87, 0.85, 0.819, 0.743, 0.684, 0.562, 0.467, 0.364, 0.298, 0.234, 0.202, 0.184, 0.176, 0.154 ]

	cc = 1.88
	n = 1.18
	h4 = 1.0

	if T == 0.0 # PGA requested

		cb_pga = PJScb2014_pga( M, Rrup, Rjb, Rx, Ztor, Zhyp, W, Dip, Fnm, Frv, Vs30, Z2p5 )

		return cb_pga
	else
		# limit the period range (note that coefficients for PGA and 0.01s are not equal, but don't interpolate below 0.01s)
		T = ( T < 0.01 ) ? 0.01 : ( T > 10.0 ) ? 10.0 : T
		tidHI = findfirst(Ti .>= T)
		tidLO = findlast(Ti .<= T)


		# if T < 0.25 we are supposed to compute PGA regardless as it may replace the Sa value if its larger
		if T < 0.25
			cb_pga = PJScb2014_pga( M, Rrup, Rjb, Rx, Ztor, Zhyp, W, Dip, Fnm, Frv, Vs30, Z2p5 ).IM
		else
			cb_pga = NaN
		end

		if tidLO == tidHI
			# we have a known period
			tid = tidHI

			# magnitude term
			if M <= 4.5
				fmag = c0[tid] + c1[tid]*M
			elseif M <= 5.5
				fmag = c0[tid] + c1[tid]*M + c2[tid]*( M - 4.5 )
			elseif M <= 6.5
				fmag = c0[tid] + c1[tid]*M + c2[tid]*( M - 4.5 ) + c3[tid]*( M - 5.5 )
			else
				fmag = c0[tid] + c1[tid]*M + c2[tid]*( M - 4.5 ) + c3[tid]*( M - 5.5 ) + c4[tid]*( M - 6.5 )
			end

			# geometric attenuation term
			fdis = ( c5[tid] + c6[tid]*M ) * log( sqrt( Rrup^2 + c7[tid]^2 ) )

			# style-of-faulting term
			fflt_F = c8[tid]*Frv + c9[tid]*Fnm
			if fflt_F == 0.0
				if M <= 4.5
					fflt_M = 0.0
				elseif M <= 5.5
					fflt_M = M - 4.5
				else
					fflt_M = 1.0
				end
				fflt = fflt_F * fflt_M
			else
				fflt = 0.0
			end

			# hanging-wall term
			if isnan(Rx)
				fhng_Rx = 0.0
				fhng = 0.0
			else
				if Rx < 0.0
					fhng_Rx = 0.0
					fhng = 0.0
				else
					R1 = W * cosd( Dip )
					if Rx < R1
						f1_Rx = h1[tid] + h2[tid]*(Rx/R1) + h3[tid]*(Rx/R1)^2
						fhng_Rx = f1_Rx
					else
						R2 = 62.0*M - 350.0
						f2_Rx = h4 + h5[tid]*((Rx-R1)/(R2-R1)) + h6[tid]*((Rx-R1)/(R2-R1))^2
						fhng_Rx = max( f2_Rx, 0.0 )
					end

					if Rrup == 0.0
						fhng_Rrup = 1.0
					else
						fhng_Rrup = ( Rrup - Rjb ) / Rrup
					end

					if M <= 5.5
						fhng_M = 0.0
					elseif M <= 6.5
						fhng_M = ( M - 5.5 )*( 1.0 + a2[tid]*( M - 6.5 ) )
					else
						fhng_M = 1.0 + a2[tid]*( M - 6.5 )
					end

					if Ztor <= 16.66
						fhng_Z = 1.0 - 0.06*Ztor
					else
						fhng_Z = 0.0
					end

					fhng_dip = ( 90.0 - Dip ) / 45.0

					fhng = c10[tid] * fhng_Rx * fhng_Rrup * fhng_M * fhng_Z * fhng_dip
				end
			end

			# shallow site response
			# note that this ignores any regional difference with Japan
			if Vs30 > k1[tid]
				fsite_G = ( c11[tid] + k2[tid]*n ) * log( Vs30 / k1[tid] )
			else
				# must compute pga1100
				cb_pga1100 = PJScb2014_pga( M, Rrup, Rjb, Rx, Ztor, Zhyp, W, Dip, Fnm, Frv, 1100.0, Z2p5 )
				pga1100 = cb_pga1100.IM
				fsite_G = c11[tid] * log( Vs30 / k1[tid] ) + k2[tid]*( log( pga1100 + cc*(Vs30/k1[tid])^n ) - log( pga1100 + cc ) )
			end
			fsite = fsite_G

			# basin response term
			# ignoring Japanese regionalisation
			if Z2p5 <= 1.0
				fsed = c14[tid] * ( Z2p5 - 1.0 )
			elseif Z2p5 <= 3.0
				fsed = 0.0
			else
				fsed = c16[tid]*k3[tid]*exp(-0.75)*( 1.0 - exp( -0.25*(Z2p5 - 3.0) ) )
			end

			# hypocentral depth term
			if Zhyp <= 7.0
				fhyp_H = 0.0
			elseif Zhyp <= 20.0
				fhyp_H = Zhyp - 7.0
			else
				fhyp_H = 13.0
			end
			if M <= 5.5
				fhyp_M = c17[tid]
			elseif M <= 6.5
				fhyp_M = c17[tid] + (c18[tid]-c17[tid])*(M - 5.5)
			else
				fhyp_M = c18[tid]
			end
			fhyp = fhyp_H * fhyp_M

			# fault dip term
			if M <= 4.5
				fdip = c19[tid]*Dip
			elseif M <= 5.5
				fdip = c19[tid]*(5.5 - M)*Dip
			else
				fdip = 0.0
			end

			# anelastic attenuation term
			# ignoring regionalisation
			if Rrup > 80.0
				fatn = c20[tid] * ( Rrup - 80.0 )
			else
				fatn = 0.0
			end

			lnSa = fmag + fdis + fflt + fhng + fsite + fsed + fhyp + fdip + fatn
			Sa = exp( lnSa )

            if T < 0.25
                if Sa < cb_pga
                    lnSa = log(cb_pga)
                    Sa = cb_pga
                end
            end

			if M <= 4.5
				τlnY = τ1[tid]
				ϕlnY = ϕ1[tid]
			elseif M <= 5.5
				τlnY = τ2[tid] + ( τ1[tid] - τ2[tid] )*( 5.5 - M )
				ϕlnY = ϕ2[tid] + ( ϕ1[tid] - ϕ2[tid] )*( 5.5 - M )
			else
				τlnY = τ2[tid]
				ϕlnY = ϕ2[tid]
			end

			if Vs30 >= k1[tid]
				α = 0.0
			else
				α = k2[tid]*pga1100 * ( 1.0/( pga1100 + cc*(Vs30/k1[tid])^n ) - 1.0/( pga1100 + cc ) )
			end

			if α == 0.0
				τ = τlnY
				ϕ = ϕlnY
			else
				τ = sqrt( τlnY^2 + α^2*cb_pga1100.τ^2 + 2.0*α*ρ[tid]*τlnY*cb_pga1100.τ )
				ϕlnYb = sqrt( ϕlnY^2 - ϕlnAF[tid]^2 )
				ϕ = sqrt( ϕlnYb^2 + ϕlnAF[tid]^2 + α^2*cb_pga1100.ϕ^2 + 2.0*α*ρ[tid]*ϕlnYb*cb_pga1100.ϕ )
			end

			σ = sqrt( τ^2 + ϕ^2 )

		else
			# we need to interpolate
			cbLO = PJScb2014( Ti[tidLO], M, Rrup, Rjb, Rx, Ztor, Zhyp, W, Dip, Fnm, Frv, Vs30, Z2p5 )
			cbHI = PJScb2014( Ti[tidHI], M, Rrup, Rjb, Rx, Ztor, Zhyp, W, Dip, Fnm, Frv, Vs30, Z2p5 )

			lnSa = cbLO.lnIM + log(T/Ti[tidLO])/log(Ti[tidHI]/Ti[tidLO])*(cbHI.lnIM-cbLO.lnIM)
			Sa = exp( lnSa )

			τ = cbLO.τ + log(T/Ti[tidLO])/log(Ti[tidHI]/Ti[tidLO])*(cbHI.τ-cbLO.τ)
			ϕ = cbLO.ϕ + log(T/Ti[tidLO])/log(Ti[tidHI]/Ti[tidLO])*(cbHI.ϕ-cbLO.ϕ)
			σ = cbLO.σ + log(T/Ti[tidLO])/log(Ti[tidHI]/Ti[tidLO])*(cbHI.σ-cbLO.σ)
		end

		gm = PJSgroundMotion(Sa, lnSa, τ, ϕ, σ)
		return gm
	end
end


function PJScb2014( Ti::Vector{U}, M::U, Rrup::U, Rjb::U, Rx::U, Ztor::U, Zhyp::U, W::U, Dip::U, Fnm::Int64, Frv::Int64, Vs30::U, Z2p5::U=0.0 ) where U <: Real
	predictions = Array{PJSgroundMotion,1}()
	for T in Ti
		push!(predictions, PJScb2014( T, M, Rrup, Rjb, Rx, Ztor, Zhyp, W, Dip, Fnm, Frv, Vs30, Z2p5 ))
	end
	return predictions
end


function PJScb2014(T::U, rup::Rupture{U}, site::Site{U}) where U<:Real
	Rrup = rupture_distance(rup, site)
	Rjb = joyner_boore_distance(rup, site)
	Rx = strike_distance(rup, site)
	Ztor = ztor(rup)
	Zhyp = zhyp(rup)
	W = rup.rupWidth
	Dip = dip(rup)
	Vs30 = vs30(site)
	Z2p5 = z2p5(site) / 1000.0 # needs to be in units of km
	Fss, Fnm, Frv, Fuk = mechanism(rup)

	gm = PJScb2014(T, rup.magnitude, Rrup, Rjb, Rx, Ztor, Zhyp, W, Dip, Fnm, Frv, Vs30, Z2p5 )
	return gm
end


function PJScb2014(Ti::Vector{U}, rup::Rupture{U}, site::Site{U}) where U<:Real
	predictions = Array{PJSgroundMotion,1}()
	for T in Ti
		push!(predictions, PJScb2014(T, rup, site))
	end
	return predictions
end



# cb = PJScb2014( 0.01, M, Rrup, Rjb, Rx, Ztor, Zhyp, W, Dip, Fnm, Frv, Vs30, Z2p5)
# cb = PJScb2014( 0.02, M, Rrup, Rjb, Rx, Ztor, Zhyp, W, Dip, Fnm, Frv, Vs30, Z2p5)
#
#
# Ti = exp10.(range(-2.0, stop=1.0, step=0.1))
#
# cbi = PJScb2014( Ti, M, Rrup, Rjb, Rx, Ztor, Zhyp, W, Dip, Fnm, Frv, Vs30, Z2p5)
#
# Sai = get_spectrum( cbi )
#
# using Plots
# pyplot()
#
# plot(Ti, Sai, lab="", xaxis=:log10, yaxis=:log10)
