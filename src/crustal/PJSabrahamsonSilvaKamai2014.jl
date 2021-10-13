
include("PJSgroundMotions.jl")
# using Main.PJSgroundMotions

# using PJSgroundMotions


function PJSask2014( T::U, M::U, Rrup::U, Rjb::U, Rx::U, Ry0::U, Ztor::U, W::U, Dip::U, Fnm::Int64, Frv::Int64, Fhw::Int64, Vs30::U, Z1p0::U, Vs30known::Int64=1, region::String="California", Fas::Int64=0, CRjb::U=NaN ) where U <: Real

    # Period-dependent coefficients
    Ti = [ -1.0, 0.0, 0.01, 0.02, 0.03, 0.05, 0.075, 0.1, 0.15, 0.2, 0.25, 0.3, 0.4, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 7.5, 10.0 ]
    Vlin = [ 330.0, 660.0, 660.0, 680.0, 770.0, 915.0, 960.0, 910.0, 740.0, 590.0, 495.0, 430.0, 360.0, 340.0, 330.0, 330.0, 330.0, 330.0, 330.0, 330.0, 330.0, 330.0, 330.0, 330 ]
    b = [ -2.02, -1.47, -1.47, -1.459, -1.39, -1.219, -1.152, -1.23, -1.587, -2.012, -2.411, -2.757, -3.278, -3.599, -3.8, -3.5, -2.4, -1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ]
    n = [ 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5, 1.5 ]
    M1 = [ 6.75, 6.75, 6.75, 6.75, 6.75, 6.75, 6.75, 6.75, 6.75, 6.75, 6.75, 6.75, 6.75, 6.75, 6.75, 6.75, 6.75, 6.75, 6.82, 6.92, 7, 7.06, 7.145, 7.25 ]
    c = [ 2400, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4, 2.4 ]
    c4 = [ 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5, 4.5 ]
    a1 = [ 5.975, 0.587, 0.587, 0.598, 0.602, 0.707, 0.973, 1.169, 1.442, 1.637, 1.701, 1.712, 1.662, 1.571, 1.299, 1.043, 0.665, 0.329, -0.06, -0.299, -0.562, -0.875, -1.303, -1.928 ]
    a2 = [ -0.919, -0.79, -0.79, -0.79, -0.79, -0.79, -0.79, -0.79, -0.79, -0.79, -0.79, -0.79, -0.79, -0.79, -0.79, -0.79, -0.79, -0.79, -0.79, -0.79, -0.765, -0.711, -0.634, -0.529 ]
    a3 = [ 0.275, 0.275, 0.275, 0.275, 0.275, 0.275, 0.275, 0.275, 0.275, 0.275, 0.275, 0.275, 0.275, 0.275, 0.275, 0.275, 0.275, 0.275, 0.275, 0.275, 0.275, 0.275, 0.275, 0.275 ]
    a4 = [ -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1 ]
    a5 = [ -0.41, -0.41, -0.41, -0.41, -0.41, -0.41, -0.41, -0.41, -0.41, -0.41, -0.41, -0.41, -0.41, -0.41, -0.41, -0.41, -0.41, -0.41, -0.41, -0.41, -0.41, -0.41, -0.41, -0.41 ]
    a6 = [ 2.3657, 2.1541, 2.1541, 2.1461, 2.1566, 2.0845, 2.0285, 2.0408, 2.1208, 2.2241, 2.3124, 2.3383, 2.4688, 2.5586, 2.6821, 2.763, 2.8355, 2.8973, 2.9061, 2.8888, 2.8984, 2.8955, 2.87, 2.8431 ]
    a7 = [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ]
    a8 = [ -0.094, -0.015, -0.015, -0.015, -0.015, -0.015, -0.015, -0.015, -0.022, -0.03, -0.038, -0.045, -0.055, -0.065, -0.095, -0.11, -0.124, -0.138, -0.172, -0.197, -0.218, -0.235, -0.255, -0.285 ]
    a10 = [ 2.36, 1.735, 1.735, 1.718, 1.615, 1.358, 1.258, 1.31, 1.66, 2.22, 2.77, 3.25, 3.99, 4.45, 4.75, 4.3, 2.6, 0.55, -0.95, -0.95, -0.93, -0.91, -0.87, -0.8 ]
    a11 = [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ]
    a12 = [ -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.1, -0.2, -0.2, -0.2 ]
    a13 = [ 0.25, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.6, 0.58, 0.56, 0.53, 0.5, 0.42, 0.35, 0.2, 0.0, 0.0, 0.0, 0.0, 0.0 ]
    a14 = [ 0.22, -0.3, -0.3, -0.3, -0.3, -0.3, -0.3, -0.3, -0.3, -0.3, -0.24, -0.19, -0.11, -0.04, 0.07, 0.15, 0.27, 0.35, 0.46, 0.54, 0.61, 0.65, 0.72, 0.8 ]
    a15 = [ 0.3, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.1, 1.03, 0.92, 0.84, 0.68, 0.57, 0.42, 0.31, 0.16, 0.05, -0.04, -0.11, -0.19, -0.3 ]
    a16 = [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ]
    a17 = [ -0.0005, -0.0072, -0.0072, -0.0073, -0.0075, -0.008, -0.0089, -0.0095, -0.0095, -0.0086, -0.0074, -0.0064, -0.0043, -0.0032, -0.0025, -0.0025, -0.0022, -0.0019, -0.0015, -0.001, -0.001, -0.001, -0.001, -0.001 ]
    a43 = [ 0.28, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.14, 0.17, 0.22, 0.26, 0.34, 0.41, 0.51, 0.55, 0.49, 0.42 ]
    a44 = [ 0.15, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.07, 0.1, 0.14, 0.17, 0.21, 0.25, 0.3, 0.32, 0.32, 0.32, 0.275, 0.22 ]
    a45 = [ 0.09, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.03, 0.06, 0.1, 0.14, 0.17, 0.2, 0.22, 0.23, 0.23, 0.22, 0.2, 0.17, 0.14 ]
    a46 = [ 0.07, -0.05, -0.05, -0.05, -0.05, -0.05, -0.05, -0.05, -0.05, -0.03, 0.0, 0.03, 0.06, 0.09, 0.13, 0.14, 0.16, 0.16, 0.16, 0.14, 0.13, 0.1, 0.09, 0.08 ]
    a25 = [ -0.0001, -0.0015, -0.0015, -0.0015, -0.0016, -0.002, -0.0027, -0.0033, -0.0035, -0.0033, -0.0029, -0.0027, -0.0023, -0.002, -0.001, -0.0005, -0.0004, -0.0002, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ]
    a28 = [ 0.0005, 0.0025, 0.0025, 0.0024, 0.0023, 0.0027, 0.0032, 0.0036, 0.0033, 0.0027, 0.0024, 0.002, 0.001, 0.0008, 0.0007, 0.0007, 0.0006, 0.0003, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ]
    a29 = [ -0.0037, -0.0034, -0.0034, -0.0033, -0.0034, -0.0033, -0.0029, -0.0025, -0.0025, -0.0031, -0.0036, -0.0039, -0.0048, -0.005, -0.0041, -0.0032, -0.002, -0.0017, -0.002, -0.002, -0.002, -0.002, -0.002, -0.002 ]
    a31 = [ -0.1462, -0.1503, -0.1503, -0.1479, -0.1447, -0.1326, -0.1353, -0.1128, 0.0383, 0.0775, 0.0741, 0.2548, 0.2136, 0.1542, 0.0787, 0.0476, -0.0163, -0.1203, -0.2719, -0.2958, -0.2718, -0.2517, -0.14, -0.0216 ]
    a36 = [ 0.377, 0.265, 0.265, 0.255, 0.249, 0.202, 0.126, 0.022, -0.136, -0.078, 0.037, -0.091, 0.129, 0.31, 0.505, 0.358, 0.131, 0.123, 0.109, 0.135, 0.189, 0.215, 0.15, 0.092 ]
    a37 = [ 0.212, 0.337, 0.337, 0.328, 0.32, 0.289, 0.275, 0.256, 0.162, 0.224, 0.248, 0.203, 0.232, 0.252, 0.208, 0.208, 0.108, 0.068, -0.023, 0.028, 0.031, 0.024, -0.07, -0.159 ]
    a38 = [ 0.157, 0.188, 0.188, 0.184, 0.18, 0.167, 0.173, 0.189, 0.108, 0.115, 0.122, 0.096, 0.123, 0.134, 0.129, 0.152, 0.118, 0.119, 0.093, 0.084, 0.058, 0.065, 0.0, -0.05 ]
    a39 = [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ]
    a40 = [ 0.095, 0.088, 0.088, 0.088, 0.093, 0.133, 0.186, 0.16, 0.068, 0.048, 0.055, 0.073, 0.143, 0.16, 0.158, 0.145, 0.131, 0.083, 0.07, 0.101, 0.095, 0.133, 0.151, 0.124 ]
    a41 = [ -0.038, -0.196, -0.196, -0.194, -0.175, -0.09, 0.09, 0.006, -0.156, -0.274, -0.248, -0.203, -0.154, -0.159, -0.141, -0.144, -0.126, -0.075, -0.021, 0.072, 0.205, 0.285, 0.329, 0.301 ]
    a42 = [ 0.065, 0.044, 0.044, 0.061, 0.162, 0.451, 0.506, 0.335, -0.084, -0.178, -0.187, -0.159, -0.023, -0.029, 0.061, 0.062, 0.037, -0.143, -0.028, -0.097, 0.015, 0.104, 0.299, 0.243 ]
    s1_VsU = [ 0.662, 0.754, 0.754, 0.76, 0.781, 0.81, 0.81, 0.81, 0.801, 0.789, 0.77, 0.74, 0.699, 0.676, 0.631, 0.609, 0.578, 0.555, 0.548, 0.527, 0.505, 0.477, 0.457, 0.429 ]
    s2_VsU = [ 0.51, 0.52, 0.52, 0.52, 0.52, 0.53, 0.54, 0.55, 0.56, 0.565, 0.57, 0.58, 0.59, 0.6, 0.615, 0.63, 0.64, 0.65, 0.64, 0.63, 0.63, 0.63, 0.63, 0.63 ]
    s3 = [ 0.38, 0.47, 0.47, 0.47, 0.47, 0.47, 0.47, 0.47, 0.47, 0.47, 0.47, 0.47, 0.47, 0.47, 0.47, 0.47, 0.47, 0.47, 0.47, 0.47, 0.47, 0.47, 0.47, 0.47 ]
    s4 = [ 0.38, 0.36, 0.36, 0.36, 0.36, 0.36, 0.36, 0.36, 0.36, 0.36, 0.36, 0.36, 0.36, 0.36, 0.36, 0.36, 0.36, 0.36, 0.36, 0.36, 0.36, 0.36, 0.36, 0.36 ]
    s1_VsK = [ 0.66, 0.741, 0.741, 0.747, 0.769, 0.798, 0.798, 0.795, 0.773, 0.753, 0.729, 0.693, 0.644, 0.616, 0.566, 0.541, 0.506, 0.48, 0.472, 0.447, 0.425, 0.395, 0.378, 0.359 ]
    s2_VsK = [ 0.51, 0.501, 0.501, 0.501, 0.501, 0.512, 0.522, 0.527, 0.519, 0.514, 0.513, 0.519, 0.524, 0.532, 0.548, 0.565, 0.576, 0.587, 0.576, 0.565, 0.568, 0.571, 0.575, 0.585 ]
    s5_JP = [ 0.58, 0.54, 0.54, 0.54, 0.55, 0.56, 0.57, 0.57, 0.58, 0.59, 0.61, 0.63, 0.66, 0.69, 0.73, 0.77, 0.8, 0.8, 0.8, 0.76, 0.72, 0.7, 0.67, 0.64 ]
    s6_JP = [ 0.53, 0.63, 0.63, 0.63, 0.63, 0.65, 0.69, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.7, 0.69, 0.68, 0.66, 0.62, 0.55, 0.52, 0.5, 0.5, 0.5, 0.5 ]

	# period-independent values
	M2 = 5.0

    # check bounds on the allowable value of T
	if T < 0.0 # PGV requested
		tidHI = tidLO = 1
	else
		# limit the period range (note that coefficients for 0.01 and 0.0 are equal)
		T = ( T < 0.01 ) ? 0.01 : ( T > 10.0 ) ? 10.0 : T
		tidHI = findfirst(Ti .>= T)
		tidLO = findlast(Ti .<= T)
	end

	if ( tidLO == tidHI )
		# we have a known period
		tid = tidHI

		# near-source saturation
		if M <= 4.0
			c4m = 1.0
		elseif M <= 5.0
			@inbounds c4m = c4[tid] - (c4[tid]-1.0)*(5.0 - M)
		else
			@inbounds c4m = c4[tid]
		end

		# effective distance
		R = sqrt( Rrup^2 + c4m^2 )

		# base scaling
		if M < M2
			@inbounds f1 = a1[tid] + a4[tid]*(M2 - M1[tid]) + a8[tid]*(8.5 - M2)^2 + a6[tid]*(M - M2) + a7[tid]*(M - M2)^2 + (a2[tid] + a3[tid]*(M2 - M1[tid]))*log(R) + a17[tid]*Rrup
		elseif M < M1[tid]
			@inbounds f1 = a1[tid] + a4[tid]*(M - M1[tid]) + a8[tid]*(8.5 - M)^2 + (a2[tid] + a3[tid]*(M - M1[tid]))*log(R) + a17[tid]*Rrup
		else
			@inbounds f1 = a1[tid] + a5[tid]*(M - M1[tid]) + a8[tid]*(8.5 - M)^2 + (a2[tid] + a3[tid]*(M - M1[tid]))*log(R) + a17[tid]*Rrup
		end

		# style of faulting terms
		if Frv == 1
			if M < 4.0
				f7 = 0.0
			elseif M <= 5.0
				@inbounds f7 = a11[tid]*(M - 4.0)
			else
				@inbounds f7 = a11[tid]
			end
		else
			f7 = 0.0
		end
		if Fnm == 1
			if M < 4.0
				f8 = 0.0
			elseif M <= 5.0
				@inbounds f8 = a12[tid]*(M - 4.0)
			else
				@inbounds f8 = a12[tid]
			end
		else
			f8 = 0.0
		end

		# site response model
		if T <= 0.5
			V1 = 1500.0
		elseif T < 3.0
			V1 = exp( -0.35*log(T/0.5) + log(1500.0) )
		else
			V1 = 800.0
		end
		Vs30_star = (Vs30 < V1) ? Vs30 : V1

		if Vs30 >= Vlin[tid]
			@inbounds f5 = (a10[tid] + b[tid]*n[tid])*log(Vs30_star/Vlin[tid])
		else
			# need nonlinear response - requires computation of the median Sa at 1180m/s
			ask1180 = PJSask2014(T, M, Rrup, Rjb, Rx, Ry0, Ztor, W, Dip, Fnm, Frv, Fhw, 1180.0, NaN, 1, region, Fas, CRjb)
			@inbounds f5 = a10[tid]*log(Vs30_star/Vlin[tid]) - b[tid]*log(ask1180.IM + c[tid]) + b[tid]*log(ask1180.IM + c[tid]*(Vs30_star/Vlin[tid])^n[tid])
		end

		# need to compute the reference Z1p0ref
		if region == "Japan"
			Z1p0ref = (1.0/1000.0) * exp( -7.67/4.0 * log((Vs30^4 + 610.0^4)/(1360.0^4 + 610.0^4)) )
		else
			Z1p0ref = (1.0/1000.0) * exp( -5.23/2.0 * log((Vs30^2 + 412.0^2)/(1360.0^2 + 412.0^2)) )
		end
		if isnan(Z1p0)
			Z1p0 = Z1p0ref
		end

		# hanging wall model
		if Fhw == 0
			f4 = 0.0
		else
			# require hanging wall terms (in some cases)
			R1 = W*cosd(Dip)
			R2 = 3*R1
			if M <= 5.5 || Rx > R2
				# event is too small, or site is too far away to be effected
				f4 = 0.0
			else
				# taper 1
				if Dip > 30.0
					T1 = (90.0 - Dip)/45.0
				else
					T1 = 60.0/45.0
				end

				# taper 2
				a2HW = 0.2
				if M < 6.5
					T2 = 1.0 + a2HW*(M - 6.5) - (1.0 - a2HW)*(M - 6.5)^2
				else
					T2 = 1.0 + a2HW*(M - 6.5)
				end

				# taper 3
				if Rx < R1
					h1 = 0.25
					h2 = 1.5
					h3 = -0.75
					T3 = h1 + h2*(Rx/R1) + h3*(Rx/R1)^2
				elseif Rx <= R2
					T3 = 1.0 - (Rx-R1)/(R2-R1)
				else
					T3 = 0.0
				end

				# taper 4
				if Ztor <= 10.0
					T4 = 1.0 - Ztor^2/100.0
				else
					T4 = 0.0
				end

				# taper 5
				if isnan(Ry0)
					# use the Rjb taper
					if Rjb == 0.0
						T5 = 1.0
					elseif Rjb < 30.
						T5 = 1.0 - Rjb/30.0
					else
						T5 = 0.0
					end
				else
					Ry1 = Rx*tand(20.0)
					if (Ry0-Ry1) <= 0.0
						T5 = 1.0
					elseif (Ry0-Ry1) < 5.0
						T5 = 1.0 - (Ry0-Ry1)/5.0
					else
						T5 = 0.0
					end
				end

				@inbounds f4 = a13[tid]*T1*T2*T3*T4*T5
			end
		end

		# depth to top of rupture
		if Ztor < 20.0
			@inbounds f6 = a15[tid]*Ztor/20.0
		else
			@inbounds f6 = a15[tid]
		end

		# soil depth model (interpolate coefficients from the binned values of the paper)
		if Vs30 <= 150.0
			@inbounds aZ = a43[tid]
		elseif Vs30 <= 250.0
			@inbounds aZ = a43[tid] + ((Vs30 - 150.0)/100.0)*(a44[tid] - a43[tid])
		elseif Vs30 <= 400.0
			@inbounds aZ = a44[tid] + ((Vs30 - 250.0)/150.0)*(a45[tid] - a44[tid])
		elseif Vs30 <= 700.0
			@inbounds aZ = a45[tid] + ((Vs30 - 400.0)/300.0)*(a46[tid] - a45[tid])
		else
			aZ = a46[tid]
		end
		f10 = aZ * log( (Z1p0 + 0.01)/(Z1p0ref + 0.01) )

		# aftershock scaling
		if Fas == 1
			if CRjb <= 5.0
				@inbounds f11 = a14[tid]
			elseif CRjb < 15.0
				@inbounds f11 = a14[tid]*( 1.0 - (CRjb-5.0)/10.0 )
			else
				f11 = 0.0
			end
		else
			f11 = 0.0
		end

		# regionalization
		if region == "Japan"
			# linear site response adjustment
			if Vs30 < 150.0
				@inbounds aV_JP = a36[tid]
			elseif Vs30 < 250.0
				@inbounds aV_JP = a36[tid] + ((Vs30 - 150.0)/100.0)*(a37[tid] - a36[tid])
			elseif Vs30 < 350.0
				@inbounds aV_JP = a37[tid] + ((Vs30 - 250.0)/100.0)*(a38[tid] - a37[tid])
			elseif Vs30 < 450.0
				@inbounds aV_JP = a38[tid] + ((Vs30 - 350.0)/100.0)*(a39[tid] - a38[tid])
			elseif Vs30 < 600.0
				@inbounds aV_JP = a39[tid] + ((Vs30 - 450.0)/150.0)*(a40[tid] - a39[tid])
			elseif Vs30 < 850.0
				@inbounds aV_JP = a40[tid] + ((Vs30 - 600.0)/250.0)*(a41[tid] - a40[tid])
			elseif Vs30 < 1150.0
				@inbounds aV_JP = a41[tid] + ((Vs30 - 850.0)/300.0)*(a42[tid] - a41[tid])
			else
				@inbounds aV_JP = a42[tid]
			end
			@inbounds Δregion = aV_JP + a29[tid]*Rrup
		elseif region == "Taiwan"
			@inbounds Δregion = a31[tid]*log(Vs30_star / Vlin[tid]) + a25[tid]
		elseif region == "China"
			@inbounds Δregion = a28[tid]*Rrup
		else
			Δregion = 0.0
		end

		lnSa = f1 + f7 + f8 + f11 + f5 + f4 + f6 + f10 + Δregion
		Sa = exp(lnSa)

		# aleatory variability
		# within-event component
		if region == "Japan"
			if Rrup < 30.0
				@inbounds ϕAL = s5[tid]
			elseif Rrup <= 80.0
				@inbounds ϕAL = s5[tid] + (s6[tid]-s5[tid])/50.0*(Rrup - 30.0)
			else
				@inbounds ϕAL = s6[tid]
			end
		else
			if Vs30known == 1
				# known Vs30 case
				@inbounds s1 = s1_VsK[tid]
				@inbounds s2 = s2_VsK[tid]
			else
				# Vs30 estimated
				@inbounds s1 = s1_VsU[tid]
				@inbounds s2 = s2_VsU[tid]
			end

			if M < 4.0
				ϕAL = s1
			elseif M <= 6.0
				ϕAL = s1 + (s2-s1)*(M - 4.0)/2.0
			else
				ϕAL = s2
			end
		end
		# between-event component
		if M < 5.0
			@inbounds τAL = s3[tid]
		elseif M <= 7.0
			@inbounds τAL = s3[tid] + (s4[tid]-s3[tid])*(M - 5.0)/2.0
		else
			@inbounds τAL = s4[tid]
		end

		ϕAmp = 0.4
		# Note that this is a correction to ASK14, as their coefficients don't make sense
		ϕB2 = ϕAL^2 - ϕAmp^2
		# ϕB = ( ϕAL < ϕAmp ) ? 0.0 : sqrt( ϕAL^2 - ϕAmp^2 )
		# ϕB = sqrt( ϕAL^2 - ϕAmp^2 )
		τB = τAL

		if Vs30 >= Vlin[tid]
			∂lnAmp = 0.0
		else
			@fastmath @inbounds ∂lnAmp = b[tid]*ask1180.IM * ( 1.0/(ask1180.IM + c[tid]*(Vs30/Vlin[tid])^n[tid]) - 1.0/(ask1180.IM + c[tid]) )
		end

		ϕ = sqrt( ϕB2*(1.0 + ∂lnAmp)^2 + ϕAmp^2 )
		# ϕ = sqrt( ϕB^2*(1.0 + ∂lnAmp)^2 + ϕAmp^2 )
		τ = τB*(1.0 + ∂lnAmp)
		σ = sqrt( τ^2 + ϕ^2 )
	else
		# we need to interpolate
		@inbounds TiLO = Ti[tidLO]
		@inbounds TiHI = Ti[tidHI]
		Tfrac = log(T/TiLO)/log(TiHI/TiLO)

		askLO = PJSask2014( TiLO, M, Rrup, Rjb, Rx, Ry0, Ztor, W, Dip, Fnm, Frv, Fhw, Vs30, Z1p0, Vs30known, region, Fas, CRjb )
		askHI = PJSask2014( TiHI, M, Rrup, Rjb, Rx, Ry0, Ztor, W, Dip, Fnm, Frv, Fhw, Vs30, Z1p0, Vs30known, region, Fas, CRjb )

		lnSa = askLO.lnIM + Tfrac*(askHI.lnIM-askLO.lnIM)
		Sa = exp( lnSa )

		τ = askLO.τ + Tfrac*(askHI.τ - askLO.τ)
		ϕ = askLO.ϕ + Tfrac*(askHI.ϕ - askLO.ϕ)
		σ = askLO.σ + Tfrac*(askHI.σ - askLO.σ)
	end

	gm = PJSgroundMotion(Sa, lnSa, τ, ϕ, σ)
	return gm
end

# T = 0.01
# M = 6.0
# Rrup = 10.0
# Rjb = 10.0
# Rx = 10.0
# Ry0 = 0.0
# Ztor = 0.0
# W = 8.0
# Dip = 90.0
# Fnm = 0
# Frv = 0
# Fhw = 0
# Vs30 = 760.0
# Z1p0 = NaN
# Vs30known = 1
# region = "California"
# Fas = 0
# CRjb = NaN
#
# @time ask = PJSask2014(T, M, Rrup, Rjb, Rx, Ry0, Ztor, W, Dip, Fnm, Frv, Fhw, Vs30, Z1p0, Vs30known, region, Fas, CRjb)
# @time ask = PJSask2014(0.015, M, Rrup, Rjb, Rx, Ry0, Ztor, W, Dip, Fnm, Frv, Fhw, Vs30, Z1p0, Vs30known, region, Fas, CRjb)

function PJSask2014( Ti::Vector{U}, M::U, Rrup::U, Rjb::U, Rx::U, Ry0::U, Ztor::U, W::U, Dip::U, Fnm::Int64, Frv::Int64, Fhw::Int64, Vs30::U, Z1p0::U, Vs30known::Int64=1, region::String="California", Fas::Int64=0, CRjb::U=NaN ) where U <: Real
	predictions = Array{PJSgroundMotion,1}()
	for T in Ti
		push!(predictions, PJSask2014( T, M, Rrup, Rjb, Rx, Ry0, Ztor, W, Dip, Fnm, Frv, Fhw, Vs30, Z1p0, Vs30known, region, Fas, CRjb ))
	end
	return predictions
end



function PJSask2014(T::U, rup::Rupture{U}, site::Site{U}) where U<:Real
	Rrup = rupture_distance(rup, site)
	Rjb = joyner_boore_distance(rup, site)
	Rx = strike_distance(rup, site)
	Ry0 = strike_parallel_distance(rup, site)
	Ztor = ztor(rup)
	W = rup.rupWidth
	Dip = dip(rup)
	Vs30 = vs30(site)
	Z1p0 = z1p0(site) / 1000.0 # depth in km
	Vs30known = site.measured
	region = site.region
	Fss, Fnm, Frv, Fuk = mechanism(rup)
	Fhw = (Rx < 0) ? 0 : 1

	gm = PJSask2014(T, rup.magnitude, Rrup, Rjb, Rx, Ry0, Ztor, W, Dip, Fnm, Frv, Fhw, Vs30, Z1p0, Vs30known, region )
	return gm
end


function PJSask2014(Ti::Vector{U}, rup::Rupture{U}, site::Site{U}) where U<:Real
	predictions = Array{PJSgroundMotion,1}()
	for T in Ti
		push!(predictions, PJSask2014(T, rup, site))
	end
	return predictions
end



# function get_spectrum( v::Vector{PJSgmpe} )
# 	Sai = zeros(length(v),1)
# 	for i = 1:length(v)
# 		@inbounds Sai[i] = v[i].Sa
# 	end
# 	return Sai
# end

# Ti = exp10.(range(-2.0, stop=1.0, step=0.05))
#
# @time aski = PJSask2014(Ti, M, Rrup, Rjb, Rx, Ry0, Ztor, W, Dip, Fnm, Frv, Fhw, Vs30, Z1p0, Vs30known, region, Fas, CRjb)
# Sai = get_spectrum(aski)
#
# using Plots
# pyplot()
# theme(:juno)
#
# plot(Ti, Sai, lab="", xaxis=:log10, yaxis=:log10)
#
# M = 7.0
# @time ask0p0 = PJSask2014(0.0, M, Rrup, Rjb, Rx, Ry0, Ztor, W, Dip, Fnm, Frv, Fhw, Vs30, Z1p0, Vs30known, region, Fas, CRjb)
# @time ask0p2 = PJSask2014(0.2, M, Rrup, Rjb, Rx, Ry0, Ztor, W, Dip, Fnm, Frv, Fhw, Vs30, Z1p0, Vs30known, region, Fas, CRjb)
#
# ask0p2.Sa / ask0p0.Sa
