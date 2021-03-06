

function PJScy2014( T::U, M::U, Rrup::U, Rjb::U, Rx::U, ΔZtor::U, Dip::U, Fnm::Int64, Frv::Int64, Fhw::Int64, Vs30::U, Vs30meas::Int64=1, ΔZ1p0::U=0.0, ΔDPP::U=0.0 ) where U <: Real
	# Period-dependent coefficients
	Ti = [ -1.0, 0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.1, 0.12, 0.15, 0.17, 0.2, 0.25, 0.3, 0.4, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 7.5, 10.0 ]
	c1 = [ 2.3549, -1.5065, -1.5065, -1.4798, -1.2972, -1.1007, -0.9292, -0.658, -0.5613, -0.5342, -0.5462, -0.5858, -0.6798, -0.8663, -1.0514, -1.3794, -1.6508, -2.1511, -2.5365, -3.0686, -3.4148, -3.9013, -4.2466, -4.5143, -5.0009, -5.3461 ]
	c1a = [ 0.165, 0.165, 0.165, 0.165, 0.165, 0.165, 0.165, 0.165, 0.165, 0.165, 0.165, 0.165, 0.165, 0.165, 0.165, 0.165, 0.165, 0.165, 0.165, 0.165, 0.1645, 0.1168, 0.0732, 0.0484, 0.022, 0.0124 ]
	c1b = [ -0.0626, -0.255, -0.255, -0.255, -0.255, -0.255, -0.255, -0.254, -0.253, -0.252, -0.25, -0.248, -0.2449, -0.2382, -0.2313, -0.2146, -0.1972, -0.162, -0.14, -0.1184, -0.11, -0.104, -0.102, -0.101, -0.101, -0.1 ]
	c1c = [ -0.165, -0.165, -0.165, -0.165, -0.165, -0.165, -0.165, -0.165, -0.165, -0.165, -0.165, -0.165, -0.165, -0.165, -0.165, -0.165, -0.165, -0.165, -0.165, -0.165, -0.1645, -0.1168, -0.0732, -0.0484, -0.022, -0.0124 ]
	c1d = [ 0.0626, 0.255, 0.255, 0.255, 0.255, 0.255, 0.255, 0.254, 0.253, 0.252, 0.25, 0.248, 0.2449, 0.2382, 0.2313, 0.2146, 0.1972, 0.162, 0.14, 0.1184, 0.11, 0.104, 0.102, 0.101, 0.101, 0.1 ]
	cn = [ 3.3024, 16.0875, 16.0875, 15.7118, 15.8819, 16.4556, 17.6453, 20.1772, 19.9992, 18.7106, 16.6246, 15.3709, 13.7012, 11.2667, 9.1908, 6.5459, 5.2305, 3.7896, 3.3024, 2.8498, 2.5417, 2.1488, 1.8957, 1.7228, 1.5737, 1.5265 ]
	cm = [ 5.423, 4.9993, 4.9993, 4.9993, 4.9993, 4.9993, 4.9993, 5.0031, 5.0172, 5.0315, 5.0547, 5.0704, 5.0939, 5.1315, 5.167, 5.2317, 5.2893, 5.4109, 5.5106, 5.6705, 5.7981, 5.9983, 6.1552, 6.2856, 6.5428, 6.7415 ]
	c3 = [ 2.3152, 1.9636, 1.9636, 1.9636, 1.9636, 1.9636, 1.9636, 1.9636, 1.9636, 1.9795, 2.0362, 2.0823, 2.1521, 2.2574, 2.344, 2.4709, 2.5567, 2.6812, 2.7474, 2.8161, 2.8514, 2.8875, 2.9058, 2.9169, 2.932, 2.9396 ]
	c5 = [ 5.8096, 6.4551, 6.4551, 6.4551, 6.4551, 6.4551, 6.4551, 6.4551, 6.8305, 7.1333, 7.3621, 7.4365, 7.4972, 7.5416, 7.56, 7.5735, 7.5778, 7.5808, 7.5814, 7.5817, 7.5818, 7.5818, 7.5818, 7.5818, 7.5818, 7.5818 ]
	chm = [ 3.0514, 3.0956, 3.0956, 3.0963, 3.0974, 3.0988, 3.1011, 3.1094, 3.2381, 3.3407, 3.43, 3.4688, 3.5146, 3.5746, 3.6232, 3.6945, 3.7401, 3.7941, 3.8144, 3.8284, 3.833, 3.8361, 3.8369, 3.8376, 3.838, 3.838 ]
	c6 = [ 0.4407, 0.4908, 0.4908, 0.4925, 0.4992, 0.5037, 0.5048, 0.5048, 0.5048, 0.5048, 0.5045, 0.5036, 0.5016, 0.4971, 0.4919, 0.4807, 0.4707, 0.4575, 0.4522, 0.4501, 0.45, 0.45, 0.45, 0.45, 0.45, 0.45 ]
	c7 = [ 0.0324, 0.0352, 0.0352, 0.0352, 0.0352, 0.0352, 0.0352, 0.0352, 0.0352, 0.0352, 0.0352, 0.0352, 0.0352, 0.0352, 0.0352, 0.0352, 0.0352, 0.0352, 0.0352, 0.0352, 0.0352, 0.016, 0.0062, 0.0029, 0.0007, 0.0003 ]
	c7b = [ 0.0097, 0.0462, 0.0462, 0.0472, 0.0533, 0.0596, 0.0639, 0.063, 0.0532, 0.0452, 0.0345, 0.0283, 0.0202, 0.009, -0.0004, -0.0155, -0.0278, -0.0477, -0.0559, -0.063, -0.0665, -0.0516, -0.0448, -0.0424, -0.0348, -0.0253 ]
	c8 = [ 0.2154, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0991, 0.1982, 0.2154, 0.2154, 0.2154, 0.2154, 0.2154, 0.2154, 0.2154, 0.2154 ]
	c8b = [ 5.0, 0.4833, 0.4833, 1.2144, 1.6421, 1.9456, 2.181, 2.6087, 2.9122, 3.1045, 3.3399, 3.4719, 3.6434, 3.8787, 4.0711, 4.3745, 4.6099, 5.0376, 5.3411, 5.7688, 6.0723, 6.5, 6.8035, 7.0389, 7.4666, 7.77 ]
	c9 = [ 0.3079, 0.9228, 0.9228, 0.9296, 0.9396, 0.9661, 0.9794, 1.026, 1.0177, 1.0008, 0.9801, 0.9652, 0.9459, 0.9196, 0.8829, 0.8302, 0.7884, 0.6754, 0.6196, 0.5101, 0.3917, 0.1244, 0.0086, 0.0, 0.0, 0.0 ]
	c9a = [ 0.1, 0.1202, 0.1202, 0.1217, 0.1194, 0.1166, 0.1176, 0.1171, 0.1146, 0.1128, 0.1106, 0.115, 0.1208, 0.1208, 0.1175, 0.106, 0.1061, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1, 0.1 ]
	c9b = [ 6.5, 6.8607, 6.8607, 6.8697, 6.9113, 7.0271, 7.0959, 7.3298, 7.2588, 7.2372, 7.2109, 7.2491, 7.2988, 7.3691, 6.8789, 6.5334, 6.526, 6.5, 6.5, 6.5, 6.5, 6.5, 6.5, 6.5, 6.5, 6.5 ]
	c11 = [ 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ]
	c11b = [ -0.3834, -0.4536, -0.4536, -0.4536, -0.4536, -0.4536, -0.4536, -0.4536, -0.4536, -0.4536, -0.4536, -0.4536, -0.444, -0.3539, -0.2688, -0.1793, -0.1428, -0.1138, -0.1062, -0.102, -0.1009, -0.1003, -0.1001, -0.1001, -0.1, -0.1 ]
	γ1 = [ -0.001852, -0.007146, -0.007146, -0.007249, -0.007869, -0.008316, -0.008743, -0.009537, -0.00983, -0.009913, -0.009896, -0.009787, -0.009505, -0.008918, -0.008251, -0.007267, -0.006492, -0.005147, -0.004277, -0.002979, -0.002301, -0.001344, -0.001084, -0.00101, -0.000964, -0.00095 ]
	γ2 = [ -0.007403, -0.006758, -0.006758, -0.006758, -0.006758, -0.006758, -0.006758, -0.00619, -0.005332, -0.004732, -0.003806, -0.00328, -0.00269, -0.002128, -0.001812, -0.001274, -0.001074, -0.001115, -0.001197, -0.001675, -0.002349, -0.003306, -0.003566, -0.00364, -0.003686, -0.0037 ]
	γ3 = [ 4.3439, 4.2542, 4.2542, 4.2386, 4.2519, 4.296, 4.3578, 4.5455, 4.7603, 4.8963, 5.0644, 5.1371, 5.188, 5.2164, 5.1954, 5.0899, 4.7854, 4.3304, 4.1667, 4.0029, 3.8949, 3.7928, 3.7443, 3.709, 3.6632, 3.623 ]
	ϕ1 = [ -0.7936, -0.521, -0.521, -0.5055, -0.4368, -0.3752, -0.3469, -0.3747, -0.444, -0.4895, -0.5477, -0.5922, -0.6693, -0.7766, -0.8501, -0.9431, -1.0044, -1.0602, -1.0941, -1.1142, -1.1154, -1.1081, -1.0603, -0.9872, -0.8274, -0.7053 ]
	ϕ2 = [ -0.0699, -0.1417, -0.1417, -0.1364, -0.1403, -0.1591, -0.1862, -0.2538, -0.2943, -0.3077, -0.3113, -0.3062, -0.2927, -0.2662, -0.2405, -0.1975, -0.1633, -0.1028, -0.0699, -0.0425, -0.0302, -0.0129, -0.0016, 0.0, 0.0, 0.0 ]
	ϕ3 = [ -0.008444, -0.00701, -0.00701, -0.007279, -0.007354, -0.006977, -0.006467, -0.005734, -0.005604, -0.005696, -0.005845, -0.005959, -0.006141, -0.006439, -0.006704, -0.007125, -0.007435, -0.00812, -0.008444, -0.007707, -0.004792, -0.001828, -0.001523, -0.00144, -0.001369, -0.001361 ]
	ϕ4 = [ 5.41, 0.102151, 0.102151, 0.10836, 0.119888, 0.133641, 0.148927, 0.190596, 0.230662, 0.253169, 0.266468, 0.26506, 0.255253, 0.231541, 0.207277, 0.165464, 0.133828, 0.085153, 0.058595, 0.031787, 0.019716, 0.009643, 0.005379, 0.003223, 0.001134, 0.000515 ]
	ϕ5 = [ 0.0202, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.001, 0.004, 0.01, 0.034, 0.067, 0.143, 0.203, 0.277, 0.309, 0.321, 0.329, 0.33 ]
	ϕ6 = [ 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0 ]
	τ1 = [ 0.3894, 0.4, 0.4, 0.4026, 0.4063, 0.4095, 0.4124, 0.4179, 0.4219, 0.4244, 0.4275, 0.4292, 0.4313, 0.4341, 0.4363, 0.4396, 0.4419, 0.4459, 0.4484, 0.4515, 0.4534, 0.4558, 0.4574, 0.4584, 0.4601, 0.4612 ]
	τ2 = [ 0.2578, 0.26, 0.26, 0.2637, 0.2689, 0.2736, 0.2777, 0.2855, 0.2913, 0.2949, 0.2993, 0.3017, 0.3047, 0.3087, 0.3119, 0.3165, 0.3199, 0.3255, 0.3291, 0.3335, 0.3363, 0.3398, 0.3419, 0.3435, 0.3459, 0.3474 ]
	σ1 = [ 0.4785, 0.4912, 0.4912, 0.4904, 0.4988, 0.5049, 0.5096, 0.5179, 0.5236, 0.527, 0.5308, 0.5328, 0.5351, 0.5377, 0.5395, 0.5422, 0.5433, 0.5294, 0.5105, 0.4783, 0.4681, 0.4617, 0.4571, 0.4535, 0.4471, 0.4426 ]
	σ2 = [ 0.36292, 0.3762, 0.3762, 0.3762, 0.3849, 0.391, 0.3957, 0.4043, 0.4104, 0.4143, 0.4191, 0.4217, 0.4252, 0.4299, 0.4338, 0.4399, 0.4446, 0.4533, 0.4594, 0.468, 0.4681, 0.4617, 0.4571, 0.4535, 0.4471, 0.4426 ]
	σ3 = [ 0.7504, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.8, 0.7999, 0.7997, 0.7988, 0.7966, 0.7792, 0.7504, 0.7136, 0.7035, 0.7006, 0.7001, 0.7, 0.7, 0.7 ]
	# Don't work with the Japan-specific or Wenchuan-specific predictions
	# σ2_Jp = [ 0.3918, 0.4528, 0.4528, 0.4551, 0.4571, 0.4642, 0.4716, 0.5022, 0.523, 0.5278, 0.5304, 0.531, 0.5312, 0.5309, 0.5307, 0.531, 0.5313, 0.5309, 0.5302, 0.5276, 0.5167, 0.4917, 0.4682, 0.4517, 0.4167, 0.3755 ]
	# γ_Jp_It = [ 2.2306, 1.5817, 1.5817, 1.574, 1.5544, 1.5502, 1.5391, 1.4804, 1.4094, 1.3682, 1.3241, 1.3071, 1.2931, 1.315, 1.3514, 1.4051, 1.4402, 1.528, 1.6523, 1.8872, 2.1348, 3.5752, 3.8646, 3.7292, 2.3763, 1.7679 ]
	# γ_Wn = [ 0.335, 0.7594, 0.7594, 0.7606, 0.7642, 0.7676, 0.7739, 0.7956, 0.7932, 0.7768, 0.7437, 0.7219, 0.6922, 0.6579, 0.6362, 0.6049, 0.5507, 0.3582, 0.2003, 0.0356, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0 ]
	# ϕ1_Jp = [ -0.7966, -0.6846, -0.6846, -0.6681, -0.6314, -0.5855, -0.5457, -0.4685, -0.4985, -0.5603, -0.6451, -0.6981, -0.7653, -0.8469, -0.8999, -0.9618, -0.9945, -1.0225, -1.0002, -0.9245, -0.8626, -0.7882, -0.7195, -0.656, -0.5202, -0.4068 ]
	# ϕ5_Jp = [ 0.9488, 0.459, 0.459, 0.458, 0.462, 0.453, 0.436, 0.383, 0.375, 0.377, 0.379, 0.38, 0.384, 0.393, 0.408, 0.462, 0.524, 0.658, 0.78, 0.96, 1.11, 1.291, 1.387, 1.433, 1.46, 1.464 ]
	# ϕ6_Jp = [ 800.0, 800.0, 800.0, 800.0, 800.0, 800.0, 800.0, 800.0, 800.0, 800.0, 800.0, 800.0, 800.0, 800.0, 800.0, 800.0, 800.0, 800.0, 800.0, 800.0, 800.0, 800.0, 800.0, 800.0, 800.0, 800.0 ]

	# Period-independent coefficients
	c2 = 1.06
	c4 = -2.1
	c4a = -0.5
	crb = 50.0
	c8a = 0.2695

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

		# need to compute the expected depth to top of rupture
		if Frv == 1
			EZtor = max(2.704 - 1.226*max(M-5.849,0.0),0.0)^2
		else
			EZtor = max(2.673 - 1.136*max(M-4.970,0.0),0.0)^2
		end
		Ztor = EZtor + ΔZtor

		@fastmath @inbounds lnYref = c1[tid] +
			( c1a[tid] + c1c[tid]/cosh(2.0*max(M-4.5,0.0)) )*Frv +
			( c1b[tid] + c1d[tid]/cosh(2.0*max(M-4.5,0.0)) )*Fnm +
			( c7[tid]  + c7b[tid]/cosh(2.0*max(M-4.5,0.0)) )*ΔZtor +
			( c11[tid] + c11b[tid]/cosh(2.0*max(M-4.5,0.0)) )*cosd(Dip)^2 +
			c2*(M-6.0) + ((c2-c3[tid])/cn[tid])*log( 1.0 + exp( cn[tid]*(cm[tid]-M) ) ) +
			c4*log( Rrup + c5[tid]*cosh( c6[tid]*max(M-chm[tid],0.0) ) ) +
			(c4a-c4)*log( sqrt( Rrup^2 + crb^2 ) ) +
			( γ1[tid] + γ2[tid]/cosh(max(M-γ3[tid],0.0)) )*Rrup +
			c8[tid]*max( 1.0 - max(Rrup-40.0,0.0)/30.0, 0.0 )*min( max(M-5.5,0.0)/0.8, 1.0 )*exp( -c8a*(M-c8b[tid])^2 )*ΔDPP +
			c9[tid]*Fhw*cosd(Dip)*( c9a[tid] + (1.0-c9a[tid])*tanh( Rx / c9b[tid] ) )*( 1.0 - sqrt(Rjb^2 + Ztor^2)/(Rrup+1.0) )
		@fastmath @inbounds lnSa = lnYref +
			ϕ1[tid]*min(log(Vs30/1130.0),0.0) +
			ϕ2[tid]*( exp(ϕ3[tid]*(min(Vs30,1130.0)-360.0)) - exp(ϕ3[tid]*(1130.0-360.0)))*log( ( exp(lnYref) + ϕ4[tid] ) / ϕ4[tid] ) +
			ϕ5[tid]*( 1.0 - exp( -ΔZ1p0/ϕ6[tid] ) )
		Sa = exp( lnSa )

		@fastmath @inbounds NL = ϕ2[tid]*( exp(ϕ3[tid]*(min(Vs30,1130.0)-360.0)) - exp(ϕ3[tid]*(1130.0-360.0))) * ( exp(lnYref) / ( exp(lnYref) + ϕ4[tid] ) )
		@inbounds τ = ( τ1[tid] + ( τ2[tid] - τ1[tid] )*( min( max( M, 5.0 ), 6.5 ) - 5.0 ) / 1.5 )
		@inbounds ϕ = ( σ1[tid] + ( σ2[tid] - σ1[tid] )*( min( max( M, 5.0 ), 6.5 ) - 5.0 ) / 1.5 ) * sqrt( ifelse( Vs30meas==1, 0.7, σ3[tid] ) + ( 1.0 + NL )^2 )
		σ = sqrt( ((1.0 + NL)*τ)^2 + ϕ^2 )

	else
		# we need to interpolate
		@inbounds TiLO = Ti[tidLO]
		@inbounds TiHI = Ti[tidHI]
		Tfrac = log(T/TiLO)/log(TiHI/TiLO)

		cyLO = PJScy2014( TiLO, M, Rrup, Rjb, Rx, ΔZtor, Dip, Fnm, Frv, Fhw, Vs30, Vs30meas, ΔZ1p0, ΔDPP )
		cyHI = PJScy2014( TiHI, M, Rrup, Rjb, Rx, ΔZtor, Dip, Fnm, Frv, Fhw, Vs30, Vs30meas, ΔZ1p0, ΔDPP )

		lnSa = cyLO.lnIM + Tfrac*(cyHI.lnIM-cyLO.lnIM)
		Sa = exp( lnSa )

		τ = cyLO.τ + Tfrac*(cyHI.τ - cyLO.τ)
		ϕ = cyLO.ϕ + Tfrac*(cyHI.ϕ - cyLO.ϕ)
		σ = cyLO.σ + Tfrac*(cyHI.σ - cyLO.σ)
	end

	gm = PJSgroundMotion(Sa, lnSa, τ, ϕ, σ)
	return gm
end




function PJScy2014( Ti::Vector{U}, M::U, Rrup::U, Rjb::U, Rx::U, ΔZtor::U, Dip::U, Fnm::Int64, Frv::Int64, Fhw::Int64, Vs30::U, Vs30meas::Int64=1, ΔZ1p0::U=0.0, ΔDPP::U=0.0 ) where U <: Real
	predictions = Array{PJSgroundMotion,1}()
	for T in Ti
		push!(predictions, PJScy2014( T, M, Rrup, Rjb, Rx, ΔZtor, Dip, Fnm, Frv, Fhw, Vs30, Vs30meas, ΔZ1p0, ΔDPP ))
	end
	return predictions
end


function PJScy2014(T::U, rup::Rupture{U}, site::Site{U}) where U<:Real
	Rrup = rupture_distance(rup, site)
	Rjb = joyner_boore_distance(rup, site)
	Rx = strike_distance(rup, site)
	Ztor = ztor(rup)
	Zhyp = zhyp(rup)
	W = rup.rupWidth
	Dip = dip(rup)
	Vs30 = vs30(site)
	Z1p0 = z1p0(site)
	Fss, Fnm, Frv, Fuk = mechanism(rup)

	if Frv == 1
		EZtor = max(2.704 - 1.226*max(rup.magnitude-5.849,0.0),0.0)^2
	else
		EZtor = max(2.673 - 1.136*max(rup.magnitude-4.970,0.0),0.0)^2
	end
	ΔZtor = Ztor - EZtor

	# get expected Z1p0 (California and non-Japan)
	EZ1p0 = exp( -7.15/4 * log( (Vs30^4 + 570.94^4)/(1360.0^4 + 570.94^4) ) )
	# EZ1p0 = exp( -5.23/2 * log( (Vs30^2 + 412.0^2)/(1360.0^2 + 412.0^2) ) )
	ΔZ1p0 = Z1p0 - EZ1p0

	Fhw = (Rx < 0.0) ? 0 : 1

	gm = PJScy2014(T, rup.magnitude, Rrup, Rjb, Rx, ΔZtor, Dip, Fnm, Frv, Fhw, Vs30, site.measured, ΔZ1p0 )
	return gm
end


function PJScy2014(Ti::Vector{U}, rup::Rupture{U}, site::Site{U}) where U<:Real
	predictions = Array{PJSgroundMotion,1}()
	for T in Ti
		push!(predictions, PJScy2014(T, rup, site))
	end
	return predictions
end
