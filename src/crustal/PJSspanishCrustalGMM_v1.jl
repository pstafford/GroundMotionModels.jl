# File implementing the Spanish Crustal GMM
# The model uses CY14 as its base and applies a stress drop correction and a distance
# correction - both for geometric spreading and anelastic attenuation



"""
	PJScy2014adjusted( T::U, M::U, Rrup::U, Rjb::U, Rx::U, ΔZtor::U, Dip::U, Fnm::Int64, Frv::Int64, Fhw::Int64, Vs30::U, Vs30meas::Int64, Δc1::U=0.0, Δcm::U=0.0, Δγ1::U=0.0, ΔZ1p0::U=0.0, ΔDPP::U=0.0 ) where U<:Real

Modification to the Chiou & Youngs (2014) gmm to account for theoretical stress drop adjustments as well as overall source and anelastic attenuation adjustment.
Note that all parameters are as in the original model with the
exception of Δc1, Δcm & Δγ1 which are:
- the change in overall source amplitude, change to c1
- the change to the cm coefficient required to account for
a particular level of stress drop.
- the change in overall anelastic attenuation rate γ1
"""
function PJScy2014adjusted( T::U, M::U, Rrup::U, Rjb::U, Rx::U, ΔZtor::U, Dip::U, Fnm::Int64, Frv::Int64, Fhw::Int64, Vs30::U, Vs30meas::Int64, Δc1::U=0.0, Δcm::U=0.0, Δγ1::U=0.0, ΔZ1p0::U=0.0, ΔDPP::U=0.0 ) where U<:Real
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

	# Period-independent coefficients
	c2 = 1.06
	c4 = -2.1
	c4a = -0.5
	crb = 50.0
	c8a = 0.2695
	c3t = 3.45

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

		# find the additional adjustment required for theoretical consistency
		Δy = (c3[tid]-c3t)*Δcm
		α = atan(c3[tid]) - atan(c2)
		γ = atan(1.0/c3[tid])
		β = π - α - γ
		z = Δy * sin(β) / sin(α)
		δm = z * cos(atan(c3[tid]))

		@inbounds lnSa = c1[tid] + Δc1 +
			( c1a[tid] + c1c[tid]/cosh(2.0*max(M-4.5,0.0)) )*Frv +
			( c1b[tid] + c1d[tid]/cosh(2.0*max(M-4.5,0.0)) )*Fnm +
			( c7[tid]  + c7b[tid]/cosh(2.0*max(M-4.5,0.0)) )*ΔZtor +
			( c11[tid] + c11b[tid]/cosh(2.0*max(M-4.5,0.0)) )*(cosd(Dip)^2) +
			c2*(M-6.0) + ((c2-c3[tid])/cn[tid])*log(1.0 + exp( cn[tid]*(cm[tid]+(Δcm-δm)-M) )) - (c2-c3[tid])*(Δcm-δm) +
			c4*log( Rrup + c5[tid]*cosh( c6[tid]*max(M-chm[tid],0.0) ) ) +
			(c4a-c4)*log( sqrt( Rrup^2 + crb^2 ) ) +
			min( 0.0, Δγ1 + γ1[tid] + γ2[tid]/cosh(max(M-γ3[tid],0.0)) )*Rrup +
			c8[tid]*max( 1.0 - max(Rrup-40.0,0.0)/30.0, 0.0 )*min( max(M-5.5,0.0)/0.8, 1.0 )*exp( -c8a*(M-c8b[tid])^2 )* ΔDPP +
			c9[tid]*Fhw*cosd(Dip)*( c9a[tid] + (1.0-c9a[tid])*tanh( Rx / c9b[tid] ) )*( 1.0 - sqrt(Rjb^2+Ztor^2)/(Rrup+1.0) ) +
			ϕ1[tid]*min(log(Vs30/1130.0),0.0)
		Sa = exp( lnSa )

		# no nonlinearity derivative as there is no nonlinearity
		NL = 0.0
		@inbounds τ = ( τ1[tid] + ( τ2[tid] - τ1[tid] )*( min( max( M, 5.0 ), 6.5 ) - 5.0 ) / 1.5 ) * (1.0 + NL)
		@inbounds ϕ = ( σ1[tid] + ( σ2[tid] - σ1[tid] )*( min( max( M, 5.0 ), 6.5 ) - 5.0 ) / 1.5 ) * sqrt( ifelse( Vs30meas==1, 0.7, σ3[tid] ) + ( 1.0 + NL )^2 )
		σ = sqrt( ((1.0 + NL)*τ)^2 + ϕ^2 )
	else
		# we need to interpolate
		@inbounds TiLO = Ti[tidLO]
		@inbounds TiHI = Ti[tidHI]
		Tfrac = log(T/TiLO)/log(TiHI/TiLO)

		cyLO = PJScy2014adjusted( TiLO, M, Rrup, Rjb, Rx, ΔZtor, Dip, Fnm, Frv, Fhw, Vs30, Vs30meas, Δc1, Δcm, Δγ1, ΔZ1p0, ΔDPP )
		cyHI = PJScy2014adjusted( TiHI, M, Rrup, Rjb, Rx, ΔZtor, Dip, Fnm, Frv, Fhw, Vs30, Vs30meas, Δc1, Δcm, Δγ1, ΔZ1p0, ΔDPP )

		lnSa = cyLO.lnSa + Tfrac*(cyHI.lnSa-cyLO.lnSa)
		Sa = exp( lnSa )

		τ = cyLO.τ + Tfrac*(cyHI.τ - cyLO.τ)
		ϕ = cyLO.ϕ + Tfrac*(cyHI.ϕ - cyLO.ϕ)
		σ = cyLO.σ + Tfrac*(cyHI.σ - cyLO.σ)
	end

	gm = PJSgroundMotion( Sa, lnSa, τ, ϕ, σ )
	return gm
end




 """
 	PJSdistanceCorrectionCoefficients(T::U, Qbranch="central") where U<:Real

 Function to provide coefficients of the distance correction model for a given period and Q/distance branch
 Argument Qbranch can be one of:
   - "lower":      gives coefficients for lower Q branch and period T
   - "central":    gives coefficients for central Q branch and period T
   - "upper":      gives coefficients for upper Q branch and period T
 """
function PJSdistanceCorrectionCoefficients(T::U, Gbranch="free") where U<:Real
 	if Gbranch == "free"
		b1 = [ 0.145, 0.1336, 0.1137, 0.1086, 0.1106, 0.1131, 0.1368, 0.2061, 0.2214, 0.1961, 0.1733, 0.1593, 0.1441, 0.1396, 0.1379, 0.1383, 0.1397, 0.142, 0.1443, 0.1462 ]
		b2 = [ -0.01839, -0.0155, -0.01123, -0.01189, -0.01423, -0.01288, -0.007528, -0.002698, -0.00553, -0.01742, -0.02715, -0.03353, -0.04154, -0.04486, -0.04751, -0.04861, -0.04994, -0.05213, -0.0543, -0.05595 ]
		b3 = [ 0.001094, 0.001048, 0.0009452, 0.0008488, 0.0007225, 0.0003673, 0.0004074, 0.001041, 0.001696, 0.002521, 0.002813, 0.002864, 0.002685, 0.002438, 0.00203, 0.001734, 0.001336, 0.0008808, 0.0005649, 0.0003679 ]
		b4 = [ -0.0005396, -0.000548, -0.0005611, -0.0005664, -0.0005741, -0.0006163, -0.0006436, -0.0006194, -0.0005363, -0.0003529, -0.0002225, -0.0001395, -0.00003534, 0.000009043, 0.00004487, 0.00005784, 0.00006509, 0.00005772, 0.00003939, 0.00002052 ]
		b5 = [ 0.7139, 0.7105, 0.6859, 0.656, 0.6588, 0.7026, 0.6656, 0.5268, 0.5104, 0.5906, 0.6593, 0.7046, 0.7621, 0.7868, 0.8075, 0.8163, 0.8247, 0.8316, 0.835, 0.8369 ]
		b6 = [ 1.011, 1.01, 1.011, 1.037, 1.083, 1.164, 1.158, 0.9825, 0.8828, 0.8851, 0.9115, 0.9285, 0.9469, 0.9528, 0.9561, 0.9569, 0.9573, 0.9577, 0.9578, 0.9579 ]
		b7 = [ -0.2875, -0.2557, -0.1964, -0.1373, -0.0275, 0.3555, 0.468, 0.1966, -0.1525, -0.65, -0.9069, -1.042, -1.18, -1.226, -1.252, -1.257, -1.252, -1.233, -1.211, -1.195 ]
	elseif Gbranch == "1p00"
		b1 = [ 0.02374, 0.01264, -0.006637, -0.01314, -0.01455, -0.01929, 0.0031, 0.07967, 0.1034, 0.08122, 0.05789, 0.04337, 0.02767, 0.02294, 0.02116, 0.02158, 0.02299, 0.02531, 0.02758, 0.02941 ]
		b2 = [ 0.01044, 0.01333, 0.01761, 0.01695, 0.01461, 0.01596, 0.02131, 0.02614, 0.02331, 0.01142, 0.001688, -0.004697, -0.0127, -0.01603, -0.01867, -0.01978, -0.02111, -0.0233, -0.02546, -0.02711 ]
		b3 = [ 0.0009968, 0.0009499, 0.0008448, 0.0007488, 0.0006257, 0.0002773, 0.0003206, 0.0009643, 0.00161, 0.002422, 0.002712, 0.002762, 0.002582, 0.002335, 0.001926, 0.001631, 0.001232, 0.0007771, 0.0004613, 0.0002645 ]
		b4 = [ -0.0006126, -0.000621, -0.0006341, -0.0006394, -0.0006471, -0.0006893, -0.0007166, -0.0006924, -0.0006093, -0.0004259, -0.0002955, -0.0002125, -0.0001083, -0.00006396, -0.00002813, -0.00001516, -0.000007917, -0.00001528, -0.00003361, -0.00005249 ]
		b5 = [ 0.8491, 0.845, 0.8192, 0.7924, 0.8026, 0.8617, 0.8276, 0.6718, 0.6369, 0.7113, 0.7815, 0.828, 0.8869, 0.912, 0.9329, 0.9418, 0.9501, 0.9571, 0.9606, 0.9625 ]
		b6 = [ 1.041, 1.039, 1.04, 1.065, 1.11, 1.188, 1.186, 1.038, 0.9389, 0.9277, 0.9474, 0.961, 0.9758, 0.9805, 0.9828, 0.9832, 0.9833, 0.9835, 0.9837, 0.9838 ]
		b7 = [ -0.1151, -0.08327, -0.02353, 0.03674, 0.1485, 0.5357, 0.6481, 0.3651, 0.01187, -0.4831, -0.7384, -0.8725, -1.01, -1.055, -1.082, -1.086, -1.081, -1.062, -1.04, -1.024 ]
	elseif Gbranch == "1p15"
		b1 = [ -0.07081, -0.08209, -0.1017, -0.1074, -0.1069, -0.1066, -0.08286, -0.01204, 0.006952, -0.01659, -0.03969, -0.05397, -0.0694, -0.07406, -0.07581, -0.07539, -0.07399, -0.07166, -0.06937, -0.06753 ]
		b2 = [ 0.03336, 0.03624, 0.04052, 0.03986, 0.03752, 0.03887, 0.04422, 0.04905, 0.04622, 0.03433, 0.0246, 0.01822, 0.01021, 0.006887, 0.004238, 0.003136, 0.001805, -0.0003872, -0.002549, -0.004201 ]
		b3 = [ 0.0009426, 0.0008962, 0.0007922, 0.0006961, 0.0005717, 0.0002198, 0.0002601, 0.0008998, 0.001552, 0.00237, 0.002661, 0.002711, 0.002531, 0.002284, 0.001876, 0.00158, 0.001181, 0.0007264, 0.0004106, 0.0002137 ]
		b4 = [ -0.0006688, -0.0006772, -0.0006904, -0.0006957, -0.0007034, -0.0007455, -0.0007728, -0.0007487, -0.0006655, -0.0004822, -0.0003517, -0.0002687, -0.0001646, -0.0001202, -0.00008436, -0.00007139, -0.00006415, -0.00007151, -0.00008984, -0.0001087 ]
		b5 = [ 0.8162, 0.8125, 0.7873, 0.7589, 0.7647, 0.8132, 0.7763, 0.6334, 0.6086, 0.6857, 0.7552, 0.8012, 0.8595, 0.8844, 0.9052, 0.9141, 0.9224, 0.9294, 0.9328, 0.9347 ]
		b6 = [ 1.025, 1.024, 1.025, 1.049, 1.092, 1.166, 1.16, 1.011, 0.9182, 0.9122, 0.9334, 0.9477, 0.9634, 0.9684, 0.971, 0.9716, 0.9718, 0.9721, 0.9722, 0.9723 ]
		b7 = [ -0.1018, -0.06998, -0.01043, 0.04913, 0.1595, 0.5432, 0.6557, 0.3801, 0.02839, -0.4675, -0.7234, -0.8578, -0.9956, -1.041, -1.067, -1.072, -1.067, -1.048, -1.026, -1.01 ]
	else
    	print("PJS::Error, argument 'Gbranch' must be one of 'free', '1p00', or '1p15'")
    	return NaN*ones(7)
 	end
 	Ti = [ 0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 5.0, 7.5, 10.0 ]
 	tid = findfirst(Ti .== T)
	if length(tid) > 0
		bi = [ b1[tid], b2[tid], b3[tid], b4[tid], b5[tid], b6[tid], b7[tid] ]
		return bi
	else
		print("PJS::Error, requested period 'T' does not have associated coefficients")
		return NaN*ones(7)
	end
end

# @time bi = PJSdistanceCorrectionCoefficients(0.01, "free")



"""

Function to compute the distance scaling correction
This function allows for interpolation between the specified response periods, but is strictly only calibrated
for the Ti values listed within the function.
Argument Gbranch can be one of:
  - "cy14":     returns 0.0, i.e., no distance correction
  - "1p15":     gives correction for 1.15 G branch
  - "1p00":		gives correction for 1.00 G branch
  - "free":     gives correction for free G branch
"""
function PJSdistanceCorrection(T::U, M::U, Rrup::U, Gbranch::String="free") where U<:Real
    if Gbranch == "cy14"
        return 0.0
    else
        # use interpolation if required for general application
		T = ( T < 0.01 ) ? 0.01 : ( ( T > 10.0 ) ? 10.0 : T )

        Ti = [ 0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 5.0, 7.5, 10.0 ]
		tidLO = findfirst(Ti .>= T)
		tidHI = findlast(Ti .<= T)

        if tidLO == tidHI # no interpolation required
            tid = tidLO
            bi = PJSdistanceCorrectionCoefficients(Ti[tid], Gbranch)
			rlim = 170.0

			lnΔRG = (bi[1]+bi[2]*(M-6.5))*log(max(Rrup,10.0)/10.0) +
					(bi[3]+bi[4]*(M-6.5))*(max(Rrup,10.0)-10.0) +
					bi[5]*exp(-((log(Rrup)-log(rlim))^2/(bi[6]^2))) +
					(bi[7]-bi[1])*log(sqrt( (Rrup^2  + rlim^2)/(10.0^2  + rlim^2) ))
        else # interpolation required
            TiLO = Ti[tidLO]
            TiHI = Ti[tidHI]
            Tfrac = log(T/TiLO)/log(TiHI/TiLO)
            # get values for bounding periods
            lnΔRG_LO = PJSdistanceCorrection(TiLO, M, Rrup, Gbranch)
            lnΔRG_HI = PJSdistanceCorrection(TiHI, M, Rrup, Gbranch)
            # interpolate
            lnΔRG = lnΔRG_LO + Tfrac*(lnΔRG_HI - lnΔRG_LO)
        end
        return lnΔRG
    end
end


"""
	PJSstressDropAdjustment(branch::String, ΔZtor::T) where T<:Real

Function to compute Δcm as a function of requested branch and deviation from expected Ztor
Argument 'branch' can be one of:
  - "cy14":       returns Δcm = 0
  - "lower":      gives Δcm for lower stress drop branch
  - "central":    gives Δcm for central stress drop branch
  - "upper":      gives Δcm for upper stress drop branch
"""
function PJSstressDropAdjustment(branch::String; Δσ_scale=1.4)
	if branch == "upper"
		Δσ_factor = Δσ_scale
	elseif branch == "central"
		Δσ_factor = 1.00
	elseif branch == "lower"
		Δσ_factor = 1.00/Δσ_scale
	elseif branch == "cy14"
		return 0.0
	else
		print("PJS::Error, argument 'branch' must be one of 'cy14', 'upper', 'central' or 'lower'")
		return NA
	end
	Δcm = (2.0/3.0) * log10( Δσ_factor )
	return Δcm
end


"""
	PJSsourceAdjustment(branch::String; ΔlnSa = 0.1)

Function to return Δc1 values
Argument 'branch' can be one of:
  - "cy14":       returns Δc1 = 0
  - "lower":      gives Δc1 for lower source branch
  - "central":    gives Δc1 for central source branch
  - "upper":      gives Δc1 for upper source branch
"""
function PJSsourceAdjustment(branch::String; ΔlnSa = 0.1)
	if branch == "upper"
		Δc1 = 1.28 * ΔlnSa
	elseif branch == "central"
		Δc1 = 0.0
	elseif branch == "lower"
		Δc1 = -1.28 * ΔlnSa
	elseif branch == "cy14"
		return 0.0
	else
		print("PJS::Error, argument 'branch' must be one of 'cy14', 'upper', 'central' or 'lower'")
		return NaN
	end
	return Δc1
end


"""
	PJSspainCrustalMedian(sourceBranch::String, distanceBranch::String, region::String, T::U, M::U, Rrup::U, Rjb::U, Rx::U, ΔZtor::U, Dip::U, Fnm::Int64, Frv::Int64, Fhw::Int64) where U<:Real)

Overall function for stress drop and distance corrected Spanish spectral ordinates
Key arguments that differ from CY14 are:
  - sourceBranch: which may be one of
          - "cy14":       which applies no stress drop adjustment
          - "lower":      giving the lower source/stress drop branch
          - "central":    giving the central source/stress drop branch
          - "upper":      giving the upper source/stress drop branch
  - distanceBranch: which may be one of
          - "cy14":       which applies no G or distance adjustment
          - "free":       giving the unconstrained distance branch
          - "1p00":       giving the 1.00 constrained distance branch
          - "1p15":       giving the 1.15 constrained distance branch
"""
function PJSspainCrustalMedian(sourceBranch::String, distanceBranch::String, T::U, M::U, Rrup::U, Rjb::U, Rx::U, ΔZtor::U, Dip::U, Fnm::Int64, Frv::Int64, Fhw::Int64) where U<:Real
	# Define Δc1 required for this source branch
    Δc1 = PJSsourceAdjustment(sourceBranch)
	# Define Δcm required for this source branch
    Δcm = PJSstressDropAdjustment(sourceBranch)

    # Compute the stress drop adjusted CY14 prediction of logarithmic Sa
    lnSaCY_src = PJScy2014adjusted(T, M, Rrup, Rjb, Rx, ΔZtor, Dip, Fnm, Frv, Fhw, 760.0, 1, Δc1, Δcm, 0.0, 0.0, 0.0).lnIM
    # Compute the logarithmic distance correction
    lnΔRG = PJSdistanceCorrection(T, M, Rrup, distanceBranch)
	# Compute derived source/stress drop and distance corrected logarithmic Sa
    lnSa = lnSaCY_src + lnΔRG
	return exp(lnSa)
end

"""
	PJSspainCrustalMedian(sourceBranch::String, distanceBranch::String, Ti::Vector{U}, M::U, Rrup::U, Rjb::U, Rx::U, ΔZtor::U, Dip::U, Fnm::Int64, Frv::Int64, Fhw::Int64) where U<:Real

Vectorised version of the Spanish generic crustal model that provides a spectrum
"""
function PJSspainCrustalMedian(sourceBranch::String, distanceBranch::String, Ti::Vector{U}, M::U, Rrup::U, Rjb::U, Rx::U, ΔZtor::U, Dip::U, Fnm::Int64, Frv::Int64, Fhw::Int64) where U<:Real
	Sai = similar(Ti)
	for i in 1:length(Ti)
		Sai[i] = PJSspainCrustalMedian(sourceBranch, distanceBranch, Ti[i], M, Rrup, Rjb, Rx, ΔZtor, Dip, Fnm, Frv, Fhw)
	end
	return Sai
end


"""
	PJSafricaCrustalMedian(sourceBranch::String, distanceBranch::String, region::String, T::U, M::U, Rrup::U, Rjb::U, Rx::U, ΔZtor::U, Dip::U, Fnm::Int64, Frv::Int64, Fhw::Int64) where U<:Real)

Overall function for source and distance corrected African-source spectral ordinates
Key arguments that differ from CY14 are:
  - sourceBranch: which may be one of
          - "cy14":       which applies no stress drop adjustment
          - "lower":      giving the lower source/stress drop branch
          - "central":    giving the central source/stress drop branch
          - "upper":      giving the upper source/stress drop branch
  - distanceBranch: which may be one of
          - "cy14":       which applies no G or distance adjustment
          - "lower":      giving relatively high attenuation rate
          - "central":    giving best fit attenuation rate
          - "upper":      giving relatively low attenuation rate
"""
function PJSafricaCrustalMedian(sourceBranch::String, distanceBranch::String, T::U, M::U, Rrup::U, Rjb::U, Rx::U, ΔZtor::U, Dip::U, Fnm::Int64, Frv::Int64, Fhw::Int64) where U<:Real
	# Define Δc1 required for this source branch
    Δc1 = PJSsourceAdjustment(sourceBranch)
	# Define Δcm required for this source branch
    Δcm = PJSstressDropAdjustment(sourceBranch)

	if distanceBranch == "cy14"
		Δγ1 = 0.0
	elseif distanceBranch == "lower"
		Δγ1 = -0.004957
	elseif distanceBranch == "central"
		Δγ1 = -0.002478
	elseif distanceBranch == "upper"
		Δγ1 = 0.0
	else
		print("distanceBranch must be one of 'cy14', 'lower', 'central', or 'upper'\n")
		return NaN
	end

    # Compute the stress drop adjusted CY14 prediction of logarithmic Sa
    lnSa = PJScy2014adjusted(T, M, Rrup, Rjb, Rx, ΔZtor, Dip, Fnm, Frv, Fhw, 760.0, 1, Δc1, Δcm, Δγ1, 0.0, 0.0).lnSa
	return exp(lnSa)
end

"""
	PJSafricaCrustalMedian(sourceBranch::String, distanceBranch::String, Ti::Vector{U}, M::U, Rrup::U, Rjb::U, Rx::U, ΔZtor::U, Dip::U, Fnm::Int64, Frv::Int64, Fhw::Int64) where U<:Real

Vectorised version of the African generic crustal model that provides a spectrum
"""
function PJSafricaCrustalMedian(sourceBranch::String, distanceBranch::String, Ti::Vector{U}, M::U, Rrup::U, Rjb::U, Rx::U, ΔZtor::U, Dip::U, Fnm::Int64, Frv::Int64, Fhw::Int64) where U<:Real
	Sai = similar(Ti)
	for i in 1:length(Ti)
		Sai[i] = PJSafricaCrustalMedian(sourceBranch, distanceBranch, Ti[i], M, Rrup, Rjb, Rx, ΔZtor, Dip, Fnm, Frv, Fhw)
	end
	return Sai
end


"""
	PJSspainCrustalStandardDeviation(sigmaBranch::String, T::U, M::U) where U<:Real

Standard deviation components for the mixture distribution
"""
function PJSspainCrustalStandardDeviation(sigmaBranch::String, T::U, M::U) where U<:Real

	Ti = [ 0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 5.0, 7.5, 10.0 ]

	if sigmaBranch == "central"
		# mixture component 1 coefficients
		σ_mix1_1 = [ 0.5643, 0.5635, 0.5627, 0.5618, 0.5611, 0.5591, 0.5572, 0.5535, 0.5499, 0.5431, 0.5369, 0.531, 0.5182, 0.5077, 0.4916, 0.4804, 0.4672, 0.4573, 0.4546, 0.4541 ]
		σ_mix1_2 = [ 0.5516, 0.5508, 0.5499, 0.549, 0.5482, 0.5462, 0.5443, 0.5405, 0.5368, 0.5298, 0.5234, 0.5174, 0.5043, 0.4934, 0.4768, 0.4653, 0.4515, 0.4413, 0.4385, 0.438 ]
		σ_mix1_3 = [ 0.4895, 0.4889, 0.4883, 0.4877, 0.4872, 0.4859, 0.4846, 0.4828, 0.4828, 0.4827, 0.4824, 0.4805, 0.4737, 0.4658, 0.4521, 0.4418, 0.4272, 0.4164, 0.4134, 0.4128 ]
		σ_mix1_4 = [ 0.3922, 0.3922, 0.3922, 0.3921, 0.3921, 0.392, 0.392, 0.3941, 0.4008, 0.4138, 0.4254, 0.4318, 0.4384, 0.4373, 0.4307, 0.4237, 0.4086, 0.3973, 0.3941, 0.3935 ]
		# mixture component 2 coefficients
		σ_mix2_1 = [ 0.8465, 0.8453, 0.844, 0.8427, 0.8416, 0.8387, 0.8359, 0.8302, 0.8249, 0.8147, 0.8053, 0.7965, 0.7774, 0.7615, 0.7373, 0.7206, 0.7007, 0.686, 0.6819, 0.6811 ]
		σ_mix2_2 = [ 0.8274, 0.8261, 0.8248, 0.8235, 0.8224, 0.8194, 0.8165, 0.8107, 0.8052, 0.7947, 0.7851, 0.7761, 0.7564, 0.7401, 0.7152, 0.6979, 0.6773, 0.662, 0.6577, 0.6569 ]
		σ_mix2_3 = [ 0.7342, 0.7334, 0.7325, 0.7316, 0.7308, 0.7288, 0.7269, 0.7242, 0.7242, 0.7241, 0.7235, 0.7207, 0.7106, 0.6987, 0.6781, 0.6626, 0.6409, 0.6246, 0.62, 0.6192 ]
		σ_mix2_4 = [ 0.5883, 0.5883, 0.5882, 0.5882, 0.5881, 0.588, 0.5879, 0.5911, 0.6012, 0.6206, 0.638, 0.6477, 0.6577, 0.656, 0.646, 0.6356, 0.6129, 0.5959, 0.5911, 0.5903 ]
	elseif sigmaBranch == "lower"
		# mixture component 1 coefficients
		σ_mix1_1 = [ 0.4751, 0.4745, 0.4737, 0.4729, 0.4724, 0.4707, 0.469, 0.4658, 0.4627, 0.4567, 0.4511, 0.4459, 0.4343, 0.4245, 0.4091, 0.398, 0.3844, 0.3738, 0.3708, 0.3702 ]
		σ_mix1_2 = [ 0.4595, 0.4588, 0.458, 0.4572, 0.4566, 0.4549, 0.4532, 0.4498, 0.4466, 0.4403, 0.4345, 0.4291, 0.417, 0.4067, 0.3906, 0.379, 0.3647, 0.3536, 0.3504, 0.3498 ]
		σ_mix1_3 = [ 0.3962, 0.3956, 0.395, 0.3944, 0.3939, 0.3926, 0.3913, 0.3897, 0.3901, 0.3908, 0.3911, 0.3897, 0.3838, 0.3765, 0.3642, 0.3541, 0.3388, 0.3267, 0.3233, 0.3226 ]
		σ_mix1_4 = [ 0.3107, 0.3103, 0.31, 0.3098, 0.3095, 0.3089, 0.3083, 0.3096, 0.3167, 0.3303, 0.3425, 0.3491, 0.356, 0.3561, 0.354, 0.3486, 0.3332, 0.3207, 0.317, 0.3164 ]
		# mixture component 2 coefficients
		σ_mix2_1 = [ 0.7127, 0.7117, 0.7106, 0.7094, 0.7085, 0.706, 0.7035, 0.6987, 0.694, 0.6851, 0.6767, 0.6689, 0.6515, 0.6367, 0.6136, 0.597, 0.5766, 0.5608, 0.5562, 0.5553 ]
		σ_mix2_2 = [ 0.6892, 0.6882, 0.6871, 0.6858, 0.6849, 0.6823, 0.6797, 0.6747, 0.6699, 0.6605, 0.6518, 0.6437, 0.6255, 0.6101, 0.5859, 0.5685, 0.547, 0.5304, 0.5256, 0.5247 ]
		σ_mix2_3 = [ 0.5942, 0.5934, 0.5926, 0.5916, 0.5909, 0.5889, 0.5869, 0.5845, 0.5851, 0.5862, 0.5867, 0.5846, 0.5756, 0.5648, 0.5463, 0.5311, 0.5081, 0.4901, 0.4849, 0.4839 ]
		σ_mix2_4 = [ 0.466, 0.4655, 0.4651, 0.4647, 0.4643, 0.4633, 0.4624, 0.4645, 0.475, 0.4954, 0.5137, 0.5236, 0.534, 0.5341, 0.531, 0.5229, 0.4998, 0.481, 0.4756, 0.4745 ]
	elseif sigmaBranch == "upper"
		# mixture component 1 coefficients
		σ_mix1_1 = [ 0.6586, 0.6576, 0.6566, 0.6557, 0.6548, 0.6525, 0.6504, 0.6461, 0.642, 0.6344, 0.6274, 0.6209, 0.607, 0.5957, 0.579, 0.5679, 0.5552, 0.5463, 0.5439, 0.5435 ]
		σ_mix1_2 = [ 0.6491, 0.6482, 0.6472, 0.6462, 0.6454, 0.6431, 0.6409, 0.6366, 0.6324, 0.6247, 0.6177, 0.6111, 0.597, 0.5855, 0.5686, 0.5574, 0.5445, 0.5354, 0.533, 0.5326 ]
		σ_mix1_3 = [ 0.5892, 0.5886, 0.5881, 0.5875, 0.587, 0.5856, 0.5844, 0.5825, 0.582, 0.581, 0.5798, 0.5774, 0.5698, 0.5613, 0.5462, 0.5358, 0.5224, 0.5131, 0.5106, 0.5102 ]
		σ_mix1_4 = [ 0.48, 0.4803, 0.4805, 0.4808, 0.481, 0.4816, 0.4822, 0.4851, 0.4914, 0.5034, 0.5141, 0.5202, 0.5265, 0.524, 0.5122, 0.5036, 0.489, 0.4792, 0.4766, 0.4762 ]
		# mixture component 2 coefficients
		σ_mix2_1 = [ 0.9878, 0.9864, 0.9849, 0.9835, 0.9822, 0.9788, 0.9756, 0.9692, 0.9631, 0.9516, 0.9412, 0.9314, 0.9105, 0.8936, 0.8685, 0.8519, 0.8328, 0.8194, 0.8159, 0.8153 ]
		σ_mix2_2 = [ 0.9737, 0.9723, 0.9708, 0.9694, 0.968, 0.9646, 0.9614, 0.9548, 0.9487, 0.9371, 0.9265, 0.9166, 0.8955, 0.8783, 0.8529, 0.836, 0.8167, 0.8031, 0.7995, 0.7989 ]
		σ_mix2_3 = [ 0.8838, 0.883, 0.8821, 0.8813, 0.8805, 0.8785, 0.8766, 0.8737, 0.873, 0.8715, 0.8698, 0.8661, 0.8548, 0.8419, 0.8192, 0.8037, 0.7836, 0.7696, 0.766, 0.7653 ]
		σ_mix2_4 = [ 0.72, 0.7204, 0.7208, 0.7212, 0.7216, 0.7224, 0.7233, 0.7277, 0.7371, 0.7551, 0.7711, 0.7803, 0.7898, 0.786, 0.7683, 0.7553, 0.7335, 0.7188, 0.7149, 0.7142 ]
	else
		print("PJS::Error, argument 'sigmaBranch' must be one of 'upper', 'central' or 'lower'")
		return NaN, NaN
	end

	tidLO = findlast(Ti .<= T)
	tidHI = findfirst(Ti .>= T)

	if tidLO == tidHI
		# we have a known period
		tid = tidLO
		if M <= 4.5
			σ_mix1 = σ_mix1_1[tid]
			σ_mix2 = σ_mix2_1[tid]
		elseif M <= 5.0
			σ_mix1 = σ_mix1_1[tid] + (σ_mix1_2[tid]-σ_mix1_1[tid])*(M - 4.5)/0.5
			σ_mix2 = σ_mix2_1[tid] + (σ_mix2_2[tid]-σ_mix2_1[tid])*(M - 4.5)/0.5
		elseif M <= 5.5
			σ_mix1 = σ_mix1_2[tid] + (σ_mix1_3[tid]-σ_mix1_2[tid])*(M - 5.0)/0.5
			σ_mix2 = σ_mix2_2[tid] + (σ_mix2_3[tid]-σ_mix2_2[tid])*(M - 5.0)/0.5
		elseif M <= 6.5
			σ_mix1 = σ_mix1_3[tid] + (σ_mix1_4[tid]-σ_mix1_3[tid])*(M - 5.5)/1.0
			σ_mix2 = σ_mix2_3[tid] + (σ_mix2_4[tid]-σ_mix2_3[tid])*(M - 5.5)/1.0
		else
			σ_mix1 = σ_mix1_4[tid]
			σ_mix2 = σ_mix2_4[tid]
		end
	else
		# we need to interpolate
		@inbounds TiLO = Ti[tidLO]
		@inbounds TiHI = Ti[tidHI]
		Tfrac = log(T/TiLO)/log(TiHI/TiLO)

		σ_mix1_LO, σ_mix2_LO = PJSspainCrustalStandardDeviation(sigmaBranch, TiLO, M)
		σ_mix1_HI, σ_mix2_HI = PJSspainCrustalStandardDeviation(sigmaBranch, TiHI, M)

		σ_mix1 = σ_mix1_LO + Tfrac*(σ_mix1_HI - σ_mix1_LO)
		σ_mix2 = σ_mix2_LO + Tfrac*(σ_mix2_HI - σ_mix2_LO)
	end

	return σ_mix1, σ_mix2
end


function PJSspainCrustal(T::U, rup::Rupture{U}, site::Site{U}, sourceBranch::String="lower", distanceBranch::String="central", sigmaBranch="central") where U<:Real
	M = rup.magnitude
	Rrup = rupture_distance(rup, site)
	Rjb = joyner_boore_distance(rup, site)
	Rx = strike_distance(rup, site)
	Ztor = ztor(rup)
	Dip = dip(rup)
	Fss, Fnm, Frv, Fuk = mechanism(rup)

	if Frv == 1
		EZtor = max(2.704 - 1.226*max(M-5.849,0.0),0.0)^2
	else
		EZtor = max(2.673 - 1.136*max(M-4.970,0.0),0.0)^2
	end
	ΔZtor = Ztor - EZtor
	ΔZ1p0 = 0.0
	Fhw = (Rx < 0.0) ? 0 : 1

	# get Sa
	Sa = PJSspainCrustalMedian(sourceBranch, distanceBranch, T, M, Rrup, Rjb, Rx, ΔZtor, Dip, Fnm, Frv, Fhw)
	lnSa = log(Sa)
	# get total standard deviation
	σ1, σ2 = PJSspainCrustalStandardDeviation(sigmaBranch, T, M)
	mm = MixtureModel(Normal, [(0.0, σ1), (0.0, σ2)])
	σ = sqrt( var(mm) )

	gm = PJSgroundMotion(Sa, lnSa, NaN, NaN, σ)
	return gm
end


function PJSspainCrustal(Ti::Vector{U}, rup::Rupture{U}, site::Site{U}, sourceBranch::String="lower", distanceBranch::String="central", sigmaBranch="central") where U<:Real
	predictions = Array{PJSgroundMotion,1}()
	for T in Ti
		push!(predictions, PJSspainCrustal(T, rup, site, sourceBranch, distanceBranch, sigmaBranch))
	end
	return predictions
end




"""
	PJSmixturePDF(z::T, μ::T, σ1::T, σ2::T) where T<:Real

Probability density function of a mixture of two normal distributions with the same mean μ and standard deviations of σ1 and σ2
"""
function PJSmixturePDF(z::T, μ::T, σ1::T, σ2::T) where T<:Real
	mm = MixtureModel(Normal, [(μ, σ1), (μ, σ2)])
	return pdf(mm, z)
end

"""
	PJSmixtureCDF(z::T, μ::T, σ1::T, σ2::T) where T<:Real

Cumulative distribution function of a mixture of two normal distributions with the same mean μ and standard deviations of σ1 and σ2
"""
function PJSmixtureCDF(z::T, μ::T, σ1::T, σ2::T) where T<:Real
	mm = MixtureModel(Normal, [(μ, σ1), (μ, σ2)])
	return cdf(mm, z)
end

"""
	PJSmixtureCCDF(z::T, μ::T, σ1::T, σ2::T) where T<:Real

Complementary cumulative distribution function of a mixture of two normal distributions with the same mean μ and standard deviations of σ1 and σ2
"""
function PJSmixtureCCDF(z::T, μ::T, σ1::T, σ2::T) where T<:Real
	mm = MixtureModel(Normal, [(μ, σ1), (μ, σ2)])
	return ccdf(mm, z)
end

# PJSmixturePDF(0.5, 0.0, σ1, σ2)
# PJSmixtureCDF(0.5, 0.0, σ1, σ2)
# PJSmixtureCCDF(0.5, 0.0, σ1, σ2)

"""
	PJSmixtureQuantile(p::T, μ::T, σ1::T, σ2::T) where T<:Real

Inverse CDF of a mixture of two normal distributions with the same mean μ and standard deviations of σ1 and σ2
"""
function PJSmixtureQuantile(p::T, μ::T, σ1::T, σ2::T) where T<:Real
	mm = MixtureModel(Normal, [(μ, σ1), (μ, σ2)])
	f(z) = p - cdf(mm, z)
	return find_zero(f, (-20.0, 20.0), Bisection())
end

# @time PJSmixtureQuantile(0.1, 0.0, σ1, σ2)


"""
	PJSmixtureEpsilonQuantile(ε::T, μ::T, σ1::T, σ2::T) where T<:Real

Inverse CDF of a mixture of two normal distributions with the same mean μ and standard deviations of σ1 and σ2. This function takes a value of ε as an input. From this an equivalent probability from a standard normal is computed, and the `PJSmixtureQuantile` function is called on that probability
"""
function PJSmixtureEpsilonQuantile(ε::T, μ::T, σ1::T, σ2::T) where T<:Real
	nd = Normal()
	p = cdf(nd, ε)
	return PJSmixtureQuantile(p, μ, σ1, σ2)
end

# @time PJSmixtureEpsilonQuantile(-1.0, 0.0, σ1, σ2)
# @time PJSmixtureEpsilonQuantile(1.0, 0.0, σ1, σ2)
