# File implementing the Spanish Crustal GMM
# The model uses CY14 as its base and applies a stress drop correction and a distance
# correction - both for geometric spreading and anelastic attenuation

using Distributions
using Roots

include("PJSgroundMotions.jl")

# """
# 	PJSgroundMotion{T<:Real}
#
# Custom type storing the median and logarithmic mean prediction from a ground
# motion model along with the relevant variance components
# """
# struct PJSgroundMotion{T<:Real}
# 	Sa::T
# 	lnSa::T
# 	τ::T
# 	ϕ::T
# 	σ::T
# end



"""
	PJScy2014theory( T::U, M::U, Rrup::U, Rjb::U, Rx::U, ΔZtor::U, Dip::U, Fnm::Int64, Frv::Int64, Fhw::Int64, Vs30::U, Vs30meas::Int64, Δcm::U=0.0, ΔZ1p0::U=0.0, ΔDPP::U=0.0 ) where U<:Real

Modification to the Chiou & Youngs (2014) gmm to account for theoretical stress
drop adjustments. Note that all parameters are as in the original model with the
exception of Δcm which is the change to the cm coefficient required to account for
a particular level of stress drop.
"""
function PJScy2014theory( T::U, M::U, Rrup::U, Rjb::U, Rx::U, ΔZtor::U, Dip::U, Fnm::Int64, Frv::Int64, Fhw::Int64, Vs30::U, Vs30meas::Int64, Δcm::U=0.0, ΔZ1p0::U=0.0, ΔDPP::U=0.0 ) where U<:Real
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

		@fastmath @inbounds lnYref = c1[tid] +
			( c1a[tid] + c1c[tid]/cosh(2.0*max(M-4.5,0.0)) )*Frv +
			( c1b[tid] + c1d[tid]/cosh(2.0*max(M-4.5,0.0)) )*Fnm +
			( c7[tid]  + c7b[tid]/cosh(2.0*max(M-4.5,0.0)) )*ΔZtor +
			( c11[tid] + c11b[tid]/cosh(2.0*max(M-4.5,0.0)) )*(cosd(Dip)^2) +
			c2*(M-6.0) + ((c2-c3[tid])/cn[tid])*log(1.0 + exp( cn[tid]*(cm[tid]+(Δcm-δm)-M) )) - (c2-c3[tid])*(Δcm-δm) +
			c4*log( Rrup + c5[tid]*cosh( c6[tid]*max(M-chm[tid],0.0) ) ) +
			(c4a-c4)*log( sqrt( Rrup^2 + crb^2 ) ) +
			( γ1[tid] + γ2[tid]/cosh(max(M-γ3[tid],0.0)) )*Rrup +
			c8[tid]*max( 1.0 - max(Rrup-40.0,0.0)/30.0, 0.0 )*min( max(M-5.5,0.0)/0.8, 1.0 )*exp( -c8a*(M-c8b[tid])^2 )* ΔDPP +
			c9[tid]*Fhw*cosd(Dip)*( c9a[tid] + (1.0-c9a[tid])*tanh( Rx / c9b[tid] ) )*( 1.0 - sqrt(Rjb^2+Ztor^2)/(Rrup+1.0) )
		@fastmath @inbounds lnSa = lnYref +
			ϕ1[tid]*min(log(Vs30/1130.0),0.0) +
			ϕ2[tid]*( exp(ϕ3[tid]*(min(Vs30,1130.0)-360.0)) - exp(ϕ3[tid]*(1130.0-360.0)))*log( ( exp(lnYref) + ϕ4[tid] ) / ϕ4[tid] ) +
			ϕ5[tid]*( 1.0 - exp( -ΔZ1p0/ϕ6[tid] ) )
		Sa = exp( lnSa )

		@fastmath @inbounds NL = ϕ2[tid]*( exp(ϕ3[tid]*(min(Vs30,1130.0)-360.0)) - exp(ϕ3[tid]*(1130.0-360.0))) * ( exp(lnYref) / ( exp(lnYref) + ϕ4[tid] ) )
		@inbounds τ = ( τ1[tid] + ( τ2[tid] - τ1[tid] )*( min( max( M, 5.0 ), 7.25 ) - 5.0 ) / 2.25 ) * (1.0+NL)
		@inbounds ϕ = ( σ1[tid] + ( σ2[tid] - σ1[tid] )*( min( max( M, 5.0 ), 7.25 ) - 5.0 ) / 2.25 ) * sqrt( ifelse( Vs30meas==1, 0.7, σ3[tid] ) + ( 1.0 + NL )^2 )
		σ = sqrt( τ*τ + ϕ*ϕ )

	else
		# we need to interpolate
		@inbounds TiLO = Ti[tidLO]
		@inbounds TiHI = Ti[tidHI]
		Tfrac = log(T/TiLO)/log(TiHI/TiLO)

		cyLO = PJScy2014theory( TiLO, M, Rrup, Rjb, Rx, ΔZtor, Dip, Fnm, Frv, Fhw, Vs30, Vs30meas, Δcm, ΔZ1p0, ΔDPP )
		cyHI = PJScy2014theory( TiHI, M, Rrup, Rjb, Rx, ΔZtor, Dip, Fnm, Frv, Fhw, Vs30, Vs30meas, Δcm, ΔZ1p0, ΔDPP )

		lnSa = cyLO.lnIM + Tfrac*(cyHI.lnIM-cyLO.lnIM)
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
function PJSdistanceCorrectionCoefficients(T::U, Qbranch="central") where U<:Real
    if Qbranch == "central"
        b1 = [ 0.2335, 0.2203, 0.2058, 0.2176, 0.2378, 0.2568, 0.2518, 0.2412, 0.2366, 0.2332, 0.232, 0.2314, 0.2308, 0.2306, 0.2306, 0.2307, 0.2307, 0.2305, 0.2306, 0.2308 ]
        b2 = [ 0.9217, 0.9377, 0.9764, 0.9895, 0.9748, 0.9498, 0.9208, 0.8847, 0.8704, 0.8602, 0.8566, 0.8549, 0.853, 0.8523, 0.8518, 0.8516, 0.8517, 0.8522, 0.8531, 0.8539 ]
        b3 = [ 0.6468, 0.6813, 0.7933, 0.8892, 0.9072, 0.7878, 0.6545, 0.5186, 0.4684, 0.4338, 0.4217, 0.4159, 0.4096, 0.4071, 0.405, 0.4043, 0.4039, 0.404, 0.4041, 0.4039 ]
        b4 = [ 0.0012, 0.0009928, 0.0004364, -0.0004544, -0.001477, -0.001594, 0.0000568, 0.00223, 0.00311, 0.003688, 0.003822, 0.003836, 0.003742, 0.003623, 0.003421, 0.003271, 0.003083, 0.002968, 0.003009, 0.003104 ]
        b5 = [ -0.000009817, 0.000005883, 0.00002577, 0.00003062, 0.00003156, 0.00002564, 0.00001679, 0.00000469, -0.000001007, -0.000006022, -0.000008621, -0.00001056, -0.00001442, -0.0000178, -0.00002382, -0.00002952, -0.00004248, -0.00007701, -0.0001217, -0.0001616 ]
        b6 = [ 0.005105, 0.005318, 0.005851, 0.006691, 0.007655, 0.008045, 0.006651, 0.004386, 0.003206, 0.002096, 0.001573, 0.001268, 0.000863, 0.0006582, 0.0004516, 0.000353, 0.0002616, 0.0001749, 0.00009157, 0.00001708 ]
        b7 = [ -0.0003478, -0.0003676, -0.0003971, -0.000415, -0.0004248, -0.0003756, -0.0002823, -0.0001612, -0.0001082, -0.00006643, -0.00004979, -0.00004092, -0.00002931, -0.000023, -0.00001611, -0.00001293, -0.000009795, -0.000002994, 0.000008319, 0.00002011 ]
        b8 = [ -0.4839, -0.4872, -0.5209, -0.5199, -0.4355, -0.3388, -0.4107, -0.5275, -0.5764, -0.6117, -0.6244, -0.6306, -0.6377, -0.6407, -0.6433, -0.6438, -0.6424, -0.637, -0.6311, -0.6263 ]
    elseif Qbranch == "lower"
        b1 = [ 0.2116, 0.1912, 0.1563, 0.1579, 0.1862, 0.2317, 0.2406, 0.238, 0.2351, 0.2326, 0.2317, 0.2312, 0.2307, 0.2305, 0.2305, 0.2306, 0.2305, 0.2303, 0.2303, 0.2305 ]
        b2 = [ 0.9439, 0.9573, 1.02, 1.081, 1.07, 0.9626, 0.9202, 0.8848, 0.8707, 0.8605, 0.8568, 0.8551, 0.8531, 0.8524, 0.8518, 0.8517, 0.8518, 0.8524, 0.8533, 0.8541 ]
        b3 = [ 0.7157, 0.7412, 0.8695, 1.053, 1.15, 0.9501, 0.7218, 0.5365, 0.4763, 0.4368, 0.4235, 0.4172, 0.4104, 0.4076, 0.4054, 0.4048, 0.4047, 0.4055, 0.406, 0.4061 ]
        b4 = [ -0.0001268, -0.0003013, -0.0006934, -0.001396, -0.002427, -0.003496, -0.002164, 0.000347, 0.001602, 0.002631, 0.003009, 0.003171, 0.003278, 0.003258, 0.003147, 0.003031, 0.002849, 0.002692, 0.002682, 0.00274 ]
        b5 = [ 0.00009378, 0.0001089, 0.0001348, 0.0001571, 0.0001758, 0.0001507, 0.00009084, 0.00003564, 0.00001764, 0.00000522, 0.0000003158, -0.000002578, -0.000007123, -0.0000103, -0.00001479, -0.00001796, -0.00002454, -0.00004733, -0.00008213, -0.0001154 ]
        b6 = [ 0.006032, 0.006209, 0.006603, 0.007302, 0.008254, 0.009274, 0.00814, 0.005703, 0.004278, 0.002855, 0.00216, 0.001748, 0.0012, 0.0009258, 0.0006575, 0.0005405, 0.0004585, 0.0004264, 0.0003973, 0.0003603 ]
        b7 = [ -0.0004249, -0.0004442, -0.0004822, -0.0005207, -0.0005431, -0.000448, -0.000314, -0.0001748, -0.0001177, -0.00007316, -0.0000556, -0.00004638, -0.0000347, -0.00002883, -0.00002369, -0.00002325, -0.00002697, -0.00003235, -0.00003088, -0.0000255 ]
        b8 = [ -0.5402, -0.526, -0.559, -0.6462, -0.6612, -0.5056, -0.4863, -0.5479, -0.585, -0.6148, -0.626, -0.6317, -0.6382, -0.6411, -0.6436, -0.6442, -0.6429, -0.6377, -0.632, -0.6273 ]
    elseif Qbranch == "upper"
        b1 = [ 0.24, 0.2301, 0.2223, 0.2353, 0.2514, 0.2617, 0.2538, 0.2418, 0.2369, 0.2333, 0.2321, 0.2315, 0.2308, 0.2307, 0.2307, 0.2307, 0.2307, 0.2306, 0.2306, 0.2309 ]
        b2 = [ 0.9141, 0.9283, 0.9555, 0.9596, 0.9562, 0.9475, 0.92, 0.8844, 0.8702, 0.8601, 0.8566, 0.8549, 0.853, 0.8523, 0.8517, 0.8516, 0.8516, 0.8522, 0.853, 0.8538 ]
        b3 = [ 0.62, 0.6542, 0.7464, 0.8084, 0.8207, 0.7507, 0.6398, 0.5144, 0.4665, 0.433, 0.4212, 0.4156, 0.4094, 0.4069, 0.4049, 0.4041, 0.4036, 0.4035, 0.4034, 0.4031 ]
        b4 = [ 0.001769, 0.001543, 0.0009177, 0.00003917, -0.0008127, -0.0006111, 0.0009678, 0.002909, 0.003639, 0.004053, 0.004102, 0.004064, 0.003902, 0.003749, 0.003517, 0.003357, 0.003169, 0.003072, 0.003131, 0.00324 ]
        b5 = [ -0.0000456, -0.00003008, -0.00001438, -0.00001785, -0.00001974, -0.000005599, -0.00000003897, -0.000002906, -0.00000599, -0.000009355, -0.00001143, -0.00001315, -0.00001693, -0.00002046, -0.00002717, -0.00003394, -0.00004942, -0.00008831, -0.0001365, -0.0001788 ]
        b6 = [ 0.004655, 0.004888, 0.005485, 0.006329, 0.007178, 0.00729, 0.005926, 0.003834, 0.002774, 0.001796, 0.001343, 0.001079, 0.0007309, 0.0005533, 0.0003702, 0.0002781, 0.0001818, 0.00007339, -0.00003073, -0.0001195 ]
        b7 = [ -0.0003177, -0.0003372, -0.0003614, -0.0003719, -0.0003827, -0.0003537, -0.0002707, -0.0001557, -0.0001044, -0.00006381, -0.00004753, -0.00003879, -0.00002719, -0.00002069, -0.00001304, -0.000008684, -0.000002747, 0.00000881, 0.00002387, 0.00003808 ]
        b8 = [ -0.4646, -0.4709, -0.4902, -0.4541, -0.3622, -0.3089, -0.3991, -0.5243, -0.5749, -0.6111, -0.6241, -0.6304, -0.6375, -0.6406, -0.6432, -0.6438, -0.6423, -0.6368, -0.6309, -0.6261 ]
    else
        print("PJS::Error, argument 'Qbranch' must be one of 'central', 'lower', or 'upper'")
        return NaN*ones(8)
    end
	Ti = [ 0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 5.0, 7.5, 10.0 ]
	tid = findfirst(Ti .== T)
    if length(tid) > 0
		bi = [ b1[tid], b2[tid], b3[tid], b4[tid], b5[tid], b6[tid], b7[tid], b8[tid] ]
        return bi
    else
        print("PJS::Error, requested period 'T' does not have associated coefficients")
        return NaN*ones(8)
    end
end


"""
	PJScy2014anelastic(T::U, M::U, Rrup::U) where U<:Real

Computes the anelastic attenuation contribution from the CY14 model
"""
function PJScy2014anelastic(T::U, M::U, Rrup::U) where U<:Real
    Ti = [ -1.0, 0.0, 0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.1, 0.12, 0.15, 0.17, 0.2, 0.25, 0.3, 0.4, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 7.5, 10.0 ]
    γ1 = [ -0.001852, -0.007146, -0.007146, -0.007249, -0.007869, -0.008316, -0.008743, -0.009537, -0.00983, -0.009913, -0.009896, -0.009787, -0.009505, -0.008918, -0.008251, -0.007267, -0.006492, -0.005147, -0.004277, -0.002979, -0.002301, -0.001344, -0.001084, -0.00101, -0.000964, -0.00095 ]
    γ2 = [ -0.007403, -0.006758, -0.006758, -0.006758, -0.006758, -0.006758, -0.006758, -0.00619, -0.005332, -0.004732, -0.003806, -0.00328, -0.00269, -0.002128, -0.001812, -0.001274, -0.001074, -0.001115, -0.001197, -0.001675, -0.002349, -0.003306, -0.003566, -0.00364, -0.003686, -0.0037 ]
    γ3 = [ 4.3439, 4.2542, 4.2542, 4.2386, 4.2519, 4.296, 4.3578, 4.5455, 4.7603, 4.8963, 5.0644, 5.1371, 5.188, 5.2164, 5.1954, 5.0899, 4.7854, 4.3304, 4.1667, 4.0029, 3.8949, 3.7928, 3.7443, 3.709, 3.6632, 3.623 ]

    # check bounds on the allowable value of T
    if T < 0.0 # PGV requested
        tidHI = 1
        tidLO = 1
    else
        # limit the period range (note that coefficients for 0.01 and 0.0 are equal)
        if T < 0.01
            T = 0.01
        elseif T > 10.0
            T = 10.0
        end
		tidLO = findfirst(Ti .>= T)
		tidHI = findlast(Ti .<= T)
    end

    if tidLO == tidHI
        # we have a known period
        tid = tidHI

        lnQ = ( γ1[tid] + γ2[tid]/cosh(max(M-γ3[tid],0.0)) )*Rrup
    else
        # we need to interpolate
        TiLO = Ti[tidLO]
        TiHI = Ti[tidHI]
        Tfrac = log(T/TiLO)/log(TiHI/TiLO)

        lnQLO = PJScy2014anelastic( TiLO, M, Rrup )
        lnQHI = PJScy2014anelastic( TiHI, M, Rrup )

        lnQ = lnQLO + Tfrac*(lnQHI-lnQLO)
    end
    return lnQ
end


"""

Function to compute the distance scaling correction
This function allows for interpolation between the specified response periods, but is strictly only calibrated
for the Ti values listed within the function.
Argument Qbranch can be one of:
  - "cy14":       returns 0.0, i.e., no distance correction
  - "lower":      gives correction for lower Q branch
  - "central":    gives correction for central Q branch
gives correction for upper Q branch
"""
function PJSdistanceCorrection(T::U, M::U, Rrup::U, Qbranch::String="central") where U<:Real
    if Qbranch == "cy14"
        return 0.0
    else
        # use interpolation if required for general application
		T = ( T < 0.01 ) ? 0.01 : ( ( T > 10.0 ) ? 10.0 : T )

        Ti = [ 0.01, 0.02, 0.03, 0.04, 0.05, 0.075, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 5.0, 7.5, 10.0 ]
		tidLO = findfirst(Ti .>= T)
		tidHI = findlast(Ti .<= T)

        if tidLO == tidHI # no interpolation required
            tid = tidLO
            bi = PJSdistanceCorrectionCoefficients(Ti[tid], Qbranch)
            lnΔQ0 = (bi[4] + bi[5]*M)*max(Rrup-10.0, 0.0) + (bi[6] + bi[7]*M)*max(Rrup-200.0, 0.0)
            lnQcy = PJScy2014anelastic(Ti[tid], M, Rrup)
            lnΔQ = min( -lnQcy, lnΔQ0 )
            lnΔRQ = bi[1]*log(max(Rrup,10.0)/10.0) + bi[2]*log(max(min(Rrup,130.0),70.0)/70.0) + bi[3]*log(max(Rrup,130.0)/130.0) + lnΔQ + bi[8]*log(sqrt(max(Rrup,10.0)^2 + 50.0^2)/sqrt(10.0^2 + 50.0^2))
        else # interpolation required
            TiLO = Ti[tidLO]
            TiHI = Ti[tidHI]
            Tfrac = log(T/TiLO)/log(TiHI/TiLO)
            # get values for bounding periods
            lnΔRQ_LO = PJSdistanceCorrection(TiLO, M, Rrup, Qbranch)
            lnΔRQ_HI = PJSdistanceCorrection(TiHI, M, Rrup, Qbranch)
            # interpolate
            lnΔRQ = lnΔRQ_LO + Tfrac*(lnΔRQ_HI - lnΔRQ_LO)
        end
        return lnΔRQ
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
function PJSstressDropAdjustment(branch::String, ΔZtor::T) where T<:Real
	if branch == "upper"
		Δσ0 = 66.0
		δΔσ = 4.0
	elseif branch == "central"
		Δσ0 = 40.0
		δΔσ = 1.4
	elseif branch == "lower"
		Δσ0 = 26.0
		δΔσ = 0.0
	elseif branch == "cy14"
	    # Δσ0 = 26.4
	    # δΔσ = 5.35
		return 0.0
	else
		print("PJS::Error, argument 'branch' must be one of 'cy14', 'upper', 'central' or 'lower'")
		return NA
	end
	Δcm = (2.0/3.0) * log10( (Δσ0 + δΔσ*ΔZtor) / max(26.4 + 5.35*ΔZtor, 5.0) )
	return Δcm
end


"""
	PJSspainCrustalMedian(stressDropBranch::String, distanceBranch::String, T::U, M::U, Rrup::U, Rjb::U, Rx::U, ΔZtor::U, Dip::U, Fnm::Int64, Frv::Int64, Fhw::Int64) where U<:Real)

Overall function for stress drop and distance corrected Spanish spectral ordinates
Key arguments that differ from CY14 are:
  - stressDropBranch: which may be one of
          - "cy14":       which applies no stress drop adjustment
          - "lower":      giving the lower stress drop branch
          - "central":    giving the central stress drop branch
          - "upper":      giving the upper stress drop branch
  - distanceBranch: which may be one of
          - "cy14":       which applies no Q or distance adjustment
          - "lower":      giving the lower Q/distance branch
          - "central":    giving the central Q/distance branch
          - "upper":      giving the upper Q/distance branch

"""
function PJSspainCrustalMedian(stressDropBranch::String, distanceBranch::String, T::U, M::U, Rrup::U, Rjb::U, Rx::U, ΔZtor::U, Dip::U, Fnm::Int64, Frv::Int64, Fhw::Int64) where U<:Real
	# Define Δcm required for this stress drop branch
    Δcm = PJSstressDropAdjustment(stressDropBranch, ΔZtor)
    # Compute the stress drop adjusted CY14 prediction of logarithmic Sa
    lnSaCY_Δσ = PJScy2014theory( T, M, Rrup, Rjb, Rx, ΔZtor, Dip, Fnm, Frv, Fhw, 760.0, 1, Δcm, 0.0, 0.0 ).lnIM
    # Compute the logarithmic distance correction
    lnΔRQ = PJSdistanceCorrection(T, M, Rrup, distanceBranch)
	# Compute derived stress drop and distance corrected logarithmic Sa
    lnSa = lnSaCY_Δσ + lnΔRQ
	return exp(lnSa)
end

"""
	PJSspainCrustalMedian(stressDropBranch::String, distanceBranch::String, Ti::Vector{U}, M::U, Rrup::U, Rjb::U, Rx::U, ΔZtor::U, Dip::U, Fnm::Int64, Frv::Int64, Fhw::Int64) where U<:Real

Vectorised version of the Spanish generic crustal model that provides a spectrum
"""
function PJSspainCrustalMedian(stressDropBranch::String, distanceBranch::String, Ti::Vector{U}, M::U, Rrup::U, Rjb::U, Rx::U, ΔZtor::U, Dip::U, Fnm::Int64, Frv::Int64, Fhw::Int64) where U<:Real
	Sai = similar(Ti)
	for i in 1:length(Ti)
		Sai[i] = PJSspainCrustalMedian(stressDropBranch, distanceBranch, Ti[i], M, Rrup, Rjb, Rx, ΔZtor, Dip, Fnm, Frv, Fhw)
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
		σ_mix1_1 = [ 0.6228, 0.6221, 0.6213, 0.6205, 0.6199, 0.6181, 0.6165, 0.6131, 0.6099, 0.6039, 0.5984, 0.5932, 0.582, 0.5726, 0.5586, 0.5488, 0.5373, 0.5288, 0.5264, 0.526 ]
		σ_mix1_2 = [ 0.6045, 0.6037, 0.603, 0.6022, 0.6015, 0.5997, 0.598, 0.5945, 0.5912, 0.585, 0.5793, 0.5739, 0.5623, 0.5526, 0.538, 0.5279, 0.5159, 0.5069, 0.5045, 0.504 ]
		σ_mix1_3 = [ 0.5373, 0.5368, 0.5363, 0.5357, 0.5353, 0.5341, 0.5329, 0.5313, 0.5314, 0.5313, 0.5311, 0.5294, 0.5234, 0.5164, 0.5042, 0.4951, 0.4822, 0.4726, 0.4699, 0.4694 ]
		σ_mix1_4 = [ 0.4415, 0.4415, 0.4414, 0.4414, 0.4413, 0.4412, 0.4411, 0.4429, 0.4489, 0.4604, 0.4709, 0.4766, 0.4827, 0.4818, 0.4762, 0.4701, 0.4566, 0.4464, 0.4436, 0.4431 ]
		# mixture component 2 coefficients
		σ_mix2_1 = [ 0.7943, 0.793, 0.7917, 0.7903, 0.7891, 0.7859, 0.7829, 0.7768, 0.7711, 0.7601, 0.75, 0.7405, 0.7198, 0.7025, 0.6762, 0.6579, 0.636, 0.6197, 0.6152, 0.6144 ]
		σ_mix2_2 = [ 0.7803, 0.779, 0.7776, 0.7762, 0.775, 0.7718, 0.7686, 0.7625, 0.7566, 0.7454, 0.7351, 0.7254, 0.7042, 0.6866, 0.6596, 0.6408, 0.6183, 0.6015, 0.5968, 0.596 ]
		σ_mix2_3 = [ 0.6914, 0.6905, 0.6896, 0.6887, 0.6879, 0.6857, 0.6836, 0.6808, 0.6808, 0.6805, 0.6799, 0.6769, 0.666, 0.6533, 0.6311, 0.6144, 0.5908, 0.5731, 0.5682, 0.5673 ]
		σ_mix2_4 = [ 0.5435, 0.5435, 0.5435, 0.5434, 0.5434, 0.5433, 0.5433, 0.5468, 0.5577, 0.5786, 0.5971, 0.6074, 0.6181, 0.6162, 0.6052, 0.5939, 0.5696, 0.5513, 0.5462, 0.5453 ]
	elseif sigmaBranch == "lower"
		# mixture component 1 coefficients
		σ_mix1_1 = [ 0.4982, 0.4977, 0.4972, 0.4965, 0.4961, 0.4948, 0.4935, 0.4911, 0.4887, 0.4841, 0.4797, 0.4757, 0.4666, 0.4587, 0.4463, 0.4373, 0.426, 0.4171, 0.4145, 0.414 ]
		σ_mix1_2 = [ 0.4752, 0.4747, 0.4741, 0.4735, 0.473, 0.4717, 0.4703, 0.4677, 0.4652, 0.4603, 0.4557, 0.4514, 0.4417, 0.4334, 0.4202, 0.4106, 0.3986, 0.3891, 0.3863, 0.3858 ]
		σ_mix1_3 = [ 0.407, 0.4065, 0.406, 0.4055, 0.4051, 0.404, 0.4028, 0.4015, 0.4022, 0.4032, 0.404, 0.403, 0.3984, 0.3927, 0.3832, 0.3748, 0.3617, 0.3512, 0.3481, 0.3475 ]
		σ_mix1_4 = [ 0.3302, 0.3297, 0.3294, 0.329, 0.3286, 0.3278, 0.3269, 0.3276, 0.3336, 0.3453, 0.356, 0.3617, 0.3677, 0.3686, 0.3696, 0.366, 0.3527, 0.3415, 0.3382, 0.3376 ]
		# mixture component 2 coefficients
		σ_mix2_1 = [ 0.6944, 0.6932, 0.692, 0.6906, 0.6896, 0.6867, 0.6839, 0.6783, 0.673, 0.6628, 0.6532, 0.6443, 0.6245, 0.6078, 0.5817, 0.5631, 0.5403, 0.5227, 0.5177, 0.5168 ]
		σ_mix2_2 = [ 0.6775, 0.6763, 0.675, 0.6737, 0.6726, 0.6697, 0.6667, 0.661, 0.6555, 0.645, 0.6352, 0.626, 0.6056, 0.5883, 0.5613, 0.5419, 0.5182, 0.4999, 0.4947, 0.4937 ]
		σ_mix2_3 = [ 0.5869, 0.586, 0.5851, 0.5841, 0.5833, 0.5811, 0.579, 0.5763, 0.5767, 0.5774, 0.5776, 0.575, 0.5648, 0.5526, 0.5315, 0.5147, 0.4896, 0.4701, 0.4646, 0.4635 ]
		σ_mix2_4 = [ 0.4504, 0.45, 0.4497, 0.4494, 0.4491, 0.4484, 0.4477, 0.4504, 0.462, 0.484, 0.5038, 0.5144, 0.5256, 0.525, 0.519, 0.5093, 0.4842, 0.4641, 0.4584, 0.4573 ]
	elseif sigmaBranch == "upper"
		# mixture component 1 coefficients
		σ_mix1_1 = [ 0.7564, 0.7555, 0.7545, 0.7536, 0.7527, 0.7504, 0.7483, 0.744, 0.74, 0.7324, 0.7256, 0.7192, 0.7056, 0.6948, 0.679, 0.6686, 0.657, 0.6491, 0.6471, 0.6467 ]
		σ_mix1_2 = [ 0.7438, 0.7429, 0.7418, 0.7409, 0.74, 0.7377, 0.7356, 0.7312, 0.7271, 0.7195, 0.7126, 0.7061, 0.6923, 0.6813, 0.6653, 0.6548, 0.643, 0.635, 0.6329, 0.6326 ]
		σ_mix1_3 = [ 0.6794, 0.6789, 0.6783, 0.6778, 0.6773, 0.676, 0.6748, 0.6729, 0.6723, 0.6709, 0.6695, 0.667, 0.6595, 0.651, 0.6361, 0.6262, 0.614, 0.6057, 0.6036, 0.6033 ]
		σ_mix1_4 = [ 0.5634, 0.5638, 0.5641, 0.5645, 0.5648, 0.5656, 0.5664, 0.5695, 0.5753, 0.5864, 0.5962, 0.6019, 0.6078, 0.6049, 0.5916, 0.5827, 0.5691, 0.5604, 0.5582, 0.5578 ]
		# mixture component 2 coefficients
		σ_mix2_1 = [ 0.8987, 0.8973, 0.8957, 0.8943, 0.893, 0.8896, 0.8863, 0.8797, 0.8735, 0.8618, 0.8511, 0.841, 0.8194, 0.8018, 0.7754, 0.7576, 0.7369, 0.7222, 0.7182, 0.7175 ]
		σ_mix2_2 = [ 0.8878, 0.8864, 0.8848, 0.8834, 0.8821, 0.8786, 0.8753, 0.8687, 0.8624, 0.8506, 0.8397, 0.8296, 0.8077, 0.7898, 0.7631, 0.7451, 0.7242, 0.7092, 0.7052, 0.7045 ]
		σ_mix2_3 = [ 0.8015, 0.8007, 0.7997, 0.7989, 0.798, 0.796, 0.794, 0.791, 0.7904, 0.7892, 0.7878, 0.7842, 0.7726, 0.7594, 0.7362, 0.7198, 0.6982, 0.6828, 0.6787, 0.678 ]
		σ_mix2_4 = [ 0.6423, 0.6427, 0.643, 0.6433, 0.6436, 0.6442, 0.6449, 0.6492, 0.6593, 0.6786, 0.6957, 0.7055, 0.7155, 0.7121, 0.6957, 0.6828, 0.6596, 0.6434, 0.6391, 0.6383 ]
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



function PJSspainCrustal(T::U, rup::Rupture{U}, site::Site{U}, stressDropBranch::String="lower", distanceBranch::String="central", sigmaBranch="central") where U<:Real
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
	Sa = PJSspainCrustalMedian(stressDropBranch, distanceBranch, T, M, Rrup, Rjb, Rx, ΔZtor, Dip, Fnm, Frv, Fhw)
	lnSa = log(Sa)
	# get total standard deviation
	σ1, σ2 = PJSspainCrustalStandardDeviation(sigmaBranch, T, M)
	mm = MixtureModel(Normal, [(0.0, σ1), (0.0, σ2)])
	σ = sqrt( var(mm) )

	gm = PJSgroundMotion(Sa, lnSa, NaN, NaN, σ)
	return gm
end




function PJSspainCrustal(Ti::Vector{U}, rup::Rupture{U}, site::Site{U}, stressDropBranch::String="lower", distanceBranch::String="central", sigmaBranch="central") where U<:Real
	predictions = Array{PJSgroundMotion,1}()
	for T in Ti
		push!(predictions, PJSspainCrustal(T, rup, site, stressDropBranch, distanceBranch, sigmaBranch))
	end
	return predictions
end




# σ1, σ2 = PJSspainCrustalStandardDeviation("central", 0.01, 6.0)
#
# mm = MixtureModel(Normal, [(0.0, σ1), (0.0, σ2)])
#
# pdf(mm, -1.0)
# @time cdf(mm, -1.0)
# quantile(mm, 0.05)
#
# n1 = Normal(0.0, σ1)
# n2 = Normal(0.0, σ2)
#
# 0.5*( pdf(n1, -1.0) + pdf(n2, -1.0) )
# @time 0.5*( cdf(n1, -1.0) + cdf(n2, -1.0) )



# p = 0.15
# f(z) = p - cdf(mm, z)
#
# @time find_zero(f, (-10.0, 10.0), Bisection())
# @time find_zero(f, (-15.0, 15.0), Bisection())
#
# quantile(n1, p)
# quantile(n2, p)
#
# 0.5*( quantile(n1, p) + quantile(n2, p) )


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
