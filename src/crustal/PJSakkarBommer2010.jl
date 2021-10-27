

function PJSab2010(T::U, M::U, Rjb::U, Sstiff::Int64, Ssoft::Int64, Fnm::Int64, Frv::Int64 ) where U<:Real

	# make sure that both site indicators are not turned on
	if Sstiff + Ssoft > 1
		return PJSgroundMotion()
	end
	# make sure that both normal and reverse indicators are not turned on
	if Fnm + Frv > 1
		return PJSgroundMotion()
	end

    if T > 3.0
        return PJSgroundMotion()
	elseif T < 0.05
		T = 0.0
    end

    # common period vector
    Ti = [  0.0, 0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 1.0, 1.05, 1.1, 1.15, 1.2, 1.25, 1.3, 1.35, 1.4, 1.45, 1.5, 1.55, 1.6, 1.65, 1.7, 1.75, 1.8, 1.85, 1.9, 1.95, 2.0, 2.05, 2.1, 2.15, 2.2, 2.25, 2.3, 2.35, 2.4, 2.45, 2.5, 2.55, 2.6, 2.65, 2.7, 2.75, 2.8, 2.85, 2.9, 2.95, 3.0 ]

    # period-dependent coefficients
    b1 = [  1.04159, 2.11528, 2.11994, 1.64489, 0.92065, 0.13978, -0.84006, -1.32207, -1.7032, -1.97201, -2.76925, -3.51672, -3.92759, -4.4949, -4.62925, -4.95053, -5.32863, -5.75799, -5.82689, -5.90592, -6.17066, -6.60337, -6.90379, -6.9618, -6.99236, -6.74613, -6.51719, -6.55821, -6.61945, -6.62737, -6.71787, -6.80776, -6.83632, -6.88684, -6.946, -7.09166, -7.22818, -7.29772, -7.35522, -7.40716, -7.50404, -7.55598, -7.53463, -7.50811, -8.09168, -8.11057, -8.16272, -7.94704, -7.96679, -7.97878, -7.88403, -7.68101, -7.72574, -7.53288, -7.41587, -7.34541, -7.24561, -7.07107, -6.99332, -6.95669, -6.92924 ]
    b2 = [  0.91333, 0.72571, 0.75179, 0.83683, 0.96815, 1.13068, 1.37439, 1.47055, 1.5593, 1.61645, 1.83268, 2.02523, 2.08471, 2.21154, 2.21764, 2.29142, 2.38389, 2.50635, 2.50287, 2.51405, 2.58558, 2.69584, 2.77044, 2.75857, 2.73427, 2.62375, 2.51869, 2.52238, 2.52611, 2.49858, 2.49486, 2.50291, 2.51009, 2.54048, 2.57151, 2.62938, 2.66824, 2.67565, 2.67749, 2.68206, 2.71004, 2.72737, 2.71709, 2.71035, 2.91159, 2.92087, 2.93325, 2.85328, 2.85363, 2.849, 2.81817, 2.7572, 2.82043, 2.74824, 2.69012, 2.65352, 2.61028, 2.56123, 2.52699, 2.51006, 2.45899 ]
    b3 = [  -0.0814, -0.07351, -0.07448, -0.07544, -0.07903, -0.08761, -0.10349, -0.10873, -0.11388, -0.11742, -0.13202, -0.14495, -0.14648, -0.15522, -0.15491, -0.15983, -0.16571, -0.17479, -0.17367, -0.17417, -0.17938, -0.18646, -0.19171, -0.1889, -0.18491, -0.17392, -0.1633, -0.16307, -0.16274, -0.1591, -0.15689, -0.15629, -0.15676, -0.15995, -0.16294, -0.16794, -0.17057, -0.17004, -0.16934, -0.16906, -0.1713, -0.17291, -0.17221, -0.17212, -0.1892, -0.19044, -0.19155, -0.18539, -0.18561, -0.18527, -0.1832, -0.17905, -0.18717, -0.18142, -0.17632, -0.17313, -0.16951, -0.16616, -0.16303, -0.16142, -0.15513 ]
    b4 = [  -2.92728, -3.33201, -3.10538, -2.75848, -2.49264, -2.33824, -2.19123, -2.12993, -2.12718, -2.16619, -2.12969, -2.04211, -1.88144, -1.79031, -1.798, -1.81321, -1.77273, -1.77068, -1.76295, -1.79854, -1.80717, -1.73843, -1.71109, -1.66588, -1.5912, -1.52886, -1.46527, -1.48223, -1.48257, -1.4331, -1.35301, -1.31227, -1.3326, -1.40931, -1.47676, -1.54037, -1.54273, -1.50936, -1.46988, -1.43816, -1.44395, -1.45794, -1.46662, -1.49679, -1.55644, -1.59537, -1.60461, -1.57428, -1.57833, -1.57728, -1.60381, -1.65212, -1.88782, -1.89525, -1.87041, -1.86079, -1.85612, -1.90422, -1.89704, -1.90132, -1.76801 ]
    b5 = [  0.2812, 0.33534, 0.30253, 0.2549, 0.2179, 0.20089, 0.18139, 0.17485, 0.17137, 0.177, 0.16877, 0.15617, 0.13621, 0.12916, 0.13495, 0.1392, 0.13273, 0.13096, 0.13059, 0.13535, 0.13599, 0.12485, 0.12227, 0.11447, 0.10265, 0.09129, 0.08005, 0.08173, 0.08213, 0.07577, 0.06379, 0.05697, 0.0587, 0.0686, 0.07672, 0.08428, 0.08325, 0.07663, 0.07065, 0.06525, 0.06602, 0.06774, 0.0694, 0.07429, 0.08428, 0.09052, 0.09284, 0.09077, 0.09288, 0.09428, 0.09887, 0.1068, 0.14049, 0.14356, 0.14283, 0.1434, 0.14444, 0.15127, 0.15039, 0.15081, 0.13314 ]
    b6 = [  7.86638, 7.74734, 8.21405, 8.31786, 8.21914, 7.20688, 6.54299, 6.24751, 6.57173, 6.78082, 7.17423, 6.7617, 6.10103, 5.19135, 4.46323, 4.27945, 4.37011, 4.62192, 4.65393, 4.8454, 4.97596, 5.04489, 5.00975, 5.08902, 5.03274, 5.08347, 5.14423, 5.29006, 5.3349, 5.19412, 5.1575, 5.27441, 5.54539, 5.93828, 6.36599, 6.82292, 7.11603, 7.31928, 7.25988, 7.25344, 7.26059, 7.4032, 7.46168, 7.51273, 7.77062, 7.87702, 7.91753, 7.61956, 7.59643, 7.50338, 7.53947, 7.61893, 8.12248, 7.92236, 7.49999, 7.26668, 7.11861, 7.36277, 7.45038, 7.60234, 7.2195 ]
    b7 = [  0.08753, 0.04707, 0.02667, 0.02578, 0.06557, 0.0981, 0.12847, 0.16213, 0.21222, 0.24121, 0.25944, 0.26498, 0.27718, 0.28574, 0.30348, 0.31516, 0.32153, 0.3352, 0.34849, 0.35919, 0.36619, 0.37278, 0.37756, 0.38149, 0.3812, 0.38782, 0.38862, 0.38677, 0.38625, 0.38285, 0.37867, 0.37267, 0.36952, 0.36531, 0.35936, 0.35284, 0.34775, 0.34561, 0.34142, 0.3372, 0.33298, 0.3301, 0.32645, 0.32439, 0.31354, 0.30997, 0.30826, 0.32071, 0.31801, 0.31401, 0.31104, 0.30875, 0.31122, 0.30935, 0.30688, 0.30635, 0.30534, 0.30508, 0.30362, 0.29987, 0.29772 ]
    b8 = [  0.01527, -0.02426, -0.00062, 0.01703, 0.02105, 0.03919, 0.0434, 0.06695, 0.09201, 0.11675, 0.13562, 0.14446, 0.15156, 0.15239, 0.15652, 0.16333, 0.17366, 0.1848, 0.19061, 0.19411, 0.19519, 0.19461, 0.19423, 0.19402, 0.19309, 0.19392, 0.19273, 0.19082, 0.19285, 0.19161, 0.18812, 0.18568, 0.18149, 0.17617, 0.17301, 0.16945, 0.16743, 0.1673, 0.16325, 0.16171, 0.15839, 0.15496, 0.15337, 0.15264, 0.1443, 0.1443, 0.14412, 0.14321, 0.14301, 0.14324, 0.14332, 0.14343, 0.14255, 0.14223, 0.14074, 0.14052, 0.13923, 0.13933, 0.13776, 0.13584, 0.13198 ]
    b9 = [  -0.04189, -0.0426, -0.04906, -0.04184, -0.02098, -0.04853, -0.05554, -0.04722, -0.05145, -0.05202, -0.04283, -0.04259, -0.03853, -0.03423, -0.04146, -0.0405, -0.03946, -0.03786, -0.02884, -0.02209, -0.02269, -0.02613, -0.02655, -0.02088, -0.01623, -0.01826, -0.01902, -0.01842, -0.01607, -0.01288, -0.01208, -0.00845, -0.00533, -0.00852, -0.01204, -0.01386, -0.01402, -0.01526, -0.01563, -0.01848, -0.02258, -0.02626, -0.0292, -0.03484, -0.03985, -0.04155, -0.04238, -0.04963, -0.0491, -0.04812, -0.0471, -0.04607, -0.05106, -0.05024, -0.04887, -0.04743, -0.04731, -0.04522, -0.04203, -0.03863, -0.03855 ]
    b10 = [  0.08015, 0.08649, 0.0791, 0.0784, 0.08438, 0.08577, 0.09221, 0.09003, 0.09903, 0.09943, 0.08579, 0.06945, 0.05932, 0.05111, 0.04661, 0.04253, 0.03373, 0.02867, 0.02475, 0.02502, 0.02121, 0.01115, 0.0014, 0.00148, 0.00413, 0.00413, -0.00369, -0.00897, -0.00876, -0.00564, -0.00215, -0.00047, -0.00006, -0.00301, -0.00744, -0.01387, -0.01492, -0.01192, -0.00703, -0.00351, -0.00486, -0.00731, -0.00871, -0.01225, -0.01927, -0.02322, -0.02626, -0.02342, -0.0257, -0.02643, -0.02769, -0.02819, -0.02966, -0.0293, -0.02963, -0.02919, -0.02751, -0.02776, -0.02615, -0.02487, -0.02469 ]
    σ1 = [  0.261, 0.272, 0.2728, 0.2788, 0.2821, 0.2871, 0.2902, 0.2983, 0.2998, 0.3037, 0.3078, 0.307, 0.3007, 0.3004, 0.2978, 0.2973, 0.2927, 0.2917, 0.2915, 0.2912, 0.2895, 0.2888, 0.2896, 0.2871, 0.2878, 0.2863, 0.2869, 0.2885, 0.2875, 0.2857, 0.2839, 0.2845, 0.2844, 0.2841, 0.284, 0.284, 0.2834, 0.2828, 0.2826, 0.2832, 0.2835, 0.2836, 0.2832, 0.283, 0.283, 0.283, 0.2829, 0.2815, 0.2826, 0.2825, 0.2818, 0.2818, 0.2838, 0.2845, 0.2854, 0.2862, 0.2867, 0.2869, 0.2874, 0.2872, 0.2876 ]
    σ2 = [  0.0994, 0.1142, 0.1167, 0.1192, 0.1081, 0.099, 0.0976, 0.1054, 0.1101, 0.1123, 0.1163, 0.1274, 0.143, 0.1546, 0.1626, 0.1602, 0.1584, 0.1543, 0.1521, 0.1484, 0.1483, 0.1465, 0.1427, 0.1435, 0.1439, 0.1453, 0.1427, 0.1428, 0.1458, 0.1477, 0.1468, 0.145, 0.1457, 0.1503, 0.1537, 0.1558, 0.1582, 0.1592, 0.1611, 0.1642, 0.1657, 0.1665, 0.1663, 0.1661, 0.1627, 0.1627, 0.1633, 0.1632, 0.1645, 0.1665, 0.1681, 0.1688, 0.1741, 0.1759, 0.1772, 0.1783, 0.1794, 0.1788, 0.1784, 0.1783, 0.1785 ]

    if T < 0.0 # return PGV
        b1 = -2.12833
        b2 = 1.21448
        b3 = -0.08137
        b4 = -2.46942
        b5 = 0.22349
        b6 = 6.41443
        b7 = 0.20354
        b8 = 0.08484
        b9 = -0.05856
        b10 = 0.01305
        τlog10 = 0.1083
        ϕlog10 = 0.2562

        logPGV = b1 + b2*M + b3*M^2 + (b4 + b5*M)*log10(sqrt( Rjb^2 + b6^2 )) + b7*Ssoft + b8*Sstiff + b9*Fnm + b10*Frv
        PGV = 10.0 .^ logPGV
        τ = τlog10 * log(10)
        ϕ = ϕlog10 * log(10)
        σ = sqrt( τ^2 + ϕ^2 )

        # PGV in units of cm/s
        return PJSgroundMotion(PGV, log(PGV), τ, ϕ, σ)
    else
        tidHI = findfirst(Ti .>= T)
        tidLO = findlast(Ti .<= T)
    end

    if tidHI == tidLO
        # we have a known period
        tid = tidHI

        μlogSa = b1[tid] + b2[tid]*M + b3[tid]*M^2 + (b4[tid] + b5[tid]*M)*log10(sqrt( Rjb^2 + b6[tid]^2 )) + b7[tid]*Ssoft + b8[tid]*Sstiff + b9[tid]*Fnm + b10[tid]*Frv - log10(980.665)
        Sa = 10.0 .^ μlogSa
		lnSa = log(Sa)
        τ = σ2[tid] * log(10)
        ϕ = σ1[tid] * log(10)
        σ = sqrt( τ^2 + ϕ^2 )
    else
        # interpolate
        @inbounds TiLO = Ti[tidLO]
        @inbounds TiHI = Ti[tidHI]

        Tfrac = log(T/TiLO)/log(TiHI/TiLO)

        abLO = PJSab2010(TiLO, M, Rjb, Sstiff, Ssoft, Fnm, Frv)
        abHI = PJSab2010(TiHI, M, Rjb, Sstiff, Ssoft, Fnm, Frv)

		# apply conversion to 'g' units from cm/s2
        lnSa = abLO.lnIM + Tfrac*(abHI.lnIM - abLO.lnIM) - log(980.665)
        Sa = exp(lnSa)

        τ = abLO.τ + Tfrac*(abHI.τ - abLO.τ)
		ϕ = abLO.ϕ + Tfrac*(abHI.ϕ - abLO.ϕ)
		σ = abLO.σ + Tfrac*(abHI.σ - abLO.σ)
	end

	gm = PJSgroundMotion(Sa, lnSa, τ, ϕ, σ)
	return gm
end



function PJSab2010(T::U, M::U, Rjb::U, Vs30::U, Fnm::Int64, Frv::Int64) where U<:Real
	Sstiff = 0
	Ssoft = 0
	if Vs30 < 360.0
		Ssoft = 1
	elseif Vs30 < 750.0
		Sstiff = 1
	end
	return PJSab2010(T, M, Rjb, Sstiff, Ssoft, Fnm, Frv)
end


function PJSab2010(Ti::Vector{U}, M::U, Rjb::U, Ssoft::Int64, Sstiff::Int64, Fnm::Int64, Frv::Int64) where U<:Real
	predictions = Vector{PJSgroundMotion}()
	for T in Ti
		push!(predictions, PJSab2010(T, M, Rjb, Ssoft, Sstiff, Fnm, Frv))
	end
	return predictions
end

function PJSab2010(Ti::Vector{U}, M::U, Rjb::U, Vs30::U, Fnm::Int64, Frv::Int64) where U<:Real
	predictions = Vector{PJSgroundMotion}()
	for T in Ti
		push!(predictions, PJSab2010(T, M, Rjb, Vs30, Fnm, Frv))
	end
	return predictions
end


"""
	PJSab2010(T::U, rup::Rupture{U}, site::Site{U}) where U<:Real

Akkar & Bommer (2010) model for PGA, PGV and response spectral ordinates (5% damped)
"""
function PJSab2010(T::U, rup::Rupture{U}, site::Site{U}) where U<:Real
	Rjb = joyner_boore_distance(rup, site)
	Fss, Fnm, Frv, Fuk = mechanism(rup)
	gm = PJSab2010(T, rup.magnitude, Rjb, vs30(site), Fnm, Frv)
	return gm
end

function PJSab2010(Ti::Vector{U},  rup::Rupture{U}, site::Site{U}) where U<:Real
	predictions = Vector{PJSgroundMotion}()
	for T in Ti
		push!(predictions, PJSab2010(T, rup, site))
	end
	return predictions
end
