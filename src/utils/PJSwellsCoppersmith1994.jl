


"""
wells_coppersmith_1994(m; sof, output, epsilon)

computes source dimensions as a function of magnitude
	- output is
"""
function wells_coppersmith_1994( m ; sof::String="SS", output::String="SRL", epsilon::Float64=0.0 )
	if sof == "SS"
		# strike-slip rupture
		if output == "SRL"
			# surface rupture
			srlMed = 10.0 .^ ( -3.55 + 0.74 * m )
			sdSRL = 0.23*log(10)
			return srlMed .* exp(epsilon * sdSRL)
		elseif output == "RLD"
			# rupture length at depth
			rlMed = 10.0 .^ ( -2.57 + 0.62 * m )
			sdRL = 0.15*log(10)
			return rlMed .* exp(epsilon * sdRL)
		elseif output == "RW"
			# rupture width
			rwMed = 10.0 .^ ( -0.76 + 0.27 * m )
			sdRW = 0.14*log(10)
			return rwMed .* exp(epsilon * sdRW)
		elseif output == "RA"
			# rupture area
			raMed = 10.0 .^ ( -3.42 + 0.90 * m )
			sdRA = 0.22*log(10)
			return raMed .* exp(epsilon * sdRA)
		end
	elseif sof == "NM"
		# normal rupture
		if output == "SRL"
			srlMed = 10.0 .^ ( -2.01 + 0.50 * m )
			sdSRL = 0.21*log(10)
			return srlMed .* exp(epsilon * sdSRL)
		elseif output == "RLD"
			rlMed = 10.0 .^ ( -1.88 + 0.50 * m )
			sdRL = 0.17*log(10)
			return rlMed .* exp(epsilon * sdRL)
		elseif output == "RW"
			rwMed = 10.0 .^ ( -1.14 + 0.35 * m )
			sdRW = 0.12*log(10)
			return rwMed .* exp(epsilon * sdRW)
		elseif output == "RA"
			raMed = 10.0 .^ ( -2.87 + 0.82 * m )
			sdRA = 0.22*log(10)
			return raMed .* exp(epsilon * sdRA)
		end
	elseif sof == "RV"
		# reverse rupture
		if output == "SRL"
			srlMed = 10.0 .^ ( -2.86 + 0.64 * m )
			sdSRL = 0.20*log(10)
			return srlMed .* exp(epsilon * sdSRL)
		elseif output == "RLD"
			rlMed = 10.0 .^ ( -2.42 + 0.58 * m )
			sdRL = 0.16*log(10)
			return rlMed .* exp(epsilon * sdRL)
		elseif output == "RW"
			rwMed = 10.0 .^ ( -1.61 + 0.41 * m )
			sdRW = 0.15*log(10)
			return rwMed .* exp(epsilon * sdRW)
		elseif output == "RA"
			raMed = 10.0 .^ ( -3.99 + 0.98 * m )
			sdRA = 0.26*log(10)
			return raMed .* exp(epsilon * sdRA)
		end
	elseif sof == "ALL"
		# all ruptures
		if output == "SRL"
			srlMed = 10.0 .^ ( -3.22 + 0.69 * m )
			sdSRL = 0.22*log(10)
			return srlMed .* exp(epsilon * sdSRL)
		elseif output == "RLD"
			rlMed = 10.0 .^ ( -2.44 + 0.59 * m )
			sdRL = 0.16*log(10)
			return rlMed .* exp(epsilon * sdRL)
		elseif output == "RW"
			rwMed = 10.0 .^ ( -1.01 + 0.32 * m )
			sdRW = 0.15*log(10)
			return rwMed .* exp(epsilon * sdRW)
		elseif output == "RA"
			raMed = 10.0 .^ ( -3.49 + 0.91 * m )
			sdRA = 0.24*log(10)
			return raMed .* exp(epsilon * sdRA)
		end
	end
end


function rupture_width(m::T, rake::T, eps::T=0.0) where T<:Real
	output = "RW"
	if isnan(rake) # unknown mechanism
		sof = "ALL"
	else
		if (rake > 30.0 && rake < 150.0)
			# reverse or reverse-obliquq
			sof = "RV"
		elseif (rake < -30.0 && rake > -150.0)
			# normal or normal oblique
			sof = "NM"
		else
			# strike-slip
			sof = "SS"
		end
	end

	return wells_coppersmith_1994(m; sof=sof, output=output, epsilon=eps)
end


function rupture_length(m::T, rake::T, eps::T=0.0) where T<:Real
	output = "RLD"
	if isnan(rake) # unknown mechanism
		sof = "ALL"
	else
		if (rake > 30.0 && rake < 150.0)
			# reverse or reverse-obliquq
			sof = "RV"
		elseif (rake < -30.0 && rake > -150.0)
			# normal or normal oblique
			sof = "NM"
		else
			# strike-slip
			sof = "SS"
		end
	end

	return wells_coppersmith_1994(m; sof=sof, output=output, epsilon=eps)
end
