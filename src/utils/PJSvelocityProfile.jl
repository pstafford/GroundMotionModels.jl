

struct VelocityProfile{T<:Real}
    thicknesses::Vector{T}
    velocities::Vector{T}
end


function vs30(profile::VelocityProfile{T}) where T<:Real
    nlayers = length(profile.thicknesses)
    if nlayers == 1
        if profile.thicknesses[1] >= 30.0
            return profile.velocities[1]
        else
            return NaN
        end
    else
        cti = cumsum(profile.thicknesses)
        zid = findall(cti .< 30.0)
        nzlo = length(zid)
        cz = 0.0
        tt = 0.0
        for i in 1:nzlo
            cz += profile.thicknesses[i]
            tt += profile.thicknesses[i] / profile.velocities[i]
        end
        # add last layer
        dz = 30.0 - cz
        tt += dz / profile.velocities[nzlo+1]
        return 30.0 / tt
    end
end

zi = [ 15.0, 5.0, 5.0, 10.0 ]
vi = [ 200.0, 400.0, 600.0, 1000.0 ]

profile = VelocityProfile(zi, vi)

@time vs30(profile)

# 30.0 / (15.0/200.0 + 5.0/400.0 + 5.0/600.0 + 5.0/1000.0)

function depth_to_velocity_horizon(profile::VelocityProfile{T}, target::T) where T<:Real
    lid = findfirst(profile.velocities .>= target)
    if lid == nothing
        print("PJS Error: profile never reaches target velocity of $target m/s\n")
        return NaN
    else
        if lid == 1
            return 0.0
        else
            ct = cumsum(profile.thicknesses)
            return ct[lid-1]
        end
    end
end

@time depth_to_velocity_horizon(profile, 1500.0)

function z1p0(profile::VelocityProfile{T}) where T<:Real
    return depth_to_velocity_horizon(profile, 1000.0)
end

@time z1p0(profile)


function z2p5(profile::VelocityProfile{T}) where T<:Real
    return depth_to_velocity_horizon(profile, 2500.0)
end

vp = VelocityProfile([5.0, 10.0, 20.0],[500.0, 700.0, 1500.0])
z1p0(vp)
z2p5(vp)
