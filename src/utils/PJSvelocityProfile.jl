

"""
    VelocityProfile{T<:Real}

Custom type representing a velocity profile. Holds vectors of layer thicknesses and velocities.
Properties are `thicknesses::Vector{T}` and `velocities::Vector{T}` where `T<:Real`.
"""
struct VelocityProfile{T<:Real}
    thicknesses::Vector{T}
    velocities::Vector{T}
end

"""
    vs30(profile::VelocityProfile{T}) where T<:Real

Compute time-averaged shear-wave velocity [m/s] over uppermost 30 m of a `VelocityProfile` instance.
"""
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


function z1p0(profile::VelocityProfile{T}) where T<:Real
    return depth_to_velocity_horizon(profile, 1000.0)
end



function z2p5(profile::VelocityProfile{T}) where T<:Real
    return depth_to_velocity_horizon(profile, 2500.0)
end
