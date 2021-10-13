
include("PJSpoint.jl")
include("PJSvelocityProfile.jl")


struct Site{T<:Real}
    position::Point{T}
    profile::VelocityProfile{T}
    measured::Int64
    region::String
end

Site() = Site(Point(), VelocityProfile([30.0],[760.0]), 1, "California")
Site(p::Point{T}, prof::VelocityProfile{T}) where T<:Real = Site(p, prof, 1, "California")
Site(p::Point{T}, prof::VelocityProfile{T}, meas::Int64) where T<:Real = Site(p, prof, meas, "California")

# utility functions
vs30(s::Site) = vs30(s.profile)
z1p0(s::Site) = z1p0(s.profile)
z2p5(s::Site) = z2p5(s.profile)
depth_to_velocity_horizon(s::Site{T}, v::T) where T<:Real = depth_to_velocity_horizon(s.profile, v)

# s = Site(Point(), VelocityProfile([20.0, 20.0],[500.0, 1000.0]), true)
# vs30(s)
# z1p0(s)
