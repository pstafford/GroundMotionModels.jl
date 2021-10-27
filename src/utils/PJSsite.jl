
"""
    Site{T<:Real}

Custom type representing a site location.

Holds properties of:
- `position::Point{T}` representing the site coordinates (Cartesian)
- `profile::VelocityProfile{T}` a velocity profile
- `measured::Int64` binary indicator representing whether the velocities are measured or estimated
- `region::String` string defining which generic region the site is in (useful for regionally varying GMMs)
"""
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
"""
    vs30(s::Site)

Computes the time-averaged shear-wave velocity over the uppermost 30 m for the site
"""
vs30(s::Site) = vs30(s.profile)
z1p0(s::Site) = z1p0(s.profile)
z2p5(s::Site) = z2p5(s.profile)
depth_to_velocity_horizon(s::Site{T}, v::T) where T<:Real = depth_to_velocity_horizon(s.profile, v)
