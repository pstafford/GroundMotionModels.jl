
import Base.+
import Base.-

"""

Adopt a coordinate system where x=East, y=North, z=Depth
Note that this isn't a right triad, but makes plotting more natural
"""
struct Point{T<:Real}
    x::T
    y::T
    z::T
end

Point() = Point(0.0, 0.0, 0.0)
Point(x, y) = Point(x, y, 0.0)


"""
    -(pl::Point, pr::Point)

Difference between two points (this defines a Vector3D whose x and y coordinates are those of the point)
"""
function -(pl::Point, pr::Point)
    # return Point(pl.x-pr.x, pl.y-pr.y, pl.z-pr.z)
    return Vector3D(pl.x-pr.x, pl.y-pr.y, pl.z-pr.z)
end

"""
    +(pl::Point, pr::Point)

Sum of two points (this defines a Vector3D whose x and y coordinates are those of the point)
"""
function +(pl::Point, pr::Point)
    # return Point(pl.x+pr.x, pl.y+pr.y, pl.z+pr.z)
    return Vector3D(pl.x+pr.x, pl.y+pr.y, pl.z+pr.z)
end


"""
    distance(p1::Point, p2::Point)

Returns the distance between two points
- `p1` is the first Point
- `p2` is the second Point
"""
function distance(p1::Point, p2::Point)
    return sqrt( (p1.x-p2.x)^2 + (p1.y-p2.y)^2 + (p1.z-p2.z)^2 )
end

"""
    squared_distance(p1::Point, p2::Point)

Returns the squared distance between two points
- `p1` is the first Point
- `p2` is the second Point
"""
function squared_distance(p1::Point, p2::Point)
    return (p1.x-p2.x)^2 + (p1.y-p2.y)^2 + (p1.z-p2.z)^2
end

"""
    surface_distance(p1::Point, p2::Point)

Returns the surface distance between two points
- `p1` is the first Point
- `p2` is the second Point
"""
function surface_distance(p1::Point, p2::Point)
    return sqrt( (p1.x-p2.x)^2 + (p1.y-p2.y)^2 )
end


"""
    point_at_cartesian_offset(p::Point{T}, Δx::T, Δy::T) where T<:Real

Finds a new point that is some offset in the x and y directions
- `Δx` is the offset in the x-direction
- `Δy` is the offset in the y-direction
"""
function point_at_cartesian_offset(p::Point{T}, Δx::T, Δy::T) where T<:Real
    return Point(p.x+Δx, p.y+Δy, p.z)
end


"""
    point_at_distance_and_azimuth(p::Point{T}, d::T, θ::T) where T<:Real

Finds a new point that is some distance and bearing from `p`
- `p` is the starting point
- `d` is the distance
- `θ` is the azimuth defined as a bearing in degrees (from North==y)
"""
function point_at_distance_and_azimuth(p::Point{T}, d::T, θ::T) where T<:Real
    Δx = d * sind(θ)
    Δy = d * cosd(θ)
    return point_at_cartesian_offset(p, Δx, Δy)
end


"""
    surface_point(p::Point)

Return the surface-projection of a Point
"""
function surface_point(p::Point)
    return Point(p.x, p.y, 0.0)
end


"""
    dot(p1::Point, p2::Point)

Dot product between the coordinates of two points
"""
function dot(p1::Point, p2::Point)
    return p1.x*p2.x + p1.y*p2.y + p1.z*p2.z
end

"""
    cross(p1::Point, p2::Point)

Cross product between two Points (using the coordinates of the points as components of a vector)
"""
function cross(p1::Point, p2::Point)
    x = p1.y*p2.z - p1.z*p2.y
    y = p1.z*p2.x - p1.x*p2.z
    z = p1.x*p2.y - p1.y*p2.x
    return Point(x, y, z)
end
