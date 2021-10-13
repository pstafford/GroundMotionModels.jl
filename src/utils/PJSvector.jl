

"""
Vector2D{T<:Real}

Basic 2D vector
"""
struct Vector2D{T<:Real}
    x::T
    y::T
end

Vector2D() = Vector2D(0.0, 0.0)

"""
Vector3D{T<:Real}

Basic 3D vector
"""
struct Vector3D{T<:Real}
    x::T
    y::T
    z::T
end

Vector3D() = Vector3D(0.0, 0.0, 0.0)

"""
    magnitude(v::Vector2D)

Compute the magnitude of a 2D vector
"""
function magnitude(v::Vector2D)
    return sqrt( v.x^2 + v.y^2 )
end

"""
    magnitude(v::Vector3D)

Compute the magnitude of a 3D vector
"""
function magnitude(v::Vector3D)
    return sqrt( v.x^2 + v.y^2 + v.z^2 )
end

"""
    dot(v1::Vector2D, v2::Vector2D)

Dot product between two 2D vectors
"""
function dot(v1::Vector2D{T}, v2::Vector2D{T}) where T<:Real
    return v1.x*v2.x + v1.y*v2.y
end

"""
    dot(v1::Vector3D, v2::Vector3D)

Dot product between two 3D vectors
"""
function dot(v1::Vector3D{T}, v2::Vector3D{T}) where T<:Real
    return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z
end

"""
    cross(v1::Vector3D, v2::Vector3D)

Cross product between two 3D vectors
"""
function cross(v1::Vector3D{T}, v2::Vector3D{T}) where T<:Real
    x = v1.y*v2.z - v1.z*v2.y
    y = v1.z*v1.x - v1.x*v2.z
    z = v1.x*v2.y - v1.y*v2.x
    return Vector3D(x, y, z)
end
