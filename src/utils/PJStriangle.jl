

struct Triangle{T<:Real}
    p1::Point{T}
    p2::Point{T}
    p3::Point{T}
end

Triangle() = Triangle(Point(), Point(), Point())

p1 = Point()
p2 = Point(1.0, 0.0, 0.0)
p3 = Point(0.0, 1.0, 0.0)

tri = Triangle(p1, p2, p3)


# """
#
# Determine if the triangle is oriented in a clockwise direction
# """
# function clockwise(tri::Triangle)
#     v1 = tri.p2 - tri.p1
#     v2 = tri.p3 - tri.p1
#
#     vc = cross(v1, v2)
#     dc = dot(vc, vc)
#     return ( dc >= 0 ) ? true : false
# end
#
# clockwise(tri)
#
# tri2 = Triangle(p1, p3, p2)
# clockwise(tri2)


"""
    distance_to_point(t::Triangle{T}, p::Point{T}) where T<:Real

Function to compute the distance from the triangle to a given point

This method computes the distance from the triangle to a point in 3D
This algorithm uses no trigonometric functions and only one square root and is thus very efficient (at the expense of being verbose)

The algorithm is based on
"David Eberly, 'Distance Between Point and Triangle in 3D', Geometric Tools, LLC, (1999)"
www.geometrictools.com/Documentation/DistancePoint3Triangle3.pdf

- parameter point is a point in 3D space
- returns a double value defining the distance to the triangle
"""
function distance_to_point(tri::Triangle{T}, p::Point{T}) where T<:Real
    # define edge vectors p2-p1
    E1 = tri.p2 - tri.p1

    # define edge p3-p1
    E2 = tri.p3 - tri.p1

    # define the vector connecting p0 and the site: p0-site
    D = tri.p1 - p

    # compute a series of dot products
    a = dot(E1, E1)
    b = dot(E1, E2)
    c = dot(E2, E2)
    d = dot(E1, D)
    e = dot(E2, D)
    f = dot(D, D)

    det = a*c - b*b
    s = b*e - c*d
    t = b*d - a*e

    if (s+t) <= det
        if s < 0.0
            if t < 0.0
                #region 4
                if d < 0.0
                    if -d >= a
                        sqrDistance = a + 2.0*d + f
                    else
                        s = -d/a
                        sqrDistance = d*s + f
                    end
                else
                    if e >= 0.0
                        sqrDistance = f
                    else
                        if -e >= c
                            sqrDistance = c + 2.0*e + f
                        else
                            t = -e/c
                            sqrDistance = e*t + f
                        end
                    end
                end #of region 4
            else
                # region 3
                if e >= 0.0
                    sqrDistance = f
                else
                    if -e >= c
                        sqrDistance = c + 2.0*e + f
                    else
                        t = -e/c
                        sqrDistance = e*t + f
                    end
                end
            end #of region 3
        else
            if t < 0.0
                # region 5
                if d >= 0.0
                    sqrDistance = f
                else
                    if -d >= a
                        sqrDistance = a + 2.0*d + f
                    else
                        s = -d/a
                        sqrDistance = d*s + f
                    end
                end
            else
                # region 0
                invDet = 1.0/det
                s = s*invDet
                t = t*invDet
                sqrDistance = s*(a*s + b*t + 2.0*d) + t*(b*s + c*t + 2.0*e) + f
            end
        end
    else
        if s < 0.0
            # region 2
            tmp0 = b + d
            tmp1 = c + e
            if tmp1 > tmp0
                # minimum on edge s+t=1
                numer = tmp1 - tmp0
                denom = a - 2.0*b + c
                if numer >= denom
                    sqrDistance = a + 2.0*d + f
                else
                    s = numer/denom
                    t = 1.0 - s
                    sqrDistance = s*(a*s + b*t + 2.0*d) + t*(b*s + c*t + 2.0*e) + f
                end
            else
                # minimum on edge s=0
                if tmp1 <= 0.0
                    sqrDistance = c + 2.0*e + f
                else
                    if e >= 0.0
                        sqrDistance = f
                    else
                        t = -e/c
                        sqrDistance = e*t + f
                    end
                end
            end #of region 2
        else
            if t < 0.0
                #region 6
                tmp0 = b + e
                tmp1 = a + d
                if tmp1 > tmp0
                    numer = tmp1 - tmp0
                    denom = a - 2.0*b + c
                    if numer >= denom
                        sqrDistance = c + 2.0*e + f
                    else
                        t = numer/denom
                        s = 1.0 - t
                        sqrDistance = s*(a*s + b*t + 2.0*d) + t*(b*s + c*t + 2.0*e) + f
                    end
                else
                    if tmp1 <= 0
                        sqrDistance = a + 2.0*d + f
                    else
                        if d >= 0.0
                            sqrDistance = f
                        else
                            s = -d/a
                            sqrDistance = d*s + f
                        end
                    end
                end #end region 6
            else
                # region 1
                numer = c + e - b - d
                if numer <= 0.0
                    sqrDistance = c + 2.0*e + f
                else
                    denom = a - 2.0*b + c
                    if numer >= denom
                        sqrDistance = a + 2.0*d + f
                    else
                        s = numer/denom
                        t = 1.0 - s
                        sqrDistance = s*(a*s + b*t + 2.0*d) + t*(b*s + c*t + 2.0*e) + f
                    end
                end #of region 1
            end
        end
    end

    # account for numerical round-off error
    sqrDistance = ( sqrDistance < 0.0 ) ? 0.0 : sqrDistance

    return sqrt(sqrDistance)
end

# p0 = Point(0.1, 0.1, 0.0)
# @time distance_to_point(tri, p0)
#
# tt(tri, p0) = @time distance_to_point(tri, p0)
# tt(tri, p0)
