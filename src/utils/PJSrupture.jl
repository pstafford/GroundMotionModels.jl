
include("PJSpoint.jl")
include("PJStriangle.jl")
include("PJSlineSegment.jl")
include("PJSwellsCoppersmith1994.jl")
include("PJSsite.jl")


struct Rupture{T<:Real}
    magnitude::T
    points::Vector{Point{T}}
    rake::T
    rupLenth::T
    rupWidth::T
    hypocentre::Point{T}
end


rup = Rupture(6.0, [Point()], 0.0, NaN, NaN, Point(NaN,NaN,NaN))

function add_point!(rup::Rupture{T}, p::Point{T}) where T<:Real
    push!(rup.points, p)
end

function remove_last_point!(rup::Rupture)
    pop!(rup.points)
end

function remove_point_at_index!(rup::Rupture, index::Int64)
    if index > 0 && index <= length(rup.points)
        deleteat!(rup.points, index)
    end
end


add_point!(rup, Point(1.0, 0.0))
rup


function ztor(rup::Rupture{T}) where T<:Real
    return rup.points[1].z
end

function zbot(rup::Rupture{T}) where T<:Real
    return rup.points[4].z
end

function zhyp(rup::Rupture{T}) where T<:Real
    if isnan(rup.hypocentre.z)
        ztor = ztor(rup)
        zbot = zbot(rup)
        return ztor + 2.0/3.0*(zbot - ztor)
    else
        return rup.hypocentre.z
    end
end

"""
mechanism(rup::Rupture{T}) where T<:Real

Returns the style-of-faulting indicator variables Fss, Fnm, Frv, Fuk based upon the
rake angle for the rupture
"""
function mechanism(rup::Rupture{T}) where T<:Real
    if isnan(rup.rake) # unknown mechanism
        return 0, 0, 0, 1
    else
        if (rup.rake > 30.0 && rup.rake < 150.0)
            # reverse or reverse-obliquq
            return 0, 0, 1, 0
        elseif (rup.rake < -30.0 && rup.rake > -150.0)
            # normal or normal oblique
            return 0, 1, 0, 0
        else
            # strike-slip
            return 1, 0, 0, 0
        end
    end
end

mechanism(rup)

"""
    strike(rup::Rupture{T}) where T<:Real

Computes the strike (in degrees) for a given rupture. This is based upon the orientation of the top edge of the fault
"""
function strike(rup::Rupture{T}) where T<:Real
    # work with p1 and p2 of the rupture
    Δp = rup.points[2] - rup.points[1]
    Δx = Δp.x
    Δy = Δp.y

    # reverse the normal ordering of inputs to this function
    θ = atand(Δx, Δy)
    return θ
end

# @time strike(rup)

"""
    dip(rup::Rupture{T}) where T<:Real

Computes the dip (in degrees) for a given rectangular rupture.
"""
function dip(rup::Rupture{T}) where T<:Real
    p2 = rup.points[2]
    p3 = rup.points[3]
    v1 = p3 - p2
    v2 = Point(p3.x, p3.y, p2.z) - p2
    mv2 = magnitude(v2)
    if mv2 == 0.0
        δ = 90.0
    else
        mv1 = magnitude(v1)
        δ = acosd( dot(v1, v2)/ (mv1*mv2) )
    end
    return δ
end

# @time dip(rup)

function hanging_wall(rup::Rupture{T}) where T<:Real
    upper_segment = LineSegment(rup.points[1], rup.points[2])
    if left_of_line_segment(upper_segment, p)
        return 0
    else
        return 1
    end
end


function rupture_from_hypocentre(hyp::Point{T}, θ::T, δ::T, L::T, W::T, xL::T, xW::T) where T<:Real
    # obtain the increments to each point
    ΔLf = (1.0 - xL) * L
    ΔLb = xL * L
    ΔWu = (1.0 - xW) * W
    ΔWd = xW * W

    Δxfu = ΔLf * sind(θ) - ΔWu * cosd(δ) * cosd(θ)
    Δxfd = ΔLf * sind(θ) + ΔWd * cosd(δ) * cosd(θ)
    Δxbu = -ΔLb * sind(θ) - ΔWu * cosd(δ) * cosd(θ)
    Δxbd = -ΔLb * sind(θ) + ΔWd * cosd(δ) * cosd(θ)

    Δyfu = ΔLf * cosd(θ) + ΔWu * cosd(δ) * sind(θ)
    Δyfd = ΔLf * cosd(θ) - ΔWd * cosd(δ) * sind(θ)
    Δybu = -ΔLb * cosd(θ) + ΔWu * cosd(δ) * sind(θ)
    Δybd = -ΔLb * cosd(θ) - ΔWd * cosd(δ) * sind(θ)

    Δzfu = -ΔWu * sind(δ)
    Δzbu = -ΔWu * sind(δ)
    Δzfd = ΔWd * sind(δ)
    Δzbd = ΔWd * sind(δ)

    pbu = Point(hyp.x+Δxbu, hyp.y+Δybu, hyp.z+Δzbu)
    pfu = Point(hyp.x+Δxfu, hyp.y+Δyfu, hyp.z+Δzfu)
    pfd = Point(hyp.x+Δxfd, hyp.y+Δyfd, hyp.z+Δzfd)
    pbd = Point(hyp.x+Δxbd, hyp.y+Δybd, hyp.z+Δzbd)

    return Rupture(NaN, [pbu, pfu, pfd, pbd], NaN, L, W, hyp)
end

hyp = Point(0.0, 0.0, 10.0)
rup = rupture_from_hypocentre(hyp, 90.0, 90.0, 10.0, 10.0, 0.5, 0.5)


function rupture_from_hypocentre_with_source_scaling(hyp::Point{T}, θ::T, δ::T, M::T, λ::T, xL::T, xW::T, scaling_relation::String="WellsCoppersmith1994") where T<:Real
    # generate source dimensions
    if scaling_relation == "WellsCoppersmith1994"
        W = rupture_width(M, λ, 0.0)
        L = rupture_length(M, λ, 0.0)
    else
        return Rupture(NaN, [Point()], 0.0, NaN, NaN, Point(NaN,NaN,NaN))
    end

    # obtain the increments to each point
    xWmax = hyp.z / sind(δ)
    if (1.0-xW) * W > xWmax
        ΔWu = xWmax
        ΔWd = W - xWmax
    else
        ΔWu = (1.0 - xW) * W
        ΔWd = xW * W
    end
    ΔLf = (1.0 - xL) * L
    ΔLb = xL * L

    Δxfu = ΔLf * sind(θ) - ΔWu * cosd(δ) * cosd(θ)
    Δxfd = ΔLf * sind(θ) + ΔWd * cosd(δ) * cosd(θ)
    Δxbu = -ΔLb * sind(θ) - ΔWu * cosd(δ) * cosd(θ)
    Δxbd = -ΔLb * sind(θ) + ΔWd * cosd(δ) * cosd(θ)

    Δyfu = ΔLf * cosd(θ) + ΔWu * cosd(δ) * sind(θ)
    Δyfd = ΔLf * cosd(θ) - ΔWd * cosd(δ) * sind(θ)
    Δybu = -ΔLb * cosd(θ) + ΔWu * cosd(δ) * sind(θ)
    Δybd = -ΔLb * cosd(θ) - ΔWd * cosd(δ) * sind(θ)

    Δzfu = -ΔWu * sind(δ)
    Δzbu = -ΔWu * sind(δ)
    Δzfd = ΔWd * sind(δ)
    Δzbd = ΔWd * sind(δ)

    pbu = Point(hyp.x+Δxbu, hyp.y+Δybu, hyp.z+Δzbu)
    pfu = Point(hyp.x+Δxfu, hyp.y+Δyfu, hyp.z+Δzfu)
    pfd = Point(hyp.x+Δxfd, hyp.y+Δyfd, hyp.z+Δzfd)
    pbd = Point(hyp.x+Δxbd, hyp.y+Δybd, hyp.z+Δzbd)

    return Rupture(M, [pbu, pfu, pfd, pbd], λ, L, W, hyp)
end

rup = rupture_from_hypocentre_with_source_scaling(Point(0.0,0.0,5.0), 90.0, 90.0, 6.0, 0.0, 0.5, 0.5)


function rupture_distance(rup::Rupture{T}, p::Point{T}) where T<:Real
    # split the rupture into two triangles
    t1 = Triangle(rup.points[1], rup.points[4], rup.points[2])
    d1 = distance_to_point(t1, p)
    if d1 > 0.0
        t2 = Triangle(rup.points[4], rup.points[3], rup.points[2])
        d2 = distance_to_point(t2, p)
        return (d1 < d2) ? d1 : d2
    else
        return 0.0
    end
end

p = Point(10.0, 5.0, 0.0)
@time r_rup = rupture_distance(rup, p)

function rupture_distance(rup::Rupture{T}, points::Vector{Point{T}}) where T<:Real
    n = length(points)
    r_rups = zeros(T, n)
    for i in 1:n
        r_rups[i] = rupture_distance(rup, points[i])
    end
    return r_rups
end

rupture_distance(rup::Rupture{T}, site::Site{T}) where T<:Real = rupture_distance(rup, site.position)
rupture_distance(rup::Rupture{T}, sites::Vector{Site{T}}) where T<:Real = rupture_distance(rup, map(s -> s.position, sites))

@time rupture_distance(rup, [Point(0.0, 5.0, 0.0), Point(0.0, 10.0, 0.0)])

function joyner_boore_distance(rup::Rupture{T}, p::Point{T}) where T<:Real
    if dip(rup) == 90.0
        # work with linesegment
        seg = LineSegment(surface_point(rup.points[1]), surface_point(rup.points[2]))

        return distance_to_line_segment(seg, p)
    else
        # map all points to the surface
        surface_rup = deepcopy(rup)
        for i in 1:length(surface_rup.points)
            surface_rup.points[i] = surface_point(surface_rup.points[i])
        end
        surface_p = surface_point(deepcopy(p))

        return rupture_distance(surface_rup, surface_p)
    end
end

@time joyner_boore_distance(rup, p)

function joyner_boore_distance(rup::Rupture{T}, points::Vector{Point{T}}) where T<:Real
    n = length(points)
    r_jbs = zeros(T, n)
    for i in 1:n
        r_jbs[i] = joyner_boore_distance(rup, points[i])
    end
    return r_jbs
end

joyner_boore_distance(rup::Rupture{T}, site::Site{T}) where T<:Real = joyner_boore_distance(rup, site.position)
joyner_boore_distance(rup::Rupture{T}, sites::Vector{Site{T}}) where T<:Real = joyner_boore_distance(rup, map(s -> s.position, sites))


function epicentral_distance(rup::Rupture{T}, p::Point{T}) where T<:Real
    Δx = p.x - rup.hypocentre.x
    Δy = p.y - rup.hypocentre.y
    return sqrt( Δx^2 + Δy^2 )
end

function epicentral_distance(rup::Rupture{T}, points::Vector{Point{T}}) where T<:Real
    n = length(points)
    r_epis = zeros(T, n)
    for i in 1:n
        r_epis[i] = epicentral_distance(rup, points[i])
    end
    return r_epis
end

epicentral_distance(rup::Rupture{T}, site::Site{T}) where T<:Real = epicentral_distance(rup, site.position)
epicentral_distance(rup::Rupture{T}, sites::Vector{Site{T}}) where T<:Real = epicentral_distance(rup, map(s -> s.position, sites))


function hypocentral_distance(rup::Rupture{T}, p::Point{T}) where T<:Real
    Δx = p.x - rup.hypocentre.x
    Δy = p.y - rup.hypocentre.y
    Δz = p.z - rup.hypocentre.z
    return sqrt( Δx^2 + Δy^2 + Δz^2 )
end

function hypocentral_distance(rup::Rupture{T}, points::Vector{Point{T}}) where T<:Real
    n = length(points)
    r_hyps = zeros(T, n)
    for i in 1:n
        r_hyps[i] = hypocentral_distance(rup, points[i])
    end
    return r_hyps
end

hypocentral_distance(rup::Rupture{T}, site::Site{T}) where T<:Real = hypocentral_distance(rup, site.position)
hypocentral_distance(rup::Rupture{T}, sites::Vector{Site{T}}) where T<:Real = hypocentral_distance(rup, map(s -> s.position, sites))


function strike_distance(rup::Rupture{T}, p::Point{T}) where T<:Real
    upper_segment = LineSegment(rup.points[1], rup.points[2])
    d = distance_to_line_containing_segment(upper_segment, p)
    if left_of_line_segment(upper_segment, p)
        return -d
    else
        return d
    end
end

@time strike_distance(rup, p)

function strike_distance(rup::Rupture{T}, points::Vector{Site{T}}) where T<:Real
    n = length(points)
    r_xs = zeros(T, n)
    for i in 1:n
        r_xs[i] = strike_distance(rup, points[i])
    end
    return r_xs
end

strike_distance(rup::Rupture{T}, site::Site{T}) where T<:Real = strike_distance(rup, site.position)
strike_distance(rup::Rupture{T}, sites::Vector{Site{T}}) where T<:Real = strike_distance(rup, map(s -> s.position, sites))

function strike_parallel_distance(rup::Rupture{T}, p::Point{T}) where T<:Real
    if dip(rup) == 90.0
        # make a small adjustment to effectively introduce a slight dip to enable the segments to be identified properly
        Δ = 0.1
        p3 = surface_point(rup.points[3])
        p4 = surface_point(rup.points[4])

        # get the strike for the rupture
        θ = strike(rup)
        # apply the increment
        sp3 = Point(p3.x + Δ*cosd(θ), p3.y - Δ*sind(θ), 0.0)
        sp4 = Point(p4.x + Δ*cosd(θ), p4.y - Δ*sind(θ), 0.0)

        forward_seg = LineSegment(surface_point(rup.points[2]), sp3)
        backward_seg = LineSegment(sp4, surface_point(rup.points[1]))
    else
        forward_seg = LineSegment(surface_point(rup.points[2]), surface_point(rup.points[3]))
        backward_seg = LineSegment(surface_point(rup.points[4]), surface_point(rup.points[1]))
    end

    if left_of_line_segment(forward_seg, p) # point is off rupture forward along strike
        return distance_to_line_containing_segment(forward_seg, p)
    elseif left_of_line_segment(backward_seg, p) # point is off rupture backward along strike
        return distance_to_line_containing_segment(backward_seg, p)
    else # the point must be within the rupture
        return 0.0
    end
end

function strike_parallel_distance(rup::Rupture{T}, points::Vector{Site{T}}) where T<:Real
    n = length(points)
    r_ys = zeros(T, n)
    for i in 1:n
        r_ys[i] = strike_parallel_distance(rup, points[i])
    end
    return r_ys
end

strike_parallel_distance(rup::Rupture{T}, site::Site{T}) where T<:Real = strike_parallel_distance(rup, site.position)
strike_parallel_distance(rup::Rupture{T}, sites::Vector{Site{T}}) where T<:Real = strike_parallel_distance(rup, map(s -> s.position, sites))
