
using LinearAlgebra

"""
    LineSegment{T<:Real}

Custom type representing a line segment connecting two points
"""
struct LineSegment{T<:Real}
    p1::Point{T}
    p2::Point{T}
end


"""
    distance_to_line_segment(seg::LineSegment{T}, p::Point{T}) where T<:Real

Function either finds the distance from (x0,y0) to its projection onto
the segment (x1,y1)-(x2,y2) if it exists, or otherwise finds the
distance to the closer of (x1,y1) and (x2,y2)
"""
function distance_to_line_segment(seg::LineSegment{T}, p::Point{T}) where T<:Real
    if dot(seg.p1-seg.p2, p-seg.p2)*dot(seg.p2-seg.p1, p-seg.p1) >= 0.0
        A = [ seg.p1.x seg.p1.y 1.0;
              seg.p2.x seg.p2.y 1.0;
              p.x p.y 1.0 ]
        d_ab = magnitude(seg.p1-seg.p2)
        d = abs(det(A))/d_ab
    else
        d_as = magnitude(seg.p1-p)
        d_bs = magnitude(seg.p2-p)
        d = min(d_as, d_bs)
    end
    return d
end

"""
distance_to_line_containing_segment(seg::LineSegment{T}, p::Point{T}) where T<:Real

Distance to the line that includes a given line segment
"""
function distance_to_line_containing_segment(seg::LineSegment{T}, p::Point{T}) where T<:Real
	if seg.p2.x == seg.p1.x
		return abs(p.x - seg.p1.x)
	else
		slope = (seg.p2.y - seg.p1.y) / (seg.p2.x - seg.p1.x)
		inter = seg.p1.y - slope*seg.p1.x
		a = -slope
		b = 1.0
		c = -inter
		d = abs(a*p.x + b*p.y + c) / sqrt(a^2 + b^2)
		return d
	end
end

# seg = LineSegment(Point(0.0,0.0), Point(1.0,1.0))
# p = Point(-2.0,2.0)
# distance_to_line_containing_segment(seg, p)

"""
left_of(seg::LineSegment{T}, p::Point{T}) where T<:Real

Determine whether the point p is 'left of' the line segment. Looking from point 1
to point 2 of the line segment, is the point p on the left?
"""
function left_of_line_segment(seg::LineSegment{T}, p::Point{T}) where T<:Real
    double_signed_area = (seg.p2.x - seg.p1.x)*(p.y - seg.p1.y) - (seg.p2.y - seg.p1.y)*(p.x - seg.p1.x)
	return ( double_signed_area >= 0.0 ) ? true : false
end
