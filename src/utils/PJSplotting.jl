# generic tools for plotting various types

include("PJSpoint.jl")
include("PJSrupture.jl")


using Plots
pyplot()
theme(:juno)

import Plots.plot
import Plots.plot3d!

function plot2d(points::Vector{Point{T}}; kwargs...) where T<:Real
    n = length(points)
    xp = zeros(T, n)
    yp = zeros(T, n)
    for i in 1:n
        xp[i] = points[i].x
        yp[i] = points[i].y
    end

    plot(xp, yp, lab="", xlabel="X Coordinate", ylabel="Y Coordinate", color=[:white], marker=:circle; kwargs...)
end

function plot2d!(points::Vector{Point{T}}; kwargs...) where T<:Real
    n = length(points)
    xp = zeros(T, n)
    yp = zeros(T, n)
    for i in 1:n
        xp[i] = points[i].x
        yp[i] = points[i].y
    end

    plot!(xp, yp, marker=:circle; kwargs...)
end


function plot3d(points::Vector{Point{T}}; kwargs...) where T<:Real
    n = length(points)
    xp = zeros(T, n)
    yp = zeros(T, n)
    zp = zeros(T, n)
    for i in 1:n
        xp[i] = points[i].x
        yp[i] = points[i].y
        zp[i] = points[i].z
    end

    plot(xp, yp, zp, lab="", xlabel="X Coordinate", ylabel="Y Coordinate", zlabel="Z Coordinate", color=[:white], seriestype = :path3d, marker=:circle; kwargs...)
end

function plot3d!(points::Vector{Point{T}}; kwargs...) where T<:Real
    n = length(points)
    xp = zeros(T, n)
    yp = zeros(T, n)
    zp = zeros(T, n)
    for i in 1:n
        xp[i] = points[i].x
        yp[i] = points[i].y
        zp[i] = points[i].z
    end

    plot!(xp, yp, zp, marker=:circle; kwargs...)
end

ps = [ Point(), Point(1.0, 0.0, 3.0), Point(1.0, 1.0, 3.0), Point(0.0, 1.0, 0.0), Point() ]

plot2d(ps)
plot3d(ps)

function plot2d(rup::Rupture{T}; kwargs...) where T<:Real
    # force the rupture closed
    plot2d([rup.points; rup.points[1]]; kwargs...)
end

function plot3d(rup::Rupture{T}; kwargs...) where T<:Real
    # force the rupture closed
    plot3d([rup.points; rup.points[1]]; kwargs...)
end

plot2d(rup)
plot3d(rup)


# test the plotting and distance calcs

hyp = Point(0.0, 0.0, 6.0)
θ = 45.0
δ = 60.0
L = 20.0
W = 10.0
xL = 0.5
xW = 0.5

rup = rupture_from_hypocentre(hyp, θ, δ, L, W, xL, xW)
rup.points

site = Point(0.0, 5.0)

plot2d(rup, aspect_ratio=:equal)
plot2d!([site], marker=:circle)

rupture_distance(rup, site)
joyner_boore_distance(rup, site)

plot3d(rup, aspect_ratio=:equal)
plot3d!([site], marker=:circle)
