using GroundMotionModels
using Test

@testset "GroundMotionModels.jl" begin
    # Write your tests here.

    @testset "Velocity Profile Tests" begin

        zi = [15.0, 15.0, 5.0, 10.0]
        vi = [200.0, 400.0, 600.0, 1000.0]

        profile = VelocityProfile(zi, vi)

        @test vs30(profile) ≈ 30.0 / (zi[1] / vi[1] + zi[2] / vi[2])

        zi = [15.0, 5.0, 5.0, 10.0]
        vi = [200.0, 400.0, 600.0, 1000.0]

        profile = VelocityProfile(zi, vi)

        @time depth_to_velocity_horizon(profile, 1500.0)

        @time z1p0(profile)

        vp = VelocityProfile([5.0, 10.0, 20.0], [500.0, 700.0, 1500.0])
        z1p0(vp)
        z2p5(vp)

    end

    @testset "Rupture tests" begin

        rup = Rupture(6.0, [Point()], 0.0, NaN, NaN, Point(NaN, NaN, NaN))

        add_point!(rup, Point(1.0, 0.0))
        rup

        mechanism(rup)

        hyp = Point(0.0, 0.0, 10.0)
        rup = rupture_from_hypocentre(hyp, 90.0, 90.0, 10.0, 10.0, 0.5, 0.5)

        W = rupture_width(6.0, 0.0)
        @code_warntype rupture_from_hypocentre_with_source_scaling(Point(0.0, 0.0, 5.0), 90.0, 90.0, 6.0, 0.0, 0.5, 0.5)

        p = Point(10.0, 5.0, 0.0)
        @code_warntype rupture_distance(rup, p)
        @code_warntype rupture_distance(rup, [Point(0.0, 5.0, 0.0), Point(0.0, 10.0, 0.0)])

        @time joyner_boore_distance(rup, p)

        @time strike_distance(rup, p)

    end

    @testset "Boore et al. (2014)" begin

        gm = PJSbssa2014(0.2, 5.75, 15.0, 1, 0, 0, 0, 350.0, 0.5, "California")
        # @time PJSbssa2014refPGA( 7.0, 10.0, 1, 0, 0, 0 )

    end


    @testset "Chiou & Youngs (2014)" begin
        T = 0.01
        M = 7.36
        Rrup = 117.75
        Rjb = 114.62
        Rx = 120.03
        Vs30 = 316.46
        Frv = 1
        Fnm = 0
        Fhw = 1
        EZ1p0 = exp(-7.15 / 4 * log((Vs30^4 + 571.0^4) / (1360.0^4 + 571.0^4)))
        Z1p0 = 260.0
        ΔZ1p0 = Z1p0 - EZ1p0
        ΔZtor = 0.0
        Dip = 75.0
        Vs30meas = 1

        gm = PJScy2014(T, M, Rrup, Rjb, Rx, ΔZtor, Dip, Fnm, Frv, Fhw, Vs30, Vs30meas, ΔZ1p0)
        gm.IM
        gm.σ

        T = 0.0
        M = 4.0
        Rrup = 5.0
        Rjb = 5.0
        Rx = 5.0
        Frv = 0
        Fnm = 0
        Fhw = 0
        ΔZtor = 0.0
        Dip = 90.0
        Vs30 = 1000.0
        Vs30meas = 1
        ΔZ1p0 = 0.0

        gm = PJScy2014(T, M, Rrup, Rjb, Rx, ΔZtor, Dip, Fnm, Frv, Fhw, Vs30, Vs30meas, ΔZ1p0)
        gm.IM

        gm = PJScy2014(T, M, Rrup, Rjb, Rx, ΔZtor, Dip, Fnm, Frv, Fhw, 300.0, Vs30meas, ΔZ1p0)
        gm.IM
        gm = PJScy2014(T, M, Rrup, Rjb, Rx, ΔZtor, Dip, Fnm, Frv, Fhw, 1000.0, Vs30meas, ΔZ1p0)
        gm.IM



    end

    @testset "Akkar & Bommer (2010)" begin

        abClass = PJSab2010(0.0, 7.0, 10.0, 0, 1, 0, 0)
        abVs30 = PJSab2010(0.0, 7.0, 10.0, 200.0, 0, 0)

        @test abClass.lnIM ≈ abVs30.lnIM

        abClass = PJSab2010(0.01, 7.0, 10.0, 1, 0, 0, 0)
        abVs30 = PJSab2010(0.01, 7.0, 10.0, 500.0, 0, 0)

        @test abClass.lnIM ≈ abVs30.lnIM

        abClass = PJSab2010(0.0, 7.0, 10.0, 0, 0, 0, 0)
        abVs30 = PJSab2010(0.0, 7.0, 10.0, 1000.0, 0, 0)

        @test abClass.lnIM ≈ abVs30.lnIM

        abNull = PJSab2010(0.0, 7.0, 10.0, 1, 1, 0, 0)
        @test isnan(abNull.lnIM)

        abNull = PJSab2010(0.0, 7.0, 10.0, 0, 0, 1, 1)
        @test isnan(abNull.lnIM)

    end

    @testset "Akkar et al. (2014)" begin
        # asb = PJSasb2014(0.0, 7.0, 10.0, 200.0, 0, 0, "Rjb")
        # asb = PJSasb2014(0.0, 7.0, 10.0, 350.0, 0, 0, "Repi")

    end

end
