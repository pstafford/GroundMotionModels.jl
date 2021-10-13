
using Statistics
using Plots
pyplot()
theme(:juno)

include("PJSwellsCoppersmith1994.jl")
include("PJSgroundMotionModels.jl")

# magnitude range
mi = collect(range(5.0, stop=7.5, step=0.5))
# start with strike-slip mechanisms
# strike, dip and rake
θ = 0.0
δ = 90.0
λ = 0.0

# mechanism flag
Frv = 0
Fnm = 0
Fss = 1

# use the CY14 expected depths, define hypcentral depths to match the
if Frv == 1
    EZtor = @. max(2.704 - 1.226*max(mi-5.849,0.0),0.0)^2
else
    EZtor = @. max(2.673 - 1.136*max(mi-4.970,0.0),0.0)^2
end

# expected widths
Wi = rupture_width.(mi, λ)

# adjust the expected depth by 2/3 of the rupture width (vertically)
# Zhyp = EZtor + 2/3*Wi * sind(δ)
Zhyp = collect(range(0.0, stop=20.0, step=0.25))

# plot(mi, Zhyp, lab="Zhyp", box=:true)
# plot!(mi, EZtor, lab="EZtor")
# plot!(mi, 2/3*Wi*sind(δ), lab="ΔW (up-dip)")
# xlabel!("Magnitude")
# ylabel!("Depth & Up-Dip Width")


# create the correspoding hypocentral locations
hypocentres = map( z -> Point(0.0, 0.0, z), Zhyp )

# create the corresponding rupture scenarios
ruptures = [ rupture_from_hypocentre_with_source_scaling(hypocentres[i], θ, δ, mi[j], λ, 0.5, 1/3) for i in 1:length(Zhyp), j in 1:length(mi) ]

# specify line of distances
# r_epi = exp10.(range(0.0, stop=3.0, step=0.05))
r_epi = [ 10.0, 70.0 ]
# relative strike (azimuth)
Δθ = 90.0

# define the x/y coordinates
x_epi = r_epi * sind(θ+Δθ)
y_epi = r_epi * cosd(θ+Δθ)
z_epi = zeros(size(r_epi))

# create the corresponding points
site_locations = [ Point(x_epi[i],y_epi[i],z_epi[i]) for i in 1:length(r_epi) ]
# create the particular sites (including the velocity profile)
# use key metrics from the Boore (2016) velocity profile
velocity_profile = VelocityProfile([30.0, 620.0, 7200.0, 10_000.0], [760.0, 1000.0, 2500.0, 3500.0])

# create the sites
sites = map(p -> Site(p, velocity_profile, 1, "California"), site_locations)

# specify periods and associated strings
Ti = [ 0.01, 0.03, 0.1, 0.3, 1.0, 3.0 ]
# Ti_str = [ "0p01", "0p02", "0p03", "0p04", "0p05", "0p075", "0p10", "0p15", "0p20", "0p30", "0p40", "0p50", "0p75", "1p00", "1p50", "2p00", "3p00", "5p00", "7p50", "10p0" ]

tid = 1 # period
mid = 3 # magnitude
sid = 1 # site location

rupture_distance(ruptures[1,mid],sites[sid])
ztor(ruptures[1,mid])
EZtor[mid]
ΔZtor = ztor(ruptures[1,mid]) - EZtor[mid]


function plot_depth(tid, mid, sid, outputpath::String="Figures/Depth/")

    Sai_ask14 = zeros(size(Zhyp))
    Sai_bssa14 = zeros(size(Zhyp))
    Sai_cb14 = zeros(size(Zhyp))
    Sai_cy14 = zeros(size(Zhyp))
    Sai_idriss14 = zeros(size(Zhyp))
    for i in 1:length(Zhyp)
        Sai_ask14[i] = PJSask2014(Ti[tid], ruptures[i,mid], sites[sid]).IM
        Sai_bssa14[i] = PJSbssa2014(Ti[tid], ruptures[i,mid], sites[sid]).IM
        Sai_cb14[i] = PJScb2014(Ti[tid], ruptures[i,mid], sites[sid]).IM
        Sai_cy14[i] = PJScy2014(Ti[tid], ruptures[i,mid], sites[sid]).IM
        Sai_idriss14[i] = PJSidriss2014(Ti[tid], ruptures[i,mid], sites[sid]).IM
    end


    # European models
    Sai_asb14_rjb = zeros(size(Zhyp))
    Sai_asb14_repi = zeros(size(Zhyp))
    Sai_asb14_rhyp = zeros(size(Zhyp))
    Sai_bindi14_rjb = zeros(size(Zhyp))
    Sai_bindi14_rhyp = zeros(size(Zhyp))
    Sai_kbc16_regional = zeros(size(Zhyp))
    Sai_kbc16_others = zeros(size(Zhyp))
    Sai_kbc16_italy = zeros(size(Zhyp))
    Sai_kbc16_turkey = zeros(size(Zhyp))
    for i in 1:length(Zhyp)
        Sai_asb14_rjb[i] = PJSasb2014(Ti[tid], ruptures[i,mid], sites[sid], "Rjb").IM
        Sai_asb14_repi[i] = PJSasb2014(Ti[tid], ruptures[i,mid], sites[sid], "Repi").IM
        Sai_asb14_rhyp[i] = PJSasb2014(Ti[tid], ruptures[i,mid], sites[sid], "Rhyp").IM
        Sai_bindi14_rjb[i] = PJSbindi2014(Ti[tid], ruptures[i,mid], sites[sid], "Rjb").IM
        Sai_bindi14_rhyp[i] = PJSbindi2014(Ti[tid], ruptures[i,mid], sites[sid], "Rhyp").IM
        Sai_kbc16_regional[i] = PJSkbc2016(Ti[tid], ruptures[i,mid], sites[sid], "Regional").IM
        Sai_kbc16_others[i] = PJSkbc2016(Ti[tid], ruptures[i,mid], sites[sid], "Others").IM
        Sai_kbc16_italy[i] = PJSkbc2016(Ti[tid], ruptures[i,mid], sites[sid], "Italy").IM
        Sai_kbc16_turkey[i] = PJSkbc2016(Ti[tid], ruptures[i,mid], sites[sid], "Turkey").IM
    end


    # Spanish models
    Sai_esp_LLC = zeros(size(Zhyp))
    Sai_esp_LCC = zeros(size(Zhyp))
    Sai_esp_LUC = zeros(size(Zhyp))
    Sai_esp_CLC = zeros(size(Zhyp))
    Sai_esp_CCC = zeros(size(Zhyp))
    Sai_esp_CUC = zeros(size(Zhyp))
    Sai_esp_ULC = zeros(size(Zhyp))
    Sai_esp_UCC = zeros(size(Zhyp))
    Sai_esp_UUC = zeros(size(Zhyp))

    for i in 1:length(Zhyp)
        Sai_esp_LLC[i] = PJSspainCrustal(Ti[tid], ruptures[i,mid], sites[sid], "lower", "lower", "central").IM
        Sai_esp_LCC[i] = PJSspainCrustal(Ti[tid], ruptures[i,mid], sites[sid], "lower", "central", "central").IM
        Sai_esp_LUC[i] = PJSspainCrustal(Ti[tid], ruptures[i,mid], sites[sid], "lower", "upper", "central").IM
        Sai_esp_CLC[i] = PJSspainCrustal(Ti[tid], ruptures[i,mid], sites[sid], "central", "lower", "central").IM
        Sai_esp_CCC[i] = PJSspainCrustal(Ti[tid], ruptures[i,mid], sites[sid], "central", "central", "central").IM
        Sai_esp_CUC[i] = PJSspainCrustal(Ti[tid], ruptures[i,mid], sites[sid], "central", "upper", "central").IM
        Sai_esp_ULC[i] = PJSspainCrustal(Ti[tid], ruptures[i,mid], sites[sid], "upper", "lower", "central").IM
        Sai_esp_UCC[i] = PJSspainCrustal(Ti[tid], ruptures[i,mid], sites[sid], "upper", "central", "central").IM
        Sai_esp_UUC[i] = PJSspainCrustal(Ti[tid], ruptures[i,mid], sites[sid], "upper", "upper", "central").IM
    end

    Sai_esp_all = [ Sai_esp_LLC Sai_esp_LCC Sai_esp_LUC Sai_esp_CLC Sai_esp_CCC Sai_esp_CUC Sai_esp_ULC Sai_esp_UCC Sai_esp_UUC ]

    Sai_esp_med = median(Sai_esp_all; dims=2)
    Sai_esp_max = maximum(Sai_esp_all; dims=2)
    Sai_esp_min = minimum(Sai_esp_all; dims=2)
    Sai_esp_sd = std(log.(Sai_esp_all); dims=2)

    dSaiu_esp = Sai_esp_max .- Sai_esp_med
    dSail_esp = Sai_esp_med .- Sai_esp_min
    dSai_esp = Sai_esp_max .- Sai_esp_min
    rSai_esp = Sai_esp_max ./ Sai_esp_min

    Sai_nga_all = [ Sai_ask14 Sai_bssa14 Sai_cb14 Sai_cy14 Sai_idriss14 ]

    Sai_nga_med = median(Sai_nga_all; dims=2)
    Sai_nga_max = maximum(Sai_nga_all; dims=2)
    Sai_nga_min = minimum(Sai_nga_all; dims=2)
    Sai_nga_sd = std(log.(Sai_nga_all); dims=2)

    dSaiu_nga = Sai_nga_max .- Sai_nga_med
    dSail_nga = Sai_nga_med .- Sai_nga_min
    dSai_nga = Sai_nga_max .- Sai_nga_min
    rSai_nga = Sai_nga_max ./ Sai_nga_min


    Sai_eur_all = [ Sai_asb14_rjb Sai_asb14_repi Sai_asb14_rhyp Sai_bindi14_rjb Sai_bindi14_rhyp Sai_kbc16_regional Sai_kbc16_others Sai_kbc16_italy Sai_kbc16_turkey ]

    Sai_eur_med = median(Sai_eur_all; dims=2)
    Sai_eur_max = maximum(Sai_eur_all; dims=2)
    Sai_eur_min = minimum(Sai_eur_all; dims=2)
    Sai_eur_sd = std(log.(Sai_eur_all); dims=2)

    dSaiu_eur = Sai_eur_max .- Sai_eur_med
    dSail_eur = Sai_eur_med .- Sai_eur_min
    dSai_eur = Sai_eur_max .- Sai_eur_min
    rSai_eur = Sai_eur_max ./ Sai_eur_min

    theme(:default)
    p1 = plot(Zhyp, Sai_esp_med, ribbon=(dSail_esp, dSaiu_esp), lab="Spain Crustal", box=:true, yaxis=:log10, color=[1], linewidth=2)
    plot!(Zhyp, Sai_esp_all, lab="", color=[1])
    plot!(Zhyp, Sai_nga_med, ribbon=(dSail_nga, dSaiu_nga), lab="NGA West2", color=[2], linewidth=2)
    plot!(Zhyp, Sai_nga_all, lab="", color=[2])
    plot!(Zhyp, Sai_eur_med, ribbon=(dSail_eur, dSaiu_eur), lab="European", color=[3], linewidth=2)
    plot!(Zhyp, Sai_eur_all, lab="", color=[3])
    xlabel!("Hypocentral depth (km)")
    ylabel!("Spectral acceleration (g)")
    title!("T = $(Ti[tid])s, M = $(mi[mid]), Repi = $(r_epi[sid])km")


    p2 = plot(Zhyp, Sai_esp_sd, lab="Spain Crustal", box=:true, color=[1], linewidth=2)
    plot!(Zhyp, Sai_nga_sd, lab="NGA West 2", color=[2], linewidth=2)
    plot!(Zhyp, Sai_eur_sd, lab="European", color=[3], linewidth=3)
    ylims!(0,0.5)
    xlabel!("Hypocentral depth (km)")
    ylabel!("Between-model standard deviation")
    title!("T = $(Ti[tid])s, M = $(mi[mid]), Repi = $(r_epi[sid])km")

    l = @layout [a b]
    p = plot(p1, p2, layout=l)
    plot!(size=(980,580))
    outfile = string(outputpath,"PJSdepth_T$(replace(string(Ti[tid]),"."=>"p"))_M$(replace(string(mi[mid]),"."=>"p"))_R$(replace(string(r_epi[sid]),"."=>"p"))_v0.pdf")
    savefig(outfile)
    return p
end

# @time plot_spectra(1, 1)

for i in 1:length(Ti)
    for j in 1:length(mi)
        for k in 1:length(r_epi)
            plot_depth(i,j,k)
        end
    end
end
