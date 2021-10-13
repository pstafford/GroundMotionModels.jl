
using Statistics
using StatsBase
using Plots
pyplot()
# theme(:juno)
theme(:default)

include("PJSwellsCoppersmith1994.jl")
include("PJSgroundMotionModels.jl")

# magnitude range
mi = collect(range(5.0, stop=8.0, step=0.1))
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
Zhyp = EZtor + 2/3*Wi * sind(δ)

# create the correspoding hypocentral locations
hypocentres = map( z -> Point(0.0, 0.0, z), Zhyp )

# create the corresponding rupture scenarios
ruptures = map(i -> rupture_from_hypocentre_with_source_scaling(hypocentres[i], θ, δ, mi[i], λ, 0.5, 1/3), 1:length(mi))


# specify line of distances
r_epi = [ 5.0, 50.0, 100.0, 150.0, 200.0, 400.0, 650.0 ]
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

tid = 1
sid = 1


function plot_magnitude(tid, sid, outputpath::String="Figures/Magnitude/")
    # NGA West 2 models
    Sai_ask14 = zeros(size(mi))
    Sai_bssa14 = zeros(size(mi))
    Sai_cb14 = zeros(size(mi))
    Sai_cy14 = zeros(size(mi))
    Sai_idriss14 = zeros(size(mi))
    for i in 1:length(mi)
        Sai_ask14[i] = PJSask2014(Ti[tid], ruptures[i], sites[sid]).IM
        Sai_bssa14[i] = PJSbssa2014(Ti[tid], ruptures[i], sites[sid]).IM
        Sai_cb14[i] = PJScb2014(Ti[tid], ruptures[i], sites[sid]).IM
        Sai_cy14[i] = PJScy2014(Ti[tid], ruptures[i], sites[sid]).IM
        Sai_idriss14[i] = PJSidriss2014(Ti[tid], ruptures[i], sites[sid]).IM
    end

    # European models
    Sai_asb14_rjb = zeros(size(mi))
    Sai_asb14_repi = zeros(size(mi))
    Sai_asb14_rhyp = zeros(size(mi))
    Sai_bindi14_rjb = zeros(size(mi))
    Sai_bindi14_rhyp = zeros(size(mi))
    Sai_kbc16_regional = zeros(size(mi))
    Sai_kbc16_others = zeros(size(mi))
    Sai_kbc16_italy = zeros(size(mi))
    Sai_kbc16_turkey = zeros(size(mi))
    for i in 1:length(mi)
        Sai_asb14_rjb[i] = PJSasb2014(Ti[tid], ruptures[i], sites[sid], "Rjb").IM
        Sai_asb14_repi[i] = PJSasb2014(Ti[tid], ruptures[i], sites[sid], "Repi").IM
        Sai_asb14_rhyp[i] = PJSasb2014(Ti[tid], ruptures[i], sites[sid], "Rhyp").IM
        Sai_bindi14_rjb[i] = PJSbindi2014(Ti[tid], ruptures[i], sites[sid], "Rjb").IM
        Sai_bindi14_rhyp[i] = PJSbindi2014(Ti[tid], ruptures[i], sites[sid], "Rhyp").IM
        Sai_kbc16_regional[i] = PJSkbc2016(Ti[tid], ruptures[i], sites[sid], "Regional").IM
        Sai_kbc16_others[i] = PJSkbc2016(Ti[tid], ruptures[i], sites[sid], "Others").IM
        Sai_kbc16_italy[i] = PJSkbc2016(Ti[tid], ruptures[i], sites[sid], "Italy").IM
        Sai_kbc16_turkey[i] = PJSkbc2016(Ti[tid], ruptures[i], sites[sid], "Turkey").IM
    end


    # Spanish models
    Sai_esp_LLC = zeros(size(mi))
    Sai_esp_LCC = zeros(size(mi))
    Sai_esp_LUC = zeros(size(mi))
    Sai_esp_CLC = zeros(size(mi))
    Sai_esp_CCC = zeros(size(mi))
    Sai_esp_CUC = zeros(size(mi))
    Sai_esp_ULC = zeros(size(mi))
    Sai_esp_UCC = zeros(size(mi))
    Sai_esp_UUC = zeros(size(mi))

    for i in 1:length(mi)
        Sai_esp_LLC[i] = PJSspainCrustal(Ti[tid], ruptures[i], sites[sid], "lower", "cy14", "central").IM
        Sai_esp_LCC[i] = PJSspainCrustal(Ti[tid], ruptures[i], sites[sid], "lower", "1p15", "central").IM
        Sai_esp_LUC[i] = PJSspainCrustal(Ti[tid], ruptures[i], sites[sid], "lower", "free", "central").IM
        Sai_esp_CLC[i] = PJSspainCrustal(Ti[tid], ruptures[i], sites[sid], "central", "cy14", "central").IM
        Sai_esp_CCC[i] = PJSspainCrustal(Ti[tid], ruptures[i], sites[sid], "central", "1p15", "central").IM
        Sai_esp_CUC[i] = PJSspainCrustal(Ti[tid], ruptures[i], sites[sid], "central", "free", "central").IM
        Sai_esp_ULC[i] = PJSspainCrustal(Ti[tid], ruptures[i], sites[sid], "upper", "cy14", "central").IM
        Sai_esp_UCC[i] = PJSspainCrustal(Ti[tid], ruptures[i], sites[sid], "upper", "1p15", "central").IM
        Sai_esp_UUC[i] = PJSspainCrustal(Ti[tid], ruptures[i], sites[sid], "upper", "free", "central").IM
    end

    Sai_esp_all = [ Sai_esp_LLC Sai_esp_LCC Sai_esp_LUC Sai_esp_CLC Sai_esp_CCC Sai_esp_CUC Sai_esp_ULC Sai_esp_UCC Sai_esp_UUC ]
    wi = [ 0.3*0.4, 0.3*0.3, 0.3*0.3, 0.4*0.4, 0.4*0.3, 0.4*0.3, 0.3*0.4, 0.3*0.3, 0.3*0.3 ]

    Sai_esp_med = median(Sai_esp_all; dims=2)
    Sai_esp_max = maximum(Sai_esp_all; dims=2)
    Sai_esp_min = minimum(Sai_esp_all; dims=2)
    Sai_esp_sd = std(log.(Sai_esp_all); dims=2)
    Sai_esp_wmu = similar(Sai_esp_sd)
    Sai_esp_wsd = similar(Sai_esp_sd)
    for i in 1:length(mi)
        Sai_esp_wmu[i] = exp.(mean(log.(Sai_esp_all[i,:]), weights(wi)))
        Sai_esp_wsd[i] = std(log.(Sai_esp_all[i,:]), weights(wi), corrected=false)
    end

    dSaiu_esp = Sai_esp_max .- Sai_esp_med
    dSail_esp = Sai_esp_med .- Sai_esp_min
    dSai_esp = Sai_esp_max .- Sai_esp_min
    rSai_esp = Sai_esp_max ./ Sai_esp_min

    Sai_nga_all = [ Sai_ask14 Sai_bssa14 Sai_cb14 Sai_cy14 Sai_idriss14 ]

    Sai_nga_med = median(Sai_nga_all; dims=2)
    Sai_nga_max = maximum(Sai_nga_all; dims=2)
    Sai_nga_min = minimum(Sai_nga_all; dims=2)
    Sai_nga_mu = exp.(mean(log.(Sai_nga_all); dims=2))
    Sai_nga_sd = std(log.(Sai_nga_all); dims=2)

    dSaiu_nga = Sai_nga_max .- Sai_nga_med
    dSail_nga = Sai_nga_med .- Sai_nga_min
    dSai_nga = Sai_nga_max .- Sai_nga_min
    rSai_nga = Sai_nga_max ./ Sai_nga_min


    Sai_eur_all = [ Sai_asb14_rjb Sai_asb14_repi Sai_asb14_rhyp Sai_bindi14_rjb Sai_bindi14_rhyp Sai_kbc16_regional Sai_kbc16_others Sai_kbc16_italy Sai_kbc16_turkey ]

    Sai_eur_med = median(Sai_eur_all; dims=2)
    Sai_eur_max = maximum(Sai_eur_all; dims=2)
    Sai_eur_min = minimum(Sai_eur_all; dims=2)
    Sai_eur_mu = exp.(mean(log.(Sai_eur_all); dims=2))
    Sai_eur_sd = std(log.(Sai_eur_all); dims=2)

    dSaiu_eur = Sai_eur_max .- Sai_eur_med
    dSail_eur = Sai_eur_med .- Sai_eur_min
    dSai_eur = Sai_eur_max .- Sai_eur_min
    rSai_eur = Sai_eur_max ./ Sai_eur_min

    theme(:default)
    p1 = plot(mi, log.(Sai_esp_wmu), ribbon=Sai_esp_wsd, lab="Spain Crustal", box=:true, color=[1], linewidth=2)
    plot!(mi, log.(Sai_nga_mu), ribbon=Sai_nga_sd, lab="NGA West2", color=[2], linewidth=2)
    plot!(mi, log.(Sai_eur_mu), ribbon=Sai_eur_sd, lab="European", color=[3], linewidth=2)
    xlabel!("Magnitude")
    ylabel!("Log Spectral acceleration")
    title!("T = $(Ti[tid]), Repi = $(r_epi[sid])km")

    p2 = plot(mi, Sai_esp_wsd, lab="Spain Crustal", box=:true, color=[1], linewidth=2)
    plot!(mi, Sai_nga_sd, lab="NGA West 2", color=[2], linewidth=2)
    plot!(mi, Sai_eur_sd, lab="European", color=[3], linewidth=3)
    ylims!(0,0.8)
    xlabel!("Magnitude")
    ylabel!("Between-model standard deviation")
    title!("T = $(Ti[tid]), Repi = $(r_epi[sid])km")

    l = @layout [a b]
    p = plot(p1, p2, layout=l)
    plot!(size=(980,580))
    outfile = string(outputpath,"PJSmagnitude_T$(replace(string(Ti[tid]),"."=>"p"))_Repi$(replace(string(r_epi[sid]),"."=>"p"))_v1.pdf")
    savefig(outfile)
    return p
end

for i in 1:length(Ti)
    for j in 1:length(r_epi)
        plot_magnitude(i,j)
    end
end
