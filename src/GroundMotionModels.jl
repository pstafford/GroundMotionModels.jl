
module GroundMotionModels

using Distributions
using Roots
using Interpolations
using LinearAlgebra


export PJSgroundMotion,
        Point,
        Rupture,
        add_point!,
        dip,
        epicentral_distance,
        hypocentral_distance,
        joyner_boore_distance,
        mechanism,
        rupture_width,
        rupture_length,
        rupture_from_hypocentre,
        rupture_from_hypocentre_with_source_scaling,
        rupture_distance,
        strike,
        strike_distance,
        strike_parallel_distance,
        zbot,
        zhyp,
        ztor,
        VelocityProfile,
        Site,
        vs30,
        depth_to_velocity_horizon,
        z1p0,
        z2p5,
        PJSab2010,
        PJSasb2014,
        PJSask2014,
        PJScb2014,
        PJScy2014,
        PJSbindi2014,
        PJSbssa2014,
        PJSkbc2016,
        PJSidriss2014,
        PJSsadigh1997,
        PJSspainCrustal


# Geometric utilities
include("utils/PJSvector.jl")
include("utils/PJSpoint.jl")
include("utils/PJStriangle.jl")
include("utils/PJSlineSegment.jl")

# Source-scaling relations (used in rupture generation)
# TODO: make a separate module for source-scaling relations
include("utils/PJSwellsCoppersmith1994.jl")


# GMM utilities
include("utils/PJSvelocityProfile.jl")
include("utils/PJSsite.jl")
include("utils/PJSrupture.jl")
include("utils/PJSgroundMotions.jl")

# Spain NPP Crustal
include("crustal/PJSspanishCrustalGMM_v1.jl")

# NGA West 2
include("crustal/PJSabrahamsonSilvaKamai2014.jl")
include("crustal/PJSbooreStewartSeyhanAtkinson2014.jl")
include("crustal/PJScampbellBozorgnia2014.jl")
include("crustal/PJSchiouYoungs2014.jl")
include("crustal/PJSidriss2014.jl")

# RESORCE
include("crustal/PJSakkarSandikayyaBommer2014.jl")
include("crustal/PJSbindi2014.jl")
include("crustal/PJSkothaBindiCotton2016.jl")

# European
include("crustal/PJSakkarBommer2010.jl")

# SRL 1997
include("crustal/PJSsadigh1997.jl")

end
