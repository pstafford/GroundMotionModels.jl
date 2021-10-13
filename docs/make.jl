using GroundMotionModels
using Documenter

DocMeta.setdocmeta!(GroundMotionModels, :DocTestSetup, :(using GroundMotionModels); recursive=true)

makedocs(;
    modules=[GroundMotionModels],
    authors="Peter Stafford <p.stafford@me.com>",
    repo="https://github.com/pstafford/GroundMotionModels.jl/blob/{commit}{path}#{line}",
    sitename="GroundMotionModels.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://pstafford.github.io/GroundMotionModels.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/pstafford/GroundMotionModels.jl",
    devbranch="main",
)
