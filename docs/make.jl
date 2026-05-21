using ExpMap
using Documenter

DocMeta.setdocmeta!(ExpMap, :DocTestSetup, :(using ExpMap); recursive=true)

makedocs(;
    modules=[ExpMap],
    authors="Michel Geradin <mgeradin@gmail.com> and contributors",
    sitename="ExpMap.jl",
    format=Documenter.HTML(;
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
