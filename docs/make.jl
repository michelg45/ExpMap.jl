using Documenter

include(joinpath(@__DIR__, "..", "src", "ExpMap.jl"))
using .ExpMap

makedocs(
    sitename = "ExpMap",
    authors  = "M. Géradin",
    modules  = [ExpMap],
    format   = Documenter.LaTeX(
        platform = "native",   # compile with system pdflatex/latexmk
        # platform = "none",   # generate .tex only, compile manually
    ),
    pages = [
        "General Description" => "index.md",

        "Types" => [
            "VEC3"      => "vec3.md",
            "RV3"       => "rv3.md",
            "MAT3"      => "mat3.md",
            "VEC6"      => "vec6.md",
            "MAT6"      => "mat6.md",
            "NodeFrame" => "nodeframe.md",
            "mat6_array" => "mat6_array.md",
        ],

        "SO(3) functions" => [
            "R_SO3"           => "R_SO3.md",
            "invR_SO3"        => "invR_SO3.md",
            "euler_to_RV3"    => "euler_to_RV3.md",
            "frame_on_line"   => "frame_on_line.md",
            "T_SO3"           => "T_SO3.md",
            "DT_SO3"          => "DT_SO3.md",
            "DinvT_SO3"       => "DinvT_SO3.md",
            "invT_SO3"        => "invT_SO3.md",
            "T_functions"     => "T_functions.md",
            "invT_functions"  => "invT_functions.md",
        ],

        "SE(3) functions" => [
            "T_SE3_data"        => "T_SE3_data.md",
            "T_SE3_input"       => "T_SE3_input.md",
            "T_SE3_input_data"  => "T_SE3_input_data.md",
            "T_SE3"             => "T_SE3.md",
            "invT_SE3_input"    => "invT_SE3_input.md",
            "invT_SE3"          => "invT_SE3.md",
            "DT_SE3"            => "DT_SE3.md",
            "DinvT_SE3"         => "DinvT_SE3.md",
        ],

        "Lie algebra / helpers" => [
            "RV3 composition" => "rv3_comp_rule.md",
            "tilde"           => "tilde.md",
            "outerp"          => "outerp.md",
            "exp_SE3"         => "exp_SE3.md",
            "log_SE3"         => "log_SE3.md",
            "sk_SE3"          => "sk_SE3.md",
            "Adj_SE3"         => "Adj_SE3.md",
            "bracket"         => "bracket.md",
        ],

        "Numerical integration" => [
            "gauss_points"        => "gauss_points.md",
            "shape_functions_1D"  => "shape_functions_1D.md",
            "shape_functions_2D"  => "shape_functions_2D.md",
        ],

        "Frame interpolation" => [
            "frame_interpol_1D"   => "frame_interpol_1D.md",
            "frame_interpol_2D"   => "frame_interpol_2D.md",
        ],

    ],
    doctest   = false,
    checkdocs = :none,
    warnonly  = true,
    remotes   = nothing,
)
