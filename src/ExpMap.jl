module ExpMap

# Two files (rv3_comp_rule, invR_SO3) are not yet duplicated in the
# documenter/src/ directory, so they are loaded from the upstream source tree.
const _src_dir = "/home/michel22/Projects/Julia/divers/src/my_algebra_new/"

using Printf
using LinearAlgebra
import Base: +
import Base: -
import Base: *
import Base: /
import Base: \
import Base: transpose
import LinearAlgebra: diag
import Base: copy
import Base: show
import Base: zero
import Base: one
import Base: getindex
import Base: setindex!
import Base: getproperty
import Base: setproperty!
import Base: inv


include(joinpath(@__DIR__, "VEC3.jl"))
include(joinpath(@__DIR__, "MAT3.jl"))
include(joinpath(@__DIR__, "RV3.jl"))
include(joinpath(@__DIR__, "VEC6.jl"))
include(joinpath(@__DIR__, "MAT6.jl"))
include(joinpath(@__DIR__, "mat6_array.jl"))
include(_src_dir * "rv3_comp_rule.jl")
include(joinpath(@__DIR__, "outerp.jl"))
include(joinpath(@__DIR__, "tilde.jl"))
include(joinpath(@__DIR__, "T_functions.jl"))
include(joinpath(@__DIR__, "invT_functions.jl"))
include(joinpath(@__DIR__, "R_SO3.jl"))
include(joinpath(@__DIR__, "T_SO3.jl"))
include(joinpath(@__DIR__, "invT_SO3.jl"))
include(_src_dir * "invR_SO3.jl")
include(joinpath(@__DIR__, "euler_to_RV3.jl"))
include(joinpath(@__DIR__, "frame_on_line.jl"))
include(joinpath(@__DIR__, "NodeFrame.jl"))
include(joinpath(@__DIR__, "exp_SE3.jl"))
include(joinpath(@__DIR__, "log_SE3.jl"))
include(joinpath(@__DIR__, "sk_SE3.jl"))
include(joinpath(@__DIR__, "Adj_SE3.jl"))
include(joinpath(@__DIR__, "bracket.jl"))
include(joinpath(@__DIR__, "DT_SO3.jl"))
include(joinpath(@__DIR__, "DinvT_SO3.jl"))
include(joinpath(@__DIR__, "T_SE3_input.jl"))
include(joinpath(@__DIR__, "T_SE3_data.jl"))
include(joinpath(@__DIR__, "T_SE3_input_data.jl"))
include(joinpath(@__DIR__, "T_SE3.jl"))
include(joinpath(@__DIR__, "invT_SE3_input.jl"))
include(joinpath(@__DIR__, "invT_SE3.jl"))
include(joinpath(@__DIR__, "DT_SE3.jl"))
include(joinpath(@__DIR__, "DinvT_SE3.jl"))
include(joinpath(@__DIR__, "gauss_points.jl"))
include(joinpath(@__DIR__, "shape_functions_1D.jl"))
include(joinpath(@__DIR__, "shape_functions_2D.jl"))
include(joinpath(@__DIR__, "frame_interpol_1D.jl"))
include(joinpath(@__DIR__, "frame_interpol_2D.jl"))

# ── Public API ────────────────────────────────────────────────────────────────
export VEC3, MAT3, RV3, VEC6, MAT6, NodeFrame
export vec_vec6, vec_mat6, mat_mat6
export R_SO3, invR_SO3, euler_to_RV3, frame_on_line, T_SO3, invT_SO3
export T_functions, invT_functions
export tilde, outerp, dotp, crossp, norm2
export VEC6_unit
export exp_SE3, log_SE3
export sk_SE3, Adj_SE3, bracket
export DT_SO3, DinvT_SO3
export T_SE3_input, T_SE3_data, T_SE3_input_data, T_SE3
export invT_SE3_input, invT_SE3
export DT_SE3, DinvT_SE3
export gauss_points
export shape_functions_1D, shape_functions_2D
export frame_interpol_1D, frame_interpol_2D

end # module ExpMap
