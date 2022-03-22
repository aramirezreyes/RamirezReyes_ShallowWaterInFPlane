module RamirezReyes_ShallowWaterInFPlane

using Oceananigans: CPU, GPU, @cuda
include(joinpath(@__DIR__,"../src/convectiveparameterization.jl"))
include(joinpath(@__DIR__,"../src/arrayutils.jl"))

export update_convective_events!, model_forcing, u_damping, v_damping


end #module
