# In this script we disable gravity, convection, large scale forcing and coriolis to oberve only the damping terms

parameters = Dict(
    "architecture" => "GPU",
    "output_filename" => "only_damping_ft_gpu",
    "f" => 0.0,
    "g" => 0.0,
    "convection_timescale" => 1000.0,
    "convection_critical_height" => 40.0,
    "heating_amplitude" => 0.0,
    "large_scale_forcing" => 0.0,
    "convective_radius" => 10_000,
    "relaxation_parameter" => 1/1000.0,
    "relaxation_height" => nothing,
    "Lx" => 30000,
    "Ly" => 30000,
    "Lz" => 40.0,
    "Nx" => 100,
    "Ny" => 100,
    "boundary_layer" => false,
    "initialization_style" => "one_convecting_point",
    "initialization_amplitude" => 1,
    "simulation_length_in_days" => 2000 / 86400,
    "output_interval_in_seconds" => 20,
    "timestep_in_seconds" => 20,
    "restart" => false,
    "checkpoint_interval_in_seconds" => Inf,
)

using RamirezReyes_ShallowWaterInFPlane

run_shallow_simulation(parameters)
