# In tis script we initialize the simulation with a rankine vortex. It did not lead to anything interesting yet.

parameters_dict = Dict(
    "architecture" => "GPU",
    "output_filename" => "rankine_vortex_run_gpu",
    "f" => 1e-5,
    "g" => 9.8,
    "convection_timescale" => 28800.0,
    "convection_critical_height" => 130.0,
    "heating_amplitude" => 1.0e9,
    "large_scale_forcing" => (1.12 / 3) * 5.0e-5,
    "convective_radius" => 30_000,
    "relaxation_parameter" => 1.0 / (0.8 * 86400),
    "relaxation_height" => nothing,
    "Lx" => 4e6,
    "Ly" => 4e6,
    "Lz" => 126.5,
    "Nx" => 500,
    "Ny" => 500,
    "boundary_layer" => true,
    "initialization_style" => "rankine_vortex",
    "rankine_radius" => 200_000,
    "rankine_amplitude" => 5.0,
    "rankine_center" => (2_000_000.0, 2_000_000.0),
    "initialization_amplitude" => 4.0,
    "simulation_length_in_days" => 20,
    "output_interval_in_seconds" => 7200,
    "timestep_in_seconds" => 40.0,
)

using RamirezReyes_ShallowWaterInFPlane

run_shallow_simulation(parameters_dict)
