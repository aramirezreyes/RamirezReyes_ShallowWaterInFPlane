# This scripts is an example of using the gaussian initialization style

parameters = Dict(
    "architecture" => "GPU",
    "output_filename" => "gaussian_packet_gpu",
    "f" => 0.0,
    "g" => 10.0,
    "convection_timescale" => 1000.0,
    "convection_critical_height" => 40.0,
    "heating_amplitude" => 0.0,
    "large_scale_forcing" => 0.0,
    "convective_radius" => 10_000,
    "relaxation_parameter" => 0.0,
    "relaxation_height" => nothing,
    "Lx" => 600_000,                                                                    
    "Ly" => 600_000,  
    "Lz" => 39.0,
    "Nx" => 200,
    "Ny" => 200,
    "boundary_layer" => true,
    "initialization_style" => "gaussian",
    "initialization_amplitude" => 0.1,
    "gaussian_sigma_x" => 100_000,
    "gaussian_sigma_y" =>20_000,
    "gaussian_rotation" => -pi/4,
    "restart" => false,
    "checkpoint_interval_in_seconds" => Inf,
    "simulation_length_in_days" => 1,
    "output_interval_in_seconds" => 20,
    "timestep_in_seconds" => 20,
)

using RamirezReyes_ShallowWaterInFPlane

run_shallow_simulation(parameters)
