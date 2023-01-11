using RamirezReyes_ShallowWaterInFPlane

parameters_dict = Dict(
    "architecture" => "GPU",
    "output_filename" => "rankine_vortex_run_GPU",
    "f" => 1e-5,# coriolis parameter5e-4 #5e-4
    "g" => 9.8,#gravitational acceleration 9.8
    "convection_timescale" => 28800.0, #convective time scale
    "convection_critical_height" => 130.0, #convection_triggering height
    "heating_amplitude" => 1.0e9,#convective heating 1.0e9 #originally 9 for heating, -8 for cooling convective heating
    "large_scale_forcing" => (1.12 / 3) * 5.0e-5,#radiative cooling rate
    "convective_radius" => 30_000, #convective radius
    "relaxation_parameter" => 1.0 / (0.8 * 86400), #relaxation timescale
    "relaxation_height" => nothing, #height to relax to, if nothing, h relax to its mean
    "Lx" => 4e6, #domain
    "Ly" => 4e6,
    "Lz" => 126.5,
    "Nx" => 500,
    "Ny" => 500,
    "boundary_layer" => true, #if true Lz=30 otherwise Lz = 41
    "initialization_style" => "rankine_vortex",
    "rankine_radius" => 200_000,
    "rankine_amplitude" => 5.0,
    "rankine_center" => (2_000_000.0, 2_000_000.0),
    "initialization_amplitude" => 4.0,
    "simulation_length_in_days" => 20,
    "output_interval_in_seconds" => 7200,
    "timestep_in_seconds" => 40.0,
)

run_shallow_simulation(parameters_dict)