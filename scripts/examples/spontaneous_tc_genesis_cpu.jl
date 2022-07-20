using RamirezReyes_ShallowWaterInFPlane

parameters = Dict(
    "architecture" => "CPU",
    "output_filename" => "spontaneous_tc_genesis_run_CPU",
    "f" => 3e-5,# coriolis parameter5e-4 #5e-4
    "g" => 9.8,#gravitational acceleration 9.8
    "convection_timescale" => 28800.0, #convective time scale
    "convection_critical_height"=> 130.0, #convection_triggering height
    "heating_amplitude"    => 1.0e9,#convective heating 1.0e9 #originally 9 for heating, -8 for cooling convective heating
    "large_scale_forcing" => (1.12/3)*1.0e-8,#radiative cooling rate
    "convective_radius"    => 30_000, #convective radius
    "relaxation_parameter" => 1.0/(2*86400), #relaxation timescale
    "relaxation_height" => nothing, #height to relax to, if nothing, h relax to its mean
    "Lx" => 2e6, #domain
    "Ly" => 2e6,
    "Lz" => 126.5,
    "Nx" => 250,
    "Ny" => 250,
    "boundary_layer" => true, #if true Lz=30 otherwise Lz = 41
    "initialization_style" => "rand",
    "initialization_amplitude" => 4.0,
    "simulation_length_in_days" => 86400.0*20,
    "output_interval_in_seconds" => 7200, 
    "timestep_in_seconds" => 40.0,
)
using RamirezReyes_ShallowWaterInFPlane

run_shallow_simulation(parameters)