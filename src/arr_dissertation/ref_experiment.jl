###
###This file creates dicts of parameters and launches slurm jobs with them to run an oceananigans simulation
###

## In here we will setup experiments for varying τd while keeping Sc, and c constant

function create_reference_experiment()
    lz(bl, hc, ampl) = bl ? hc - ampl : hc + ampl 
    params = Dict(
    "architecture" => "GPU",
    "f" => 0.0,# coriolis parameter5e-4 #5e-4
    "g" => 10.0,#gravitational acceleration 9.8
    "convection_timescale" => 14400.0, #convective time scale
    "convection_critical_height"=> 40.0, #convection_triggering height
    "heating_amplitude"    => 5e8,#convective heating 1.0e9 #originally 9 for heating, -8 for cooling convective heating
    "large_scale_forcing" => Derived("heating_amplitude", x -> x/1.3392857142857142e14),#radiative cooling rate corresponding to pair q=5e7, r = 1.12/3*1.0e-6
    "convective_radius"    => 30_000.0, #convective radius
    "relaxation_parameter" => 1.0/(0.8*86400.0), #relaxation timescale
    "Lx" => 48e6, #domain
    "Ly" => 48e6,
    "Nx" => 6000,
    "Ny" => 6000,
    "boundary_layer" => true,
    "initialization_style" => "rand",
    "initialization_amplitude" => 0.01,
    "Lz" => Derived(["boundary_layer","convection_critical_height","initialization_amplitude"], lz),
    "simulation_length_in_days" => 200,
    "output_interval_in_seconds" => 43200.0, 
    "timestep_in_seconds" => 60.0,
    "relaxation_height" => nothing,
    "restart" => false,
    "checkpoint_interval_in_seconds" => 50*86400.0,
)
end
