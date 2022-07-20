"""
    This is intended to be launched from scripts/read_parameter_file_and_launch_15d_simulation.jl
    """
function run_shallow_simulation_debug(arch; ultrashort = false)

    #Setup physicsl parameters
    parameters_dict = Dict{String, Union{Nothing, Real, String}}(
    "architecture" => arch,
    "output_filename" => "debug_run_"*arch,
    "f" => 5e-4, #Coriolis
    "g" => 9.8, #Gravity
    "convection_timescale" => 10800.0, #Duration of convective events (in seconds)
    "convection_critical_height" => 130.0, #Critical height that triggers convection (in meters)
    "heating_amplitude" => 3e9, #The amplitude of the convective event (q0 in paper)
    "large_scale_forcing" => (1.12/3)*1.0e-8, #The amplitude of the large scale forcing
    "convective_radius"    => 30000.0, #Radius of convective event (in meters)
    "relaxation_parameter" => 1.0/(2*86400), # 1/τ where tau is the relaxation timescale for friction and h recovery
    "relaxation_height" => nothing, #131.0 #Target for the recovery of the h field if nothing it will relax to the mean
    "Lx" => 1.5e6, #Size of the domain (in meters)
    "Ly" => 1.5e6,
    "Lz" => 126.5, # A characteristic height
    "Nx" => 500, #Number of points
    "Ny" => 500,
    "boundary_layer" => true, #If true, convection is a mass sink, otherwise is false,
    "initialization_style" => "rand",
    "initialization_amplitude" => 4.0,
    "simulation_length_in_days" => ultrashort ? 100 / 86400 : 10000/86400,
    "output_interval_in_seconds" => 5, 
    "timestep_in_seconds" => 5,
    )
    
    run_shallow_simulation(parameters_dict)

end #runsimulation
