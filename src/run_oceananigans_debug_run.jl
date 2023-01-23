"""
    This is intended to be launched from scripts/read_parameter_file_and_launch_15d_simulation.jl
    """
function run_shallow_simulation_debug(arch; ultrashort = false, restart = false)

    #Setup physicsl parameters
    parameters_dict = create_debug_parameters(arch; ultrashort, restart)

    run_shallow_simulation(parameters_dict)

end #runsimulation
