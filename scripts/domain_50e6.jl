###
###This file creates dicts of parameters and launches slurm jobs with them to run an oceananigans simulation
###
using DrWatson
@quickactivate "RamirezReyes_ShallowWaterInFPlane"
#using RamirezReyes_ShallowWaterInFPlane

parameter_space = Dict(
    "architecture" => "GPU",
    "f" => 0.0,# coriolis parameter5e-4 #5e-4
    "g" => 10.0,#gravitational acceleration 9.8
    "convection_timescale" => [14400.0,7200.0], #convective time scale
    "convection_critical_height"=> 40.0, #convection_triggering height
    "heating_amplitude"    => [5e7, 5e6, 5e5],#convective heating 1.0e9 #originally 9 for heating, -8 for cooling convective heating
    "large_scale_forcing" => Derived("heating_amplitude", x -> x/1.3392857142857142e14),#radiative cooling rate corresponding to pair q=5e7, r = 1.12/3*1.0e-6
    "convective_radius"    => 30_000.0, #convective radius
    "relaxation_parameter" => 1.0/(0.8*86400.0), #relaxation timescale
    "Lx" => 48e6, #domain
    "Ly" => 48e6,
    "Nx" => 6000,
    "Ny" => 6000,
    "boundary_layer" => [true,false],
    "initialization_style" => "rand",
    "initialization_amplitude" => 4.0,
    "simulation_length_in_days" => 200,
    "output_interval_in_seconds" => 28800.0, 
    "timestep_in_seconds" => 60.0,
    "restart" => false,
    "checkpoint_interval_in_seconds" => 10*86400.0,

)


dicts = dict_list(parameter_space)

map(dicts) do dict
    dict["Lz"] = dict["boundary_layer"] ? dict["convection_critical_height"] - 3.9 :  dict["convection_critical_height"] + 3.9 
    dict["relaxation_height"] = nothing
end

res = tmpsave(dicts)
for r in Iterators.partition(res,4)
    arg1 = r[1]
    arg2 = r[2]
    arg3 = r[3]
    arg4 = r[4]
    submit = `sbatch submit_simulations_sbatch.bash $arg1 $arg2 $arg3 $arg4`
    #submit = `bash print_jobs.bash $arg1 $arg2 $arg3 $arg4`
    run(submit)
end
