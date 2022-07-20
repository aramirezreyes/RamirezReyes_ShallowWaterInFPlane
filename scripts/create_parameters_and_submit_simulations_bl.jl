###
###This file creates dicts of parameters and launches slurm jobs with them to run an oceananigans simulation
###
using DrWatson
@quickactivate "RamirezReyes_ShallowWaterInFPlane"
#using RamirezReyes_ShallowWaterInFPlane

parameter_space = Dict(
    "architecture" => "GPU",
    "f" => [5f-3, 1f-4, 5f-5, 3f-5, 1f-5, 5f-6],# coriolis parameter5e-4 #5e-4
    "g" => 10.0f0,#gravitational acceleration 9.8
    "convection_timescale" => [28800.0f0, 10800.0f0], #convective time scale
    "convection_critical_height"=> 130.0f0, #convection_triggering height
    "heating_amplitude"    => [5f11,1f11,5f10,1f10],#convective heating 1.0e9 #originally 9 for heating, -8 for cooling convective heating
    "large_scale_forcing" => (1.12f0/3.0f0)*[1.0f-8, 1.0f-9],#radiative cooling rate
    "convective_radius"    => [30000.0f0,20000.0f0,10000.0f0], #convective radius
    "relaxation_parameter" => 1.0/(2*86400.0), #relaxation timescale
    "relaxation_height" => 131.0f0, #height to relax to
    "Lx" => 2f6, #domain
    "Ly" => 2f6,
    "Lz" => 126.0f0,
    "Nx" => 500,
    "Ny" => 500,
    "boundary_layer" => true,
    "initialization_style" => "rand",
    "initialization_amplitude" => 4.0,
    "simulation_length_in_days" => 100,
    "output_interval_in_seconds" => 7200, 
    "timestep_in_seconds" => 5,
)


dicts = dict_list(parameter_space)
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
