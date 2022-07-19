###
###This file creates dicts of parameters and launches slurm jobs with them to run an oceananigans simulation
###
using DrWatson
@quickactivate "RamirezReyes_ShallowWaterInFPlane"
#using RamirezReyes_ShallowWaterInFPlane

parameter_space = Dict(
    "architecture" => "GPU",
    "f" => [0.0f0, 1f-3, 5f-4, 1f-4, 5f-5],# coriolis parameter5e-4 #5e-4
    "g" => 10.0f0,#gravitational acceleration 9.8
    "convection_timescale" => [28800.0f0, 10800.0f0], #convective time scale
    "convection_critical_height"=> 40.0f0, #convection_triggering height
    "heating_amplitude"    => [5f10,1f10],#convective heating 1.0e9 #originally 9 for heating, -8 for cooling convective heating
    "large_scale_forcing" => (1.12f0/3.0f0)*[1.0f-8, 1.0f-9],#radiative cooling rate
    "convective_radius"    => [30000.0f0,20000.0f0,10000.0f0], #convective radius
    "relaxation_parameter" => 1.0/(2*86400.0), #relaxation timescale
    "relaxation_height" => nothing, #height to relax to, if nothing, h relax to its mean
    "Lx" => 6f6, #domain
    "Ly" => 6f6,
    "Lz" => 45.0f0,
    "Nx" => 2000,
    "Ny" => 2000,
    "boundary_layer" => false,
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
