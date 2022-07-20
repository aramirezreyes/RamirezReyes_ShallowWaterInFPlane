###
###This file creates dicts of parameters and launches slurm jobs with them to run an oceananigans simulation
###
using DrWatson
@quickactivate "RamirezReyes_ShallowWaterInFPlane"
#using RamirezReyes_ShallowWaterInFPlane

parameter_space = Dict(
    "architecture" => "GPU",
    "f" => [5e-3, 1e-4, 5e-5, 3e-5, 1e-5, 5e-6, 1e-6, 0.0],# coriolis parameter5e-4 #5e-4
    "g" => 10.0,#gravitational acceleration 9.8
    "convection_timescale" => [28800.0, 10800.0], #convective time scale
    "convection_critical_height"=> 130.0, #convection_triggering height
    "heating_amplitude"    => [5e11,1e11,5e10,1e10,5e9,1e9],#convective heating 1.0e9 #originally 9 for heating, -8 for cooling convective heating
    "large_scale_forcing" => (1.12/3.0)*[1.0e-8],#radiative cooling rate
    "convective_radius"    => [30_000.0,20_000.0,10_000.0, 5000.0], #convective radius
    "relaxation_parameter" => 1.0/(2*86400.0), #relaxation timescale
    "relaxation_height" => 131.0, #height to relax to
    "Lx" => 6e6, #domain
    "Ly" => 6e6,
    "Lz" => 126.0e0,
    "Nx" => 1500,
    "Ny" => 1500,
    "boundary_layer" => true,
    "initialization_style" => "rand",
    "initialization_amplitude" => 4.0,
    "simulation_length_in_days" => 100,
    "output_interval_in_seconds" => 7200.0, 
    "timestep_in_seconds" => 20.0,

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
