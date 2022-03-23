###
###This file creates dicts of parameters and launches slurm jobs with them to run an oceananigans simulation
###
using DrWatson
@quickactivate "RamirezReyes_ShallowWaterInFPlane"
#using RamirezReyes_ShallowWaterInFPlane

parameter_space = Dict(
    "architecture" => "GPU",
    "coriolis_parameter" => 0.0,#5e-4 #5e-4
    "gravitational_acceleration" => 10,# 9.8
    "convection_timescale" => 28800,
    "critical_height"=> 40,
    "heating_amplitude"    => 5.0e11,#1.0e9 #originally 9 for heating, -8 for cooling
    "radiative_cooling_rate" => (1.12/3)*1.0e-8,
    "convective_radius"    => 30000.0,
    "relaxation_timescale" => 1.0/7200.0,
    "relaxation_height" => 38,
    "Lx" => 1.5e6,
    "Ly" => 1.5e6,
    "Lz" => 45,
    "Nx" => 512,
    "Ny" => 512
)


dicts = dict_list(parameter_space)
res = tmpsave(dicts)
for r in res
    println(r)
    submit = `sbatch submit_simulations_sbatch.bash $r`
    run(submit)
end
