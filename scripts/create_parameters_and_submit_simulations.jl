###
###This file creates dicts of parameters and launches slurm jobs with them to run an oceananigans simulation
###
using DrWatson
@quickactivate "RamirezReyes_ShallowWaterInFPlane"
#using RamirezReyes_ShallowWaterInFPlane

parameter_space = Dict(
    "architecture" => "GPU",
    "f" => 0.0,# coriolis parameter5e-4 #5e-4
    "g" => 10,#gravitational acceleration 9.8
    "convec_t" => (28800, 10800), #convective time scale
    "h_c"=> 40, #convection_triggering height
    "q0"    => (5.0e11,1e11,5e10,1e10),#convective heating 1.0e9 #originally 9 for heating, -8 for cooling convective heating
    "r" => (1.12/3)*1.0e-8,#radiative cooling rate
    "convec_r"    => (30000,20000,10000), #convective radius
    "relax_t" => 1.0/(2*86400.0), #relaxation timescale
    "relax_h" => 39, #height to relax to
    "Lx" => 1.5e6, #domain
    "Ly" => 1.5e6,
    "Lz" => 45,
    "Nx" => 200,
    "Ny" => 200
)


dicts = dict_list(parameter_space)
res = tmpsave(dicts)
for r in res
    println(r)
    submit = `sbatch submit_simulations_sbatch.bash $r`
    run(submit)
end
