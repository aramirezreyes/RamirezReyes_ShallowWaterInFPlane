###
### This file reads a parameter file and runs a simulation. It is meant to be called from submit_simulations_sbatch.bash or submit_simulations_pbs.bash
###
using DrWatson
@quickactivate "RamirezReyes_ShallowWaterInFPlane"
using RamirezReyes_ShallowWaterInFPlane

params_file = ARGS[1]
@info "I received:", ARGS
@info "Reading file:", params_file
params = load(projectdir("_research", "tmp", params_file), "params")
@info params
run_shallow_simulation(params)
