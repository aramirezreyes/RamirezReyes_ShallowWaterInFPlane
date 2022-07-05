###
###This file reads a parameter file and runs a simulation. It is meant to be called from submit_simulations_sbatch.bash
###
using DrWatson
@quickactivate "RamirezReyes_ShallowWaterInFPlane"
using RamirezReyes_ShallowWaterInFPlane

params_file = ARGS[1]
params = load(projectdir("_research", "tmp", params_file), "params")
@info params
run_shallow_simulation_100d(params)
