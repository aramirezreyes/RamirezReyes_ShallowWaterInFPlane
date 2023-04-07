# This is an example script to launch a simulation in a computer that uses the PBS manager. 
#If you use it you should also take a look at submit_simulations_pbs.bash
using DrWatson
@quickactivate "RamirezReyes_ShallowWaterInFPlane"                                                                                                   

lz(bl, hc, ampl) = bl ? hc - ampl : hc + ampl

# I am creating a dict of parameters. 
# I use Derived when I want a parameter to be a function of another one (in this case Lz depends on wether we use the boundary layer formulation or not)
# I use an array when I want to explore different parameters. In this case I am only varying the Coriolis parameter

parameter_space = Dict(
    "architecture" => "GPU",
    "f" => [1e-6, 2e-6, 4e-6, 8e-6, 1e-5, 2e-5, 4e-5, 8e-5],# coriolis parameter5e-4 #5e-4                                                                                                    
    "g" => 10.0,#gravitational acceleration 9.8                                                                                                  
    "convection_timescale" => 14400.0, #convective time scale                                                                                    
    "convection_critical_height"=> 40.0, #convection_triggering height                                                                           
    "heating_amplitude"    => 5e8,#convective heating 1.0e9 #originally 9 for heating, -8 for cooling convective \heating                                                                                                                                                                   
    "large_scale_forcing" => Derived("heating_amplitude", x -> x/1.3392857142857142e14),#radiative cooling rate cresponding to pair q=5e7, r = 1.12/3*1.0e-6                                                                                                                                     
    "convective_radius"    => 30_000.0, #convective radius                                                                                       
    "relaxation_parameter" => 1.0/(0.8*86400.0), #relaxation timescale                                                                           
    "Lx" => 48e6, #domain                                                                                                                        
    "Ly" => 48e6,
    "Nx" => 6000,
    "Ny" => 6000,
    "boundary_layer" => true,
    "initialization_style" => "rand",
    "initialization_amplitude" => 0.01,
    "Lz" => Derived(["boundary_layer","convection_critical_height","initialization_amplitude"], lz),
    "simulation_length_in_days" => 10,
    "output_interval_in_seconds" => 42300.0,
    "timestep_in_seconds" => 120.0,
    "relaxation_height" => nothing,
    "restart" => false,
    "checkpoint_interval_in_seconds" => 200*86400.0,
)

# dict_list is a DrWatson function. It will create a dict for each combination of parameters (in this case for every Coriolis value)
dicts = dict_list(parameter_space)
                                                                                                                                       
# This will submit jobs, each job expects to run 4 jobs on 4 different GPUS
res = tmpsave(dicts)
for r in Iterators.partition(res,4)
   arg1 = r[1]
   arg2 = r[2]
   arg3 = r[3]
   arg4 = r[4]
   submit = `sbatch submit_simulations_sbatch.bash $arg1 $arg2 $arg3 $arg4`                                                                     
   #submit =`qsub -v f1=$arg1,f2=$arg2,f3=$arg3,f4=$arg4 submit_simulations_pbs.bash`
   #submit = `bash print_jobs.bash $arg1 $arg2 $arg3 $arg4`                                                                                      
   run(submit)
end

