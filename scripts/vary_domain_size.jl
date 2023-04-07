
using DrWatson
@quickactivate "RamirezReyes_ShallowWaterInFPlane"
#using RamirezReyes_ShallowWaterInFPlane                                                                          

lz(bl, hc, ampl) = bl ? hc - ampl : hc + ampl

reference_experiment = Dict(
    "architecture" => "GPU",
    "f" => 0.0,# coriolis parameter5e-4 #5e-4                                                                     
    "g" => 10.0,#gravitational acceleration 9.8                                                                   
    "convection_timescale" => 14400.0, #convective time scale                                                     
    "convection_critical_height"=> 40.0, #convection_triggering height                                            
    "heating_amplitude"    => 5e7,#convective heating 1.0e9 #originally 9 for heating, -8 for cooling convective \heating                                                                                                           
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
    "simulation_length_in_days" => 200,
    "output_interval_in_seconds" => 42300.0,
    "timestep_in_seconds" => 60.0,
    "relaxation_height" => nothing,
    "restart" => false,
    "checkpoint_interval_in_seconds" => 50*86400.0,
)


parameter_space = copy(reference_experiment)

Δx = parameter_space["Lx"]/parameter_space["Nx"] #We will change the domain size but maintaining this constant.   

parameter_space["Lx"] = parameter_space["Lx"] .÷ [2, 4, 8, 16, 24, 48, 96, 100]
parameter_space["Ly"] = Derived("Lx", identity)
parameter_space["Nx"] = Derived("Lx", x -> floor(Int,x / Δx))
parameter_space["Ny"] = Derived("Nx", identity)


dicts = dict_list(parameter_space)
#for dict in dicts                                                                                                
#    @info dict                                                                                                   
#end                                                                                                              

res = tmpsave(dicts)
for r in Iterators.partition(res,4)
   arg1 = r[1]
   arg2 = r[2]
   arg3 = r[3]
   arg4 = r[4]
   #submit = `sbatch submit_simulations_sbatch.bash $arg1 $arg2 $arg3 $arg4`
   submit =`qsub -v f1=$arg1,f2=$arg2,f3=$arg3,f4=$arg4 submit_simulations_pbs.bash`
   #submit = `bash print_jobs.bash $arg1 $arg2 $arg3 $arg4`                                                       
   run(submit)
end


