# This file creates the simulations that vary the Coriolis parameter in my doctoral dissertation
using DrWatson
@quickactivate "RamirezReyes_ShallowWaterInFPlane"                                                                                                   

# 1. We pull the parameters of the reference simulation
# 2. We substitute the Coriolis parameter by a list of coriolis parameters
# 3. We create parameters for all each combinations
# 4. We launch jobs to do simulations

parameter_space = create_reference_experiment()
parameter_space["f"] = [1e-6, 2e-6, 4e-6, 8e-6, 1e-5, 2e-5, 4e-5, 8e-5]

# dict_list is a DrWatson function. It will create a dict for each combination of parameters (in this case for every Coriolis value)
dicts = dict_list(parameter_space)
                                                                                                                                       
# This will submit jobs, each job expects to run 4 jobs on 4 different GPUS
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