###############################################################################
import time, os
from math import sqrt
import subprocess, traceback

from write_run_data import (
	write_initial_data, write_run_data, write_slurm_script
)
from make_read_write_files_dirs import (
	make_directory_name_current_time_and_Nx
)
###############################################################################
### launch N sims with fixed initial data
### subprocess.call: wait uNtil subprocess over to coNtinue
###############################################################################
#def convergence_test(
#	sleep_time:float, run_data:dict, sim_number:int
#) -> None:
def convergence_test(
	sleep_time, run_data, sim_number
):
	Nx  = run_data["Nx"]
	Nt  = run_data["Nt"]
	tss = run_data["t_step_save"]

	for iC in range(0,sim_number): 
		run_data["output_dir"] = make_directory_name_current_time_and_Nx(run_data)
		print("\n",run_data["output_dir"],"\n")
		os.makedirs(run_data["output_dir"])

		write_initial_data(run_data)
		write_run_data(run_data)

		if (run_data["computer"]=="Feynman_cluster"):
			write_slurm_script(run_data)
			subprocess.call("sbatch run_TEdGB_collapse.slurm", shell="True")  
		else:
			subprocess.call(
				"./sim >{}/output_{}.txt 2>&1 &".format(
					run_data["output_dir"],run_data["theory"]),
			 shell="True")  

		time.sleep(sleep_time)

		run_data["Nx"] = str(2*(int(run_data["Nx"])-1)+1)
		run_data["t_step_save"] = str(2*(int(run_data["t_step_save"])))
		run_data["Nt"] = str(2*(int(run_data["Nt"])-1)+1)
###
###	return with the original Nx, etc. values
###	
	run_data["Nx"] = Nx
	run_data["Nt"] = Nt
	run_data["t_step_save"] = tss

	return
