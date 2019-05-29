#############################################################################
###
### This formats the heredoc to feed into ssc_polarSlicing_main_run
###
#############################################################################

import shutil

#############################################################################
def write_initial_data(run_data: dict) -> None:

	with open("{}initial_data.txt".format(run_data["home_dir"]), "w") as f:
		f.write("r2Exp_amp={}\n".format(run_data["r2Exp_amp"]))
		f.write("r2Exp_width={}\n".format(run_data["r2Exp_width"]))
		f.write("r2Exp_r0={}\n".format(run_data["r2Exp_r0"]))
		f.write("r2Exp_power={}\n".format(run_data["r2Exp_power"]))
		f.write("r2Exp_direction={}\n".format(run_data["r2Exp_direction"]))

		f.write("initial_black_hole_mass={}\n".format(run_data["initial_black_hole_mass"]))

	shutil.copyfile(
		"{}initial_data.txt".format(run_data["home_dir"]),
		"{}initial_data.txt".format(run_data["output_dir"])
	)

	return
#############################################################################
def write_run_data(run_data: dict) -> str:

	with open("{}basic_run_params.txt".format(run_data["home_dir"]), "w") as f:
		f.write("output_dir={}\n".format(run_data["output_dir"]))
		f.write("param_search={}\n".format(run_data["param_search"]))
		f.write("theory={}\n".format(run_data["theory"]))
		f.write("solver_Ze={}\n".format(run_data["solver_Ze"]))
		f.write("use_excision={}\n".format(run_data["use_excision"]))
		f.write("initial_data={}\n".format(run_data["initial_data"]))
		f.write("Nx={}\n".format(run_data["Nx"]))
		f.write("Nt={}\n".format(run_data["Nt"]))
		f.write("t_step_save={}\n".format(run_data["t_step_save"]))
		f.write("stereographic_L={}\n".format(run_data["stereographic_L"]))
		f.write("courant_n={}\n".format(run_data["courant_n"]))
		f.write("coupling_gbs={}\n".format(run_data["coupling_gbs"]))
		f.write("errlim={}\n".format(run_data["errlim"]))

	shutil.copyfile(
		"{}basic_run_params.txt".format(run_data["home_dir"]),
		"{}basic_run_params.txt".format(run_data["output_dir"])
	)
	
	return
#############################################################################
def write_slurm_script(run_data: dict) -> None:

	run_data["varName"] = "output"

	outputName = "{}{}".format(run_data["output_dir"],"output.out") 	

	with open("{}run_TEdGB_collapse.slurm".format(run_data["home_dir"]), "w") as f:
		f.write("#!/bin/sh\n")
		f.write("#SBATCH -N 1\t\t# nodes=1\n")
		f.write("#SBATCH --ntasks-per-node=1\t\t# ppn=1\n")
		f.write("#SBATCH -J g{:.2f}\t\t# job name\n".format(float(run_data["coupling_gbs"])))
		f.write("#SBATCH -t {}\t\t# walltime (dd:hh:mm:ss)\n".format(run_data["walltime"]))
		f.write("#SBATCH -p dept\t\t# partition/queue name\n")
		f.write("#SBATCH --mem={}MB\t\t# memory in MB\n".format(run_data["memory_usage_MB"]))
		f.write("#SBATCH --output={}\t\t# file for STDOUT\n".format(outputName))
		f.write("#SBATCH --mail-user=jripley@princeton.edu\t\t# Mail  id of the user\n")
#		f.write("#SBATCH --mail-type=begin\t\t# Slurm will send mail at the beginning of the job\n")
#		f.write("#SBATCH --mail-type=end\t\t# Slurm will send at the completion of your job\n")
		f.write("\n./collapse\n\n")

	shutil.copyfile(
		"{}run_TEdGB_collapse.slurm".format(run_data["home_dir"]),
		"{}run_TEdGB_collapse.slurm".format(run_data["output_dir"])
	)		
	return 
