#############################################################################

import shutil

#############################################################################
def write_initial_data(run_data: dict) -> None:

	with open("{}/initial_data.txt".format(run_data["home_dir"]), "w") as f:
		f.write("amp={}\n".format(run_data["amp"]))
		f.write("width={}\n".format(run_data["width"]))
		f.write("center={}\n".format(run_data["center"]))
		f.write("direction={}\n".format(run_data["direction"]))
		f.write("initial_data={}\n".format(run_data["initial_data"]))

		f.write("initial_black_hole_mass={}\n".format(run_data["initial_black_hole_mass"]))

	shutil.copyfile(
		"{}/initial_data.txt".format(run_data["home_dir"]),
		"{}/initial_data.txt".format(run_data["output_dir"])
	)

	return
#############################################################################
def write_run_data(run_data: dict) -> str:

	with open("{}/run_data.txt".format(run_data["home_dir"]), "w") as f:
		f.write("output_dir={}\n".format(run_data["output_dir"]))
		f.write("param_search={}\n".format(run_data["param_search"]))
		f.write("theory={}\n".format(run_data["theory"]))
		f.write("solver_Ze={}\n".format(run_data["solver_Ze"]))
		f.write("use_excision={}\n".format(run_data["use_excision"]))
		f.write("Nx={}\n".format(run_data["Nx"]))
		f.write("Nt={}\n".format(run_data["Nt"]))
		f.write("t_step_save={}\n".format(run_data["t_step_save"]))
		f.write("stereographic_L={}\n".format(run_data["stereographic_L"]))
		f.write("cfl_num={}\n".format(run_data["cfl_num"]))
		f.write("coupling_gbs={}\n".format(run_data["coupling_gbs"]))
		f.write("err_tolerance={}\n".format(run_data["err_tolerance"]))

	shutil.copyfile(
		"{}/run_data.txt".format(run_data["home_dir"]),
		"{}/run_data.txt".format(run_data["output_dir"])
	)
	
	return
#############################################################################
def write_slurm_script(run_data: dict) -> None:

	run_data["var_name"] = "output"

	outputName = "{}/{}".format(run_data["output_dir"],"output.out") 	

	with open("{}/run_TEdGB_collapse.slurm".format(run_data["home_dir"]), "w") as f:
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
		"{}/run_TEdGB_collapse.slurm".format(run_data["home_dir"]),
		"{}/run_TEdGB_collapse.slurm".format(run_data["output_dir"])
	)		
	return 
