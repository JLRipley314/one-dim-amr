#############################################################################
### to write initial data and run param files. Depending on
### where you are running this (e.g. the Feynman cluster)
### also sets up slurm scripts, etc.  
#############################################################################

import subprocess, sys, time, os
from write_run_data import (
	write_initial_data, write_run_data, write_slurm_script
)
from format_run_data import format_run_data

from make_read_write_files_dirs import (
	make_directory_name_current_time
)

from convergence_test import convergence_test 

from critical_param_search import critical_param_search

##############################################################################
assert len(sys.argv) == 2, (
	"argv[1] is empty-meed a run_type to run!"
	)
	
run_type = sys.argv[1]
##############################################################################
### Initialize run data ##############################################################################
run_data = {
	"computer"	: "jlr_laptop",#"Feynman_cluster",#

	"param_search"	: "no",#"yes",#

	"theory"	: "EdGB",#"massless_scalar_GR",#
###
###	EdGB: coupling for Gauss-Bonnet scalar (gbs)
###
	"coupling_gbs"	: 5.0,

	"stereographic_L" : 100,
	
	"use_excision"	: "yes",#"no",#

	"solver_QP"	: "", 
	"solver_Al"	: "",
	"solver_Ze"	: "constrained",#"free",# 
	"solver_Ps"	: "",

	"characteristics_calculator" : "",
###
###	Nx, Nt should be of the form 2**n + 1 with n an integer
###	cfl_num: Courant Friedrichs Lewy number (for hyperbolics)
###
	"Nx"		: 2**7+1,
	"Nt"		: 2**17+1,#2**13+1,
	"t_step_save"	: 2**5,
	"cfl_num"	: 0.25,  
	"err_tolerance"	: 1.0e-10, 
	
	"initial_data" : "initial_black_hole",#"r4Exp",#"initial_black_hole_with_r4Exp",#"presetPQPhi",# 
###
###	if initial_data is r4Exp
###
	"amp"		: 0.029983671875,
	"width"		: 10.0,
	"center"	: 10.0,
	"direction"	: "ingoing",#"stationary",#"outgoing",#
###
###	if initial_data is initial_black_hole
###
	"initial_black_hole_mass" : 20.0,

	"var_name"	: "output",
###
### Information for the slurm scripts if running of Feynman
###
### walltime format: (dd:hh:mm:ss)
###
	"walltime"        : "8:00:00",
	"memory_usage_MB" : "2000"
}
##############################################################################
### output directory for the computer we are running on
##############################################################################
if (run_data["computer"]=="jlr_laptop"):
	run_data["base_output_dir"] = "/home/jripley/one-dim-amr/examples/massless_scalar_collapse/output"
elif (run_data["computer"]=="feynman_cluster"):
	run_data["base_output_dir"] = "/group/grtheory/EdGB_SSCollapse_output"
else:
	raise SystemError("'computer' value not set") 

run_data["home_dir"] = os.getcwd()
##############################################################################
### derived paramters
##############################################################################
dx = float(run_data["stereographic_L"]) / (float(run_data["Nx"])-1.0)
dt = float(float(run_data["cfl_num"]) * dx) 

run_data["dx"] = dx
run_data["dt"] = dt
##############################################################################
### machine eps: machine epsilon
##############################################################################
from numpy import finfo
run_data["machine_eps"] = finfo(float).eps
##############################################################################
### correct formatting for printing
##############################################################################
run_data = format_run_data(run_data)
##############################################################################
##############################################################################
if run_type == "basic_run":

	run_data["output_dir"] = make_directory_name_current_time(run_data)
	os.makedirs(run_data["output_dir"])

	write_initial_data(run_data)
	write_run_data(run_data)

	if (run_data["computer"]=="Feynman_cluster"):
		write_slurm_script(run_data)
		subprocess.call("sbatch run_TEdGB_collapse.slurm", shell="True")  
	else:
		subprocess.call(
			"./sim >{}/output.txt 2>&1 &".format(
				run_data["output_dir"],run_data["theory"]),
		 shell="True")  
##############################################################################
### with fixed initial data we run sim with N different resolutions:
### Nx, (2*(Nx-1)+1) , etc. 
### Note we assume Nx of the for 2**n + 1 with n integer
### sleep time: how long to wait before launching the next sim (in seconds)
##############################################################################
elif (run_type == "convergence_test"): 
	sim_number = int(input("Launch how many sims with different resolutions? (enter integer) "))
	sleep_time = 5  
	convergence_test(sleep_time, run_data, sim_number)
##############################################################################
### specify the param and range to search for edge of black hole formation
##############################################################################
elif (run_type) == "critical_param_search":

	param = "amp"	
	param_range = [0.02998,0.03]	
	critical_param_search(run_data, param, param_range)
##############################################################################
else:
	print("ERROR: run_type = {}".format(run_type))
