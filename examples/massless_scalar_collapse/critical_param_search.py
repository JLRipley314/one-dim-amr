##############################################################################
### search for critical param for black hole formation
### search output.out for outer_null > 0
### 	if >0 then decrease param
###		if neer ==0 for whole sim then increase param
##############################################################################
import subprocess, sys, os, time 

from write_run_data import (
	write_initial_data, write_run_data, 
	write_slurm_script
)
from format_run_data import format_run_data

from make_read_write_files_dirs import (
	make_new_file,
	append_to_file,
	make_directory_name_current_time
)

CRITICAL_PARAM_ERROR = 1e-6
OUTPUT_DIR = "/home/jripley/one-dim-amr/examples/massless_scalar_collapse/output"

###############################################################################
def compute_normalized_difference(lower:float,upper:float)->float:
	return (
		2*(upper-lower)/(upper+lower)
	)
###############################################################################
### if finds black hole in output then True, otherwise False
###############################################################################
def check_output_black_hole_formed(output_file: str) -> bool:
	with open(output_file, 'r') as f:
		for line in f:
			if (line.split(':')[0]  == "outermost trapped jC"):
				return True
	return False
###############################################################################
### if run finsihed then True, otherwise False
###############################################################################
def check_output_run_finished(output_file: str) -> bool:
	with open(output_file, 'r') as f:
		for line in f:
			if (line.split(':')[0] == "time evolving sim (sec)"):
				return True
	return False
###############################################################################
### if code crashed then True, otherwise False
###############################################################################
def check_output_code_crashed(output_file: str) -> bool:
	with open(output_file, 'r') as f:
		for line in f:
			if (line.split(':')[0] == 'ERROR'):
				return True
	return False
###############################################################################
###############################################################################
def critical_param_search(
	run_data: dict, param: str, param_range:list ### List[float]
) -> None:

	found_critical_param = False

	param_search_file_output = "{}/critical_param_search.txt".format(OUTPUT_DIR)
	make_new_file(param_search_file_output)

	while found_critical_param == False:

		run_data[param] = (param_range[0] + param_range[1]) / 2
		run_data["output_dir"] = make_directory_name_current_time(run_data)
		os.makedirs(run_data["output_dir"])
		output_file = run_data['output_dir'] + 'output.txt'

		normalized_difference = compute_normalized_difference(param_range[0],param_range[1])

		append_to_file(param_search_file_output, param + ' value: ' + str(run_data[param]) + '\n') 
		append_to_file(param_search_file_output, 'range: [{},{}] \n'.format(param_range[0], param_range[1])) 
		append_to_file(param_search_file_output, 'difference: {:.2e} \n'.format(normalized_difference)) 

		write_initial_data(run_data)
		write_run_data(run_data)
		this_run = subprocess.call(
			"./sim >{}/output.txt 2>&1 &".format(
				run_data["output_dir"],run_data["theory"]),
		 shell="True")  

		black_hole_formed = False 
		code_crashed = False 
		run_finished = False 
		while True:
			time.sleep(10) ### check every n seconds
			black_hole_formed = check_output_black_hole_formed(output_file)
			code_crashed = check_output_code_crashed(output_file)
			run_finished = check_output_run_finished(output_file)

			if ((black_hole_formed == True)
			or  (code_crashed == True)
			or  (run_finished == True)
			):
				break

		if (black_hole_formed == True):
			param_range[1] = (param_range[0] + param_range[1]) / 2
			append_to_file(param_search_file_output, 'black_hole_formed\n') 

		elif (code_crashed == True):
			param_range[1] = (param_range[0] + param_range[1]) / 2
			append_to_file(param_search_file_output, 'code_crashed\n') 

		elif (run_finished == True):
			param_range[0] = (param_range[0] + param_range[1]) / 2
			append_to_file(param_search_file_output, 'run_finished\n') 

		else:
			raise ValueError('black_hole_formed, code_crashed, run_finished all are False!')

		if (normalized_difference < CRITICAL_PARAM_ERROR):
			found_critical_param = True
		else:
			continue 

	return
