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

CRITICAL_PARAM_SEARCH_OUTPUT = (
	'/home/jripley/EdGB_HorizonPenetrating/critical_param_search_output/'
)

###############################################################################
### if finds black hole in output then True, otherwise False
###############################################################################
def check_output_black_hole_formed(output_file: str) -> bool:
	with open(output_file, 'r') as f:
		for line in f:
			line = line.split('\t')
			for component in line:
				component = component.split(':')
				if ((str(component[0]) == 'outer_null')
				and (int(component[1]) > 0)
				):
					return True
	return False
###############################################################################
### if run finsihed then True, otherwise False
###############################################################################
def check_output_run_finished(output_file: str) -> bool:
	with open(output_file, 'r') as f:
		for line in f:
			if line.split(':')[0] == 'Total sim time (sec)':
				return True
	return False
###############################################################################
### if code crashed then True, otherwise False
###############################################################################
def check_output_code_crashed(output_file: str) -> bool:
	with open(output_file, 'r') as f:
		for line in f:
			if ((line.split(':')[0] == 'ERROR')
			or  (line.split(':')[0] == 'jC')
			):
				return True
	return False
###############################################################################
###############################################################################
def critical_param_search(
	run_data: dict, param: str, param_range:list ### List[float]
) -> None:

	critical_param_error = 1e-8
	found_critical_param = False

	param_search_file_output = CRITICAL_PARAM_SEARCH_OUTPUT + 'critical_value_{}_{}_{}.txt'.format(param,run_data['Ny'],run_data['coupling_gbs'])
	make_new_file(param_search_file_output)

	while found_critical_param == False:

		run_data[param] = (param_range[0] + param_range[1]) / 2

		append_to_file(param_search_file_output, param + ' value: ' + str(run_data[param]) + '\n') 
		append_to_file(param_search_file_output, 'range: [{},{}] \n'.format(param_range[0], param_range[1])) 
		append_to_file(param_search_file_output, 'difference: {:.2e} \n'.format(param_range[1]-param_range[0])) 

		run_data['output_dir'] = make_directory_name_current_time(run_data['output_dir'])
		os.makedirs(run_data["output_dir"])

		output_file = run_data['output_dir'] + 'output.out'

		write_txtFiles_initialData(run_data)
		write_txtFiles_basicRunParams(run_data)
		write_slurmScript(run_data)
		subprocess.call("sbatch run_TEdGB_collapse.slurm", shell="True")  

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

		if 2*(param_range[1] - param_range[0])/(param_range[1] + param_range[0]) < critical_param_error:
			found_critical_param = True
		else:
			continue 

	return
