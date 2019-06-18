###############################################################################
###
### appending to, reading, writing, appending, making, etc files and directories 
###
###############################################################################

import sys, os, fnmatch, time

###############################################################################
def make_directory_name_current_time(run_data:dict) -> str:

	dir_name = "_".join("_".join(time.asctime().split(" ")).split(":"))

	dir_name_path = "{}/{}/".format(run_data["base_output_dir"],dir_name)

	return dir_name_path
###############################################################################
def make_directory_name_current_time_and_Nx(run_data:dict) -> str:

	dir_name = "_".join("_".join(time.asctime().split(" ")).split(":"))

	dir_name += "_Nx{}".format(run_data["Nx"])
	dir_name += "_Nt{}".format(run_data["Nt"])

	dir_name_path = "{}/{}/".format(run_data["base_output_dir"],dir_name)

	return dir_name_path
###############################################################################
def make_new_file(file_name: str) -> None:
	open(file_name, "w").close()
	return 
###############################################################################
def append_to_file(file_name: str, message: str) -> None:
	with open(file_name, "a") as f:
		f.write(message) 
	return 
