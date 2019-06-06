##############################################################################
### format for writing to initial data and run params for
### initial data file, run param file, and slurm script
##############################################################################
def format_run_data(run_data: dict) -> dict:

	run_data["cfl_num"] = "{:.10e}".format(float(run_data["cfl_num"]))
	run_data["err_tolerance"] = "{:.10e}".format(float(run_data["err_tolerance"]))
	run_data["coupling_gbs"]= "{:.10e}".format(float(run_data["coupling_gbs"]))

	run_data["amp"] = "{:.10e}".format(float(run_data["amp"]))
	run_data["width"] = "{:.10e}".format( float(run_data["width"]))
	run_data["center"] = "{:.10e}".format( float(run_data["center"]))

	run_data["initial_black_hole_mass"] = "{:.10e}".format(float(run_data["initial_black_hole_mass"]))

	for param in run_data:
		run_data[param] = str(run_data[param])

	return run_data 
