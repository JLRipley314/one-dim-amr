##############################################################################
### format for writing to initial data and run params for
### initial data file, run param file, and slurm script
##############################################################################
def format_run_data(run_data: dict) -> dict:

	run_data["courant_n"]    = "{:.10e}".format(float(run_data["courant_n"]))
	run_data["errlim"]       = "{:.10e}".format(float(run_data["errlim"]))
	run_data["coupling_gbs"] = "{:.10e}".format(float(run_data["coupling_gbs"]))

	run_data["r2Exp_amp"]   = "{:.10e}".format(float(run_data["r2Exp_amp"]))
	run_data["r2Exp_width"] = "{:.10e}".format( float(run_data["r2Exp_width"]))
	run_data["r2Exp_r0"]    = "{:.10e}".format( float(run_data["r2Exp_r0"]))
	run_data["r2Exp_power"] = "{:.10e}".format( float(run_data["r2Exp_power"]))

	run_data["initial_black_hole_mass"] = "{:.10e}".format(float(run_data["initial_black_hole_mass"]))

	for param in run_data:
		run_data[param] = str(run_data[param])

	return run_data 
