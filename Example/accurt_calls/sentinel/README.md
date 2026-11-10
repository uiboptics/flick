To get an overview of all input parameters available for configuring
Flick, run the Python script

`make_config.py`

and inspect the resulting `flick_tmp/config` text file. In this example,
the configuration file is overwritten on every run and should therefore
be modified only through the `set` functions in the Python scripts.

Next, run

`compute.py`

and inspect the resulting `output/` directory.

## Jupyter notebooks

The same workflow can be run interactively with the Jupyter notebooks in
this directory. Start Jupyter Lab or Jupyter Notebook from this directory so
that the notebooks can find the Python scripts and their generated files.

Run the notebooks in this order:

1. `run_make_config.ipynb` runs `make_config.py` and displays the generated
   `flick_tmp/config` file. The configuration is overwritten each time the
   notebook runs, so make configuration changes through the `set` functions
   in the Python scripts.
2. `run_compute.ipynb` runs `compute.py`, displays Matplotlib output in the
   notebook, and shows the saved optimization plot and best fitted parameter
   values. The computation invokes Flick-AccuRT several times and may take a
   while.

The computation notebook writes its results to `output/`, including
`optimization_result.png` and `best_parameter_values.txt`.



