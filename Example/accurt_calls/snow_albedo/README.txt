
To get an overview of all input parameters that can be set to
configure Flick, run the Python script

  `make_config.py`

and then inspect the generated `flick_tmp/config` text file. In this
example, the config file is overwritten for each run, so it should
only be modified via the `set` functions in the Python scripts.

Run

  `compute_albedo.py`

and inspect the generated `output/computed_albedo.txt` file, which can
be plotted with

  `plot_albedo.py`

Modify these scripts as desired. Note that any parameter listed in
`flick_tmp/config` can be added to the Python scripts and set to the
desired value.




