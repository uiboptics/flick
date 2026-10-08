To get an overview of all input parameters that can be set to
configure Flick, run the Python script

  `make_config.py`

and inspect the resulting `flick_tmp/config` text file. In this
example, the config file is overwritten on each run and should
therefore only be modified via the `set` functions in the Python
scripts.

Run

  `compute_albedo.py`

and inspect the resulting `output/computed_albedo.txt` file, which can
be plotted using

  `plot_albedo.py`.

Modify these scripts as desired. Note that any parameter listed in
`flick_tmp/config` can be added to the Python scripts and set to the
desired value. To avoid code conflicts with future Flick releases, it
is recommended to make a copy of this directory before making any
modifications.





