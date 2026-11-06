To get an overview of all input parameters available for configuring
Flick, run the Python script

`make_config.py`

and inspect the resulting `flick_tmp/config` text file. In this example,
the configuration file is overwritten on every run and should therefore
be modified only through the `set` functions in the Python scripts.

Next, run

`compute.py`

and inspect the resulting data in the `output/` directory. The data can
be plotted by running

`plot.py`

Modify these scripts as desired. Any parameter listed in
`flick_tmp/config` can be added to the Python scripts and set to the
desired value.

To avoid conflicts with changes introduced by future Flick releases,
it is recommended that you copy this example directory before making
any modifications.




