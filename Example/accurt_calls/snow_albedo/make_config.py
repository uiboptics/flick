"""
Generate flick_tmp/config

See also README.txt
"""

import os
import sys
sys.path.append(os.environ['FLICK_PATH']+'/python_script')
import flick

from pathlib import Path
os.chdir(Path(__file__).resolve().parent)

maximum_snow_thickness = 5

f = flick.relative_radiation()
f.generate_config("default", "flick_tmp")
f.set('detector_orientation','down')
f.set('detector_type','plane_irradiance')
f.set('detector_height', maximum_snow_thickness)
f.set('reference_detector_height', maximum_snow_thickness)
f.set('print_iops','true')
f.set('pure_water_volume_fraction',0)
