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

c = flick.accurt_config()
c.set_streams(8)
c.set('detector_orientation','up')
c.set('detector_height', maximum_snow_thickness)
c.set('reference_detector_height', maximum_snow_thickness)
c.set('print_iops','true')
