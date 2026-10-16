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

c = flick.accurt_config()
c.set('detector_orientation','down')
c.set('detector_type','plane_irradiance')
c.set('detector_height', 2)
c.set('reference_detector_height', 2)
c.set('print_iops','true')
c.set('pure_water_volume_fraction',0)
