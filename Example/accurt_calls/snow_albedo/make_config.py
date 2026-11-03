"""
Generate flick_tmp/config

See also README.txt
"""

import os
import flick
flick.setup_script(__file__)


c = flick.accurt_config()
c.set('detector_orientation','down')
c.set('detector_type','plane_irradiance')
c.set('detector_height', 2)
c.set('reference_detector_height', 2)
c.set('print_iops','true')
c.set('pure_water_volume_fraction',0)
