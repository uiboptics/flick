"""
Generate flick_tmp/config

See also README.txt
"""

import os
import flick
flick.setup_script(__file__)


c = flick.accurt_config()
c.set('detector_orientation','up')
c.set('detector_type','plane_irradiance')
c.set('reference_detector_height', flick.top_of_atmosphere)
c.set('print_iops','true')
