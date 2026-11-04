"""
Generate flick_tmp/config

See also README.txt
"""

import flick
flick.setup_script(__file__)

c = flick.accurt_config()
c.set('print_iops','true')
