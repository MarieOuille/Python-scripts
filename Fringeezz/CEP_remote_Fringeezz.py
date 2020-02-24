# -*- coding: utf-8 -*-
"""
Created on Fri Nov 22 15:48:34 2019

Fastlite's code for CEP remote control
@author: ouille
"""

import Pyro4
proxy = Pyro4.Proxy('PYRO:fringeezz@PCname:54310') # connect to Fringeezz, use the correct hostname instead of PCname
proxy.get_params_dict() # not required: check currently set parameters
# next part can be run in a loop
proxy.stop_acquisition() # important: stop acquisition before changing the target
proxy.set_params_dict(dict(target=1)) # set the target
proxy.start_acquisition() # restart acquisition
proxy.get_stats_dict() # not required: check the current CEP stats (mean, std, min, max)