from acolite import pleiades
from acolite import plotting
from acolite import ac
from acolite import aerlut
from acolite import dem

from acolite import acolite
from acolite import planetscope
from acolite import worldview

from .shared import *
from . import output 

#
#import os
#from acolite.shared import import_config
#pp_path = os.path.dirname(__file__)
#config=import_config(pp_path+'/../config/config.txt')

#
import os
path = os.path.dirname(os.path.abspath(__file__))

if not os.path.exists('{}{}config'.format(path, os.path.sep)):
    path = os.path.split(path)[0]

    ## check if binary distribution
    if '{}dist{}acolite'.format(os.path.sep,os.path.sep) in path:
        ## two levels for this file
        for i in range(2): path = os.path.split(path)[0]

cfile='{}{}config{}config.txt'.format(path,os.path.sep,os.path.sep)
config = import_config(cfile)

## test whether we can find the relative paths
for t in config:
    if os.path.exists(config[t]): continue
    tmp = path + os.path.sep + config[t]
    config[t] = os.path.abspath(tmp)
