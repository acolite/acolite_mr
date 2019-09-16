from .acolite_mr_ac import acolite_mr_ac

#import os
#path=os.path.dirname(os.path.abspath(__file__))
#for i in range(2): path = os.path.split(path)[0]

### check if binary distribution
#if '{}dist{}acolite'.format(os.path.sep,os.path.sep) in path:
#    ## two more levels for this file
#    for i in range(2): path = os.path.split(path)[0]

#cfile='{}{}config{}acolite_config.txt'.format(path,os.path.sep,os.path.sep)
#config = import_config(cfile)

### test whether we can find the relative paths
#for t in config:
#    if t in ['version']: continue
#    if os.path.exists(config[t]): continue
#    tmp = path + os.path.sep + config[t]
#    config[t] = os.path.abspath(tmp)
