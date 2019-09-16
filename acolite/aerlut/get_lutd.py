## def get_lutd
## imports LUTs into dictionary
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2019-04-30
## modifications:

def get_lutd(luts = ['PONDER-LUT-201704-MOD1-1013mb', 'PONDER-LUT-201704-MOD2-1013mb'],
             pressures = ['0500mb', '1100mb']):
    import acolite as ac
    
    ## read luts
    lutdir= '{}/LUT/'.format(ac.config['pp_data_dir'])
    lutd={}
    for lutid in luts:
        ld, lm = ac.aerlut.import_lut(lutid,lutdir, override=0)
        lutd[lutid] = {'lut':ld, 'meta':lm}
        ## read luts at other pressures
        lut_split = lutid.split('-')
        for p in pressures:
            lut0 = '-'.join(lut_split[0:-1]+[p])
            lut_data, lut_meta = ac.aerlut.import_lut(lut0,lutdir, override=0)
            lutd[lut0] = {'lut':lut_data, 'meta':lut_meta}            
    return(lutd)
