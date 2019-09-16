## def select_model2
## new model selection based on generic LUTs
## written by Quinten Vanhellemont, RBINS for the PONDER project
## QV 2019-04-29
## modifications: QV 2019-04-30 changed rhod to dict with raa, vza, sza, wave and tt_gas per band
##                QV 2019-06-12 added rhod_model_selection min_tau

def select_model2(rhod, sensor,
                  pressure = None, 
                  rhod_list_selection = 'intercept', 
                  rhod_fit_bands = 2, rhod_fit_selection = 'min_drmsd',
                  rhod_model_selection = 'min_drmsd', 
                  rhod_tgas_cutoff = 0.95, rhod_min_wave = 400.,
                  lowess_frac = 0.5, 
                  luts = ['PONDER-LUT-201704-MOD1-1013mb', 'PONDER-LUT-201704-MOD2-1013mb'], 
                  lutd=None):
    
    import numpy as np
    import scipy.interpolate
    from statsmodels.nonparametric.smoothers_lowess import lowess
    import acolite as ac
    
    ## read luts
    lutdir= '{}/LUT/'.format(ac.config['pp_data_dir'])
    if lutd is None: lutd=ac.aerlut.get_lutd()
    
    ## get sensor RSR
    pp_path = ac.config['pp_data_dir']
    rsr_file = pp_path+'/RSR/'+sensor+'.txt'
    rsr, rsr_bands = ac.shared.rsr_read(file=rsr_file)

    ## select rhod
    rhod_sel = {}
    for b in rhod:
        if rhod[b]['rhod'].shape != (): ## if multiple rhod are given
            if rhod_list_selection == "smooth": ## use lowess smoothed min
                tmp = lowess(rhod[b]['rhod'],np.arange(0,len(rhod[b]['rhod'])),frac=lowess_frac)[:,1]
                rhod_sel[b] = tmp[0]
            elif rhod_list_selection == "list_smooth":  ## use lowess smoothed
                tmp = lowess(rhod[b]['rhod'],np.arange(0,len(rhod[b]['rhod'])),frac=lowess_frac)[:,1]
                rhod_sel[b] = tmp
            elif rhod_list_selection == "intercept": ##  use OLS intercept
                reg = ac.shared.regression.lsqfity(np.arange(0,len(rhod[b]['rhod'])), rhod[b]['rhod']) 
                rhod_sel[b] = reg[1]
            elif rhod_list_selection == "darkest": ## absolute darkest
                rhod_sel[b] = np.nanmin(rhod[b]['rhod'])
                
        ## other cases just copy the given rhod
        if b not in rhod_sel: rhod_sel[b] = rhod[b]['rhod']
            
    ## run through luts
    model_results = {}
    for lutid in luts: 
        if pressure is not None: ## interpolate LUTs to given pressure
            lut_data, lut_meta = ac.aerlut.aerlut_pressure_hyper(lutid, lutdir, pressure, lut_data_dict=lutd)
            lutidp = lutid.replace('1013', str(pressure))
            lutd[lutidp] = {'lut':lut_data, 'meta':lut_meta}
        else:
            lut_data, lut_meta = lutd[lutid]['lut'], lutd[lutid]['meta']

        ## set up LUT dimensions for interpolation
        # lut data is par, wave, azi, thv, ths, wind, tau
        lut_dim = [[i for i in range(len(lut_meta['par']))]]
        lut_dim += [lut_meta[k] for k in ['wave', 'azi', 'thv', 'ths', 'tau']]

        ## set up LUT interpolation
        rgi = scipy.interpolate.RegularGridInterpolator(lut_dim, lut_data[:,:,:,:,:,0,:], 
                                                      bounds_error=False, fill_value=np.nan)

        taua_steps = lut_meta['tau']
        waves_mu = lut_meta['wave']
        ip = [i for i,value in enumerate(lut_meta['par']) if value == 'romix']

        ## get romix for each band corresponding to taua steps
        ratm_bands = {b:[] for b in rhod}
        for ta in taua_steps:
            #ret = rgi((ip, waves_mu, raa, vza, sza, ta))
            for b in rhod:
                ret = rgi((ip, waves_mu, rhod[b]['raa'],rhod[b]['vza'], rhod[b]['sza'], ta))
                bv = ac.shared.rsr_convolute(ret, waves_mu, rsr[b]['response'], rsr[b]['wave'])
                ratm_bands[b].append(bv)

        ## interpolate rhod to taua steps in each band
        taua_all = []
        rhod_list = False
        for b in rhod_sel:
            cur_tau = np.interp(rhod_sel[b], ratm_bands[b], taua_steps)
            if False:
                if rhod_sel[b].shape != (): ## if multiple rhod are given
                    rhod_list = True
                    cur_tau[cur_tau == min(taua_steps)] = np.nan
                else: ## single value per band
                    if cur_tau == min(taua_steps): cur_tau = np.nan
            taua_all.append(cur_tau)
            
        if rhod_sel[b].shape != (): rhod_list = True
        if rhod_list:
            print('rhod list not yet implemented')
            print(ip)
            tmp = [rgi((ip, waves_mu[ti], raa, vza, sza, t)) for ti, t in enumerate(taua_all)]
            #return(taua_all, rgi, ip, waves_mu, raa, vza, sza)
            return(tmp, rhod_sel, taua_all)
        
            sub_idx = np.argsort(taua_all_sel)
            sel_idx = sub_idx[0]
            taua_sel = taua_all_sel[sel_idx]
        else:
            ## exclude bands and find lowest taua
            rhod_sub_idx = [iw for iw, b in enumerate(rhod) if 
                            (rhod[b]['tt_gas'] > rhod_tgas_cutoff) & (rhod[b]['wave'] >= rhod_min_wave)]
            taua_sub = [taua_all[iw] for iw in rhod_sub_idx]
            #sub_idx = np.argsort(taua_all)
            #sel_idx = sub_idx[0]
            sub_idx = np.argsort(taua_sub)
            sel_idx = rhod_sub_idx[sub_idx[0]]
            taua_sel = taua_all[sel_idx]

        ## in case only 1 band passes the exclusion criteria
        ## rhod_tgas_cutoff & rhod_min_wave
        if len(sub_idx) == 1:
            rhod_fit_bands = 1
            print('Setting rhod fit bands to 1')

        ## get ratm for this taua and convolute to the sensor to get rmsd fit
        #tmp = rgi((ip, waves_mu, raa, vza, sza, taua_sel))
        #ratm_sel = [ac.shared.rsr_convolute(tmp, waves_mu, rsr[b]['response'], rsr[b]['wave']) for b in rhod]
        ratm_sel = []
        for b in rhod:
            tmp  = rgi((ip, waves_mu, rhod[b]['raa'],rhod[b]['vza'], rhod[b]['sza'], taua_sel)) 
            ratm_sel.append(ac.shared.rsr_convolute(tmp, waves_mu, rsr[b]['response'], rsr[b]['wave']))
            
        if rhod_fit_bands == 2:
            ## find best fit for two bands
            if rhod_fit_selection == 'min_tau':
                ## get rmsd fit with 2 bands giving lowest (non-zero) aot
                sub_fit = [sub_idx[1]]
            elif rhod_fit_selection == 'min_drmsd':
                ## get rmsd fit with 2 best fitting bands
                sub_fit = sub_idx[1:]

            ## find rmsd bands
            sel_rmsd = 1.0
            rhod_list = [rhod_sel[b] for b in rhod]
            for si, iw in enumerate(sub_fit):
                rmsd = ac.shared.rmsd([rhod_list[sel_idx], rhod_list[iw]], 
                                      [ratm_sel[sel_idx], ratm_sel[iw]])
                if rmsd < sel_rmsd:
                    sel_idx2 = rhod_sub_idx[si]
                    sel_rmsd = rmsd*1.0
        else:
            rhod_list = [rhod_sel[b] for b in rhod]
            ratm_list = [ratm_sel[ib] for ib, b in enumerate(rhod)]
            sel_rmsd = ac.shared.rmsd(rhod_list, ratm_list)
            sel_idx2 = -1

        ## store fit
        model_results[lutid] = {'taua': taua_sel, 'rmsd':sel_rmsd, 
                                'sel_idx': sel_idx, 'sel_idx2': sel_idx2, 
                                'pressure':pressure, 'lutid':lutid if pressure is None else lutidp,
                                'rhod_fit_bands':rhod_fit_bands, 'rhod_fit_selection':rhod_fit_selection,
                                'lut_dim':lut_dim, 'rgi':rgi, 'lut_meta':lut_meta, 'ratm_sel':ratm_sel}
        print(lutid, taua_sel, sel_rmsd, sel_idx, sel_idx2)
        
    ## choose model
    sel_model = None
    ## model with best fit
    if rhod_model_selection == 'min_drmsd':
        sel_rmsd = 1.0
        for lut in model_results:
            if model_results[lut]['rmsd'] <= sel_rmsd:
                sel_rmsd = model_results[lut]['rmsd']
                sel_model = lut
    ## or model giving the lowest tau
    elif rhod_model_selection == "min_tau":
        sel_tau = 10.
        for lut in model_results:
            if model_results[lut]['taua'] <= sel_tau:
                sel_tau = model_results[lut]['taua']
                sel_model = lut

    if sel_model is None: sel_model = luts[0]

    ## get ac parameters for selected model
    result = model_results[sel_model]
    result['rhod']=rhod
    result['rhod_sel']=rhod_sel
    
    ## compute all parameters for selected model
    if False: 
        waves_mu = lutd[sel_model]['meta']['wave']
        for ip, par in enumerate(lutd[sel_model]['meta']['par']):
            tmp = result['rgi']((ip, waves_mu, raa, vza, sza, result['taua']))
            par_sel = [ac.shared.rsr_convolute(tmp, waves_mu, rsr[b]['response'], rsr[b]['wave']) for b in rsr_bands]
            result[par] = {b:par_sel[ib] for ib, b in enumerate(rsr_bands)}
        ## remove interpolator from result
        for k in ['rgi', 'lut_dim']: del result[k]
    
    return(result)
