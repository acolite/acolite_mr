## def plot_dark_spectrum
## plots ACOLITE retrieved dark spectrum - to be improved
## written by Quinten Vanhellemont, RBINS for the PONDER project
## 2018-03-05 (moved from acolite_ac)
## modifications: QV 2018-03-13 simplified

def plot_dark_spectrum(metadata, ds_plot, waves, ratm_s, rorayl_s, rdark, dark_idx, tau550,sel_model_lut_meta,
                      xlim = (400,2300), ylim = (0,0.14)):
                import matplotlib.pyplot as plt
                ratm = []
                rorayl = []
                dark = []
                ds_waves = []
            
                sensor = metadata['SATELLITE_SENSOR']
                if '_' in sensor:
                    sat, sen = sensor.split('_')
                else:
                    sat = sensor
                    sen = ''
                    
                for b, band in enumerate(rdark): 
                    if band not in ratm_s: continue
                    ratm.append(ratm_s[band])
                    rorayl.append(rorayl_s[band])                        
                    ds_waves.append(float(waves[b]))
                    dark.append(rdark[band])
                
                band_title = dark_idx
                                          
                plt.plot(ds_waves, dark, 'o--', color='black', label=r'$\rho_{dark}$')
                plt.plot(ds_waves, rorayl, 'o--', color='blue', label=r'$\rho_{Rayleigh}$')
                plt.plot(ds_waves, ratm, 'o--', color='red', label=r'$\rho_{path}$')
                plt.xlim(xlim)
                plt.ylim(ylim)
                plt.xlabel('Wavelength (nm)')
                plt.ylabel(r'$\rho$')

                plt.title('{}/{} {}\n{}'.format(sat, sen, metadata['TIME'].strftime('%Y-%m-%d (%H:%M UTC)'),
                                               r'$\theta_s$='+ '{:.1f} '.format(metadata['THS'])+ r'$\tau_{a}550$'+'={:.3f} (mod{}, {})'.format(tau550,sel_model_lut_meta['aermod'][0], band_title)))
                plt.legend(loc='upper right')
                plt.savefig(ds_plot, dpi=150)
                plt.close()

