
import os,timeit
import numpy as np
import pandas as pd
from pyqsofit.PyQSOFit import QSOFit
from astropy.io import fits
import matplotlib.pyplot as plt
import warnings
import emcee # optional, for MCMC


#preset parameters
rej_abs_line = False
rej_abs_conti = False
decompose_host = True
wave_mask = None
poly = True

warnings.filterwarnings("ignore")

QSOFit.set_mpl_style()

path_ex = os.getcwd()
saveDirectory = "results/pyQSOFit_set2_test"
pd.set_option('display.max_columns', 100)

DESI_HSC = pd.read_csv(r'data/DESI_HSC_QSO_noSDSS_z_0.8_20pc.csv')
DESI_HSC['HSC_ID'] = DESI_HSC['HSC_ID'].astype(str)
DESI_HSC['DESI_ID'] = DESI_HSC['DESI_ID'].astype(str)

#path_out = '/Users/tang/Documents/Binary_quasar/data/targets/sp_data/rescale/fitting/'
#os.chdir('/Users/tang/Documents/red_quasars/DESI/DESI_QSO') # not sure where to set this, leaving out for now

obj_list = DESI_HSC['DESI_ID'][:10]
for DESI_id in obj_list:
    obj_df = DESI_HSC[DESI_HSC['DESI_ID'] == DESI_id]
    data = pd.read_csv('data/spectra_set2/DESI_HSC_noSDSS_z_0.8_20pc/'+'%s_spec.csv'%DESI_id,index_col=0)
    
    #decide the threshold of line width to generate line results
    broad_fwhm = 1000
    
    
    lam = data['wave']       # OBS wavelength [A]
    flux = data['flux']            # OBS flux [erg/s/cm^2/A]
    err = data['flux_err']   # 1 sigma error
    z = np.float64(obj_df['Z_DESI'])
        
    #let's do decomposition for all the sources for now
    decompose_host = z<0.8
    
    # Prepare data
    q_mle = QSOFit(lam, flux, err, z, path=path_ex) #plateid=plateid, mjd=mjd, fiberid=fiberid,
    
    # Double check the installation path with the PCA / Fe template files
    print('install path:', q_mle.install_path)
    # Change it if you installed them somewhere else
    #q_mle.install_path = '...'
    
    start = timeit.default_timer()
    # Do the fitting
    try:
        q_mle.Fit(name='%s' %DESI_id,  # customize the name of given targets. Default: plate-mjd-fiber
                # prepocessing parameters
                nsmooth=1,  # do n-pixel smoothing to the raw input flux and err spectra
                and_mask=False,  # delete the and masked pixels
                or_mask=False,  # delete the or masked pixels
                reject_badpix=True,  # reject 10 most possible outliers by the test of pointDistGESD
                deredden=True,  # correct the Galactic extinction
                wave_range=None,  # trim input wavelength
                wave_mask=wave_mask,  # 2-D array, mask the given range(s)
        
                # host decomposition parameters
                decompose_host=decompose_host,  # If True, the host galaxy-QSO decomposition will be applied
                host_prior=True, # If True, the code will adopt prior-informed method to assist decomposition. Currently, only 'CZBIN1' and 'DZBIN1' model for QSO PCA are available. And the model for galaxy must be PCA too.
                host_prior_scale=0.2, # scale of prior panelty. Usually, 0.2 works fine for SDSS spectra. Adjust it smaller if you find the prior affect the fitting results too much.
        
                host_line_mask=True, # If True, the line region of galaxy will be masked when subtracted from original spectra.
                decomp_na_mask=True, # If True, the narrow line region will be masked when perform decomposition
                qso_type='CZBIN1', # PCA template name for quasar
                npca_qso=10, # numebr of quasar templates
                host_type='PCA', # template name for galaxy
                npca_gal=5, # number of galaxy templates
                
                # continuum model fit parameters
                Fe_uv_op=True,  # If True, fit continuum with UV and optical FeII template
                poly=poly,  # If True, fit continuum with the polynomial component to account for the dust reddening
                BC=False,  # If True, fit continuum with Balmer continua from 1000 to 3646A
                initial_guess=None,  # Initial parameters for continuum model, read the annotation of this function for detail
                rej_abs_conti=rej_abs_conti,  # If True, it will iterately reject 3 sigma outlier absorption pixels in the continuum
                n_pix_min_conti=100,  # Minimum number of negative pixels for host continuuum fit to be rejected.
        
                # emission line fit parameters
                linefit=True,  # If True, the emission line will be fitted
                rej_abs_line=rej_abs_line,
                # If True, it will iterately reject 3 sigma outlier absorption pixels in the emission lines
        
                # fitting method selection
                MC=False,
                # If True, do Monte Carlo resampling of the spectrum based on the input error array to produce the MC error array
                MCMC=False,
                # If True, do Markov Chain Monte Carlo sampling of the posterior probability densities to produce the error array
                nsamp=200,
                # The number of trials of the MC process (if MC=True) or number samples to run MCMC chain (if MCMC=True)
        
                # advanced fitting parameters
                param_file_name='scripts/qsopar.fits',  # Name of the qso fitting parameter FITS file.
                nburn=20,  # The number of burn-in samples to run MCMC chain
                nthin=10,  # To set the MCMC chain returns every n samples
                epsilon_jitter=0.,
                # Initial jitter for every initial guass to avoid local minimum. (Under test, not recommanded to change)
        
                # customize the results
                # save_result=True,  # If True, all the fitting results will be saved to a fits file
                # save_fits_name=obj.file_name+'_'+suffix,  # The output name of the result fits
                save_fits_path=saveDirectory,  # The output path of the result fits
                plot_fig=True,  # If True, the fitting results will be plotted
                save_fig=True,  # If True, the figure will be saved
                plot_corner=True,  # Whether or not to plot the corner plot results if MCMC=True
        
                # debugging mode
                verbose=True,  # turn on (True) or off (False) debugging output
                
                # sublevel parameters for figure plot and emcee
                kwargs_plot={
                    'save_fig_path': saveDirectory,  # The output path of the figure
                    'broad_fwhm'   : broad_fwhm , # km/s, lower limit that code decide if a line component belongs to broad component
                    #'plot_br_prop' : True # Plot broad line properties
                },
                kwargs_conti_emcee={},
                kwargs_line_emcee={})
        qso = q_mle.qso
        host = q_mle.host
        wavelengths = q_mle.wave
        
        np.savetxt(saveDirectory+f'/{DESI_id}_data.csv',np.transpose([wavelengths,qso,host]),delimiter=",")
    except:
        continue
    end = timeit.default_timer()
                    
    print(f'Fitting finished in {np.round(end - start, 1)}s')
    df_q_mle = pd.DataFrame({'idx':q_mle.line_result_name,'value':q_mle.line_result})
    df_q_mle.to_csv(saveDirectory+'/%s_line.csv' %DESI_id)

    df_cont_q_mle = pd.DataFrame({'idx':q_mle.conti_result_name,'value':q_mle.conti_result})
    df_cont_q_mle.to_csv(saveDirectory+'/%s_cont.csv' %DESI_id,index='idx')
    
    