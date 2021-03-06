"""
exosim_n
Read and display results
"""

import pickle
import numpy as np
import matplotlib.pyplot as plt
plt.style.use('classic')
import exosim_n
from exosim_n.classes.params import Params
from exosim_n.classes.options import Options
from scipy import interpolate
import os

def run(results_file):
    
    
    DEBUG= Options.DEBUG
    results_file='FIXED_NAME.pickle'
    if DEBUG:
        if DEBUG:plt.show() 
        showSignal=Options.showSignal
        showExposure=Options.showExposure

        showBadPixels=Options.showBadPixels
        results_file='FIXED_NAME.pickle'
    # exosim_n_path =  os.path.dirname((os.path.dirname(exosim_n.__file__)))
    exosim_n_path = os.path.abspath(os.path.join(os.path.dirname( __file__ ), '../..'))
    
    paths_file ='%s/exosim_n/input_files/exosim_n_paths.txt'%(exosim_n_path)
        
    params_to_opt = Params(0, paths_file, 3)  
    paths = params_to_opt.params
    output_directory = paths['output_directory']
    if output_directory == '':  
        output_directory ='%s/output'%(exosim_n_path)
      
    results_file =  '%s/%s'%(output_directory, results_file)
 
    with open(results_file, 'rb') as handle:
       res_dict = pickle.load(handle)
        
    no_list = np.array([ 'All noise','All photon noise','Source photon noise','Dark current noise',
                        'Zodi noise','Emission noise','Read noise','Spatial jitter noise',
                        'Spectral jitter noise','Combined jitter noise','No noise - no background','No noise - all background'])  
    color = ['0.5','r', 'b','k','orange','pink', 'y','g','purple','r', '0.8','c']
              
    ch =res_dict['ch']
    
    sim_text = '%s.txt'%(results_file)
    
    with open(sim_text) as f:
        content = f.readlines()
        content = [x.strip() for x in content]        
        for i in range(len(content)):       
            if content[i] != '' and content[i][0] != '#':
                aa = content[i].split()
                if aa[0] == 'Wavelength:':
                    wavlim=[np.float(aa[1]), np.float(aa[2])]

        
    if res_dict['simulation_mode'] == 2:

            wl = res_dict['wl']     
            idx = np.argwhere ((res_dict['wl']>=wavlim[0])&(res_dict['wl']<=wavlim[1])).T[0]
            wl =res_dict['wl'][idx]
            if res_dict['simulation_realisations'] > 1:   
                p_stack = res_dict['p_stack'][:,idx]
            else:
                p_stack = res_dict['p_stack'][idx]
            p_mean = res_dict['p_mean'][idx]
            p_std = res_dict['p_std'][idx]
            
            
            cr = res_dict['input_spec']
            cr_wl = res_dict['input_spec_wl']
         
            idx0 = np.argwhere ((np.array(cr_wl)>=wavlim[0])&(np.array(cr_wl)<=wavlim[1])).T[0]   
            cr = cr[idx0]
            cr_wl = cr_wl[idx0]      
            
       
            plt.figure('spectrum %s'%(res_dict['time_tag']))
            if res_dict['simulation_realisations'] > 1:    
                for i in range(p_stack.shape[0]):
                    plt.plot(wl,p_stack[i], '.', color='0.5', alpha =0.2)         
        
            plt.plot(wl,p_mean, 'o-', color='0.5', label = 'mean recovered spectrum')
            plt.errorbar(wl,p_mean,p_std, ecolor='0.5')
            f = interpolate.interp1d(cr_wl,cr, bounds_error=False)
            in_spec = np.array(f(wl))
            plt.plot(wl,in_spec, 'o-', color='r', linewidth=2, label='input spectrum')
            plt.legend(loc='best')
            plt.ylabel('Contrast ratio')
            plt.xlabel('Wavelength ($\mu m$)')
            plt.grid(True)
            
            plt.figure('Residual bias between mean recover spectrum and input spectrum %s'%(res_dict['time_tag']))
            f = interpolate.interp1d(cr_wl,cr, bounds_error=False)
            in_spec = np.array(f(wl))
            plt.plot(wl,1e6*(p_mean-in_spec), 'o-', color='b', label = 'residual bias')
            plt.legend(loc='best')
            plt.ylabel('Residual bias (ppm)')
            plt.xlabel('Wavelength ($\mu m$)')
            plt.grid(True) 
           
            
            plt.figure('percent residual bias between mean recover spectrum and input spectrum %s'%(res_dict['time_tag']))
            f = interpolate.interp1d(cr_wl,cr, bounds_error=False)
            in_spec = np.array(f(wl))
            plt.plot(wl,100*(p_mean-in_spec)/in_spec, 'o-', color='b', label = 'residual bias')
            plt.legend(loc='best')
            plt.ylabel('Percent difference from input value')
            plt.xlabel('Wavelength ($\mu m$)')
            plt.grid(True)
                     
            p_std[np.isnan(p_std)] = 0          
            plt.figure('precision %s'%(res_dict['time_tag']))
            plt.ylabel('1 sigma error on transit depth (ppm)')
            plt.xlabel('Wavelength ($\mu m$)')
            plt.plot(wl, p_std*1e6, 'bo', alpha=0.5)       
                     
             
            if len(wl)>1:
                r = 4 
                z = np.polyfit(wl, p_std*1e6, r)
                p= np.poly1d(z)
                # yhat = p(wav)
                # ybar = sum(p_std)/len(p_std)
                # SST = sum((p_std - ybar)**2)
                # SSreg = sum((yhat - ybar)**2)
                # R2 = SSreg/SST  
                y =0
                for i in range (0,r+1):
                    y = y + z[i]*wl**(r-i) 
            else: 
                y = p_std*1e6
                
            plt.figure('precision %s'%(res_dict['time_tag']))
            plt.plot(wl, y, '-', color='r', linewidth=2) 
            plt.grid(True)
            
                
            for ntransits in [1,10,100]:
                f = interpolate.interp1d(cr_wl,cr, bounds_error=False)
                rand_spec = np.array(f(wl))
                if ntransits == 1:
                    plt.figure('sample spectrum for %s transit %s'%(ntransits, res_dict['time_tag']))  
                else:
                    plt.figure('sample spectrum for %s transits %s'%(ntransits, res_dict['time_tag']))  
                plt.plot(cr_wl,cr, '-', color='r', linewidth=2, label='input spectrum')
                for i in range(len(wl)):
                    rand_spec[i] = np.random.normal(rand_spec[i], y[i]/1e6/np.sqrt(ntransits))
                plt.plot(wl, rand_spec, 'o-', color='b', label = 'randomized spectrum')
                plt.errorbar(wl, rand_spec, y/1e6/np.sqrt(ntransits), ecolor='b')
                plt.legend(loc='best')
                plt.ylabel('Contrast ratio')
                plt.xlabel('Wavelength ($\mu m$)')
                plt.grid(True)



            if not DEBUG or (DEBUG and showBadPixels):
                plt.figure('bad pixels %s'%(res_dict['time_tag']))          
                plt.imshow(res_dict['bad_map'], interpolation='none', aspect='auto')
                ticks = np.arange(res_dict['bad_map'].shape[1])[0::int(res_dict['bad_map'].shape[1]/10)]  
                ticklabels =  np.round(res_dict['pixel wavelengths'], 2)[0::int(res_dict['bad_map'].shape[1]/10)]  
                plt.xticks(ticks=ticks, labels = ticklabels)
                plt.ylabel('Spatial pixel')
                plt.xlabel('Wavelength ($\mu m$)')
                
            
            plt.figure('example integration image %s'%(res_dict['time_tag']))          
            plt.imshow(res_dict['example_exposure_image'], interpolation='none', aspect='auto')
            ticks = np.arange(res_dict['example_exposure_image'].shape[1])[0::int(res_dict['example_exposure_image'].shape[1]/10)]  
            ticklabels =  np.round(res_dict['pixel wavelengths'], 2)[0::int(res_dict['example_exposure_image'].shape[1]/10)]  
            plt.xticks(ticks=ticks, labels = ticklabels)
            plt.ylabel('Spatial pixel')
            plt.xlabel('Wavelength ($\mu m$)')           
            cbar = plt.colorbar() 
            cbar.set_label('Count (e$^-$)',size=12)
            plt.xlabel('Wavelength ($\mu m$)')
            
            if 'realization_0_binned_lc' in res_dict.keys():
                plt.figure('example light curve  from first realisation')          
                lc = res_dict['realization_0_binned_lc'][:, int(res_dict['realization_0_binned_lc'].shape[1]/2)]
                wav = res_dict['wl'][int(res_dict['realization_0_binned_lc'].shape[1]/2)]
                time = res_dict['exp_end_time'] 
                plt.plot(time, lc, label = f'{np.round(wav,2)} microns')
                plt.legend(loc='best')
                plt.grid()
                plt.xlabel('Time (sec)')
                plt.ylabel('Signal (e$^-$)')   
                            
     
                
    elif res_dict['simulation_mode'] == 1:
       
        no_dict =  res_dict['noise_dic']  
  
        for key in no_dict.keys():
   
            idx = np.argwhere(no_list==key)[0].item()
            col = color[idx]
            
            idx = np.argwhere ((no_dict[key]['wl']>=wavlim[0])&(no_dict[key]['wl']<=wavlim[1])).T[0]          

            noise_type = key
            wl = no_dict[key]['wl'][idx]
            
            print (idx, wl)
          
            if res_dict['simulation_realisations'] == 1:          
                sig_stack = no_dict[key]['signal_mean_stack'][idx]
                no_stack = no_dict[key]['signal_std_stack'][idx]
                
              
                if 'fracNoT14_mean' in no_dict[key].keys():
                    fracNoT14_stack = no_dict[key]['fracNoT14_stack'][idx]
            else:
                sig_stack = no_dict[key]['signal_mean_stack'][:,idx]
                no_stack = no_dict[key]['signal_std_stack'][:,idx]
                if 'fracNoT14_mean' in no_dict[key].keys():
                    fracNoT14_stack = no_dict[key]['fracNoT14_stack'][:,idx]
            
            sig_mean = no_dict[key]['signal_mean_mean'][idx]
            no_mean = no_dict[key]['signal_std_mean'][idx]
            
            
            
            
            if 'fracNoT14_mean' in no_dict[key].keys():
                fracNoT14_mean = no_dict[key]['fracNoT14_mean'][idx]
            if  not DEBUG or (DEBUG and showSignal):
                plt.figure('signal %s'%(res_dict['time_tag']))
                plt.plot(wl,sig_mean, 'o-', color = col, label = noise_type)
                if res_dict['simulation_realisations'] > 1:    
                    for i in range(sig_stack.shape[0]):
                        plt.plot(wl,sig_stack[i], ':', color = col, alpha=0.5)           
                plt.legend(loc='best', ncol = 3, borderpad =0.3, fontsize=10)
                plt.ylabel('Signal (e$^-$)')
                plt.xlabel('Wavelength ($\mu m$)')
                plt.grid(True)

                
                plt.figure('noise %s'%(res_dict['time_tag']))

                plt.plot(wl,no_mean, 'o-', color = col, label = noise_type)
                if res_dict['simulation_realisations'] > 1:
                    for i in range(no_stack.shape[0]):
                        plt.plot(wl,no_stack[i], '.', color = col, alpha=0.5)           
                plt.legend(loc='best', ncol = 3, borderpad =0.3, fontsize=10)
                plt.ylabel('Noise (e$^-$)')
                plt.xlabel('Wavelength ($\mu m$)')
                plt.grid(True)
                #if DEBUG:plt.show()      
                if res_dict['simulation_realisations'] > 1:
                    for i in range(no_stack.shape[0]):
                        plt.plot(wl,no_stack[i], '.', color = col, alpha=0.5)                 

                

                #if DEBUG:plt.show()   
                if 'fracNoT14_mean' in no_dict[key].keys():
                    plt.figure('fractional noise %s'%(res_dict['time_tag']))
                    plt.plot(wl,fracNoT14_mean, 'o-', color = col, label = noise_type)
                    if res_dict['simulation_realisations'] > 1:
                        for i in range(fracNoT14_stack.shape[0]):
                            plt.plot(wl,fracNoT14_stack[i], '.', color = col, alpha=0.5)   
                                
                    plt.legend(loc='best', ncol = 3, borderpad =0.3, fontsize=10)
                    plt.ylabel('Fractional noise at T14 (ppm)')
                    plt.xlabel('Wavelength ($\mu m$)')
                    plt.grid(True)
                    if DEBUG:plt.show()      #SHOW 1: Signal + Noise
                    plt.ylim(fracNoT14_mean.min() - fracNoT14_mean.min()*0.2
                            , fracNoT14_mean.max() + fracNoT14_mean.max()*0.2)
                
                                
                    plt.figure('precision %s'%(res_dict['time_tag']))
                    plt.plot(wl,fracNoT14_mean*np.sqrt(2), 'o', color = col, label = noise_type, alpha=0.5)
                    plt.legend(loc='best', ncol = 3, borderpad =0.3, fontsize=10)
                    plt.ylabel('1$\sigma$ error on transit depth (ppm)')
                    plt.xlabel('Wavelength ($\mu m$)')
                    plt.ylim(fracNoT14_mean.min()*np.sqrt(2) - fracNoT14_mean.min()*np.sqrt(2)*0.2
                            , fracNoT14_mean.max()*np.sqrt(2) + fracNoT14_mean.max()*np.sqrt(2)*0.2)
        
                    if len(wl)>1:
                        r = 4 
                        
                        z = np.polyfit(wl, fracNoT14_mean*np.sqrt(2), r)
                        p= np.poly1d(z)
                        # yhat = p(wav)
                        # ybar = sum(p_std)/len(p_std)
                        # SST = sum((p_std - ybar)**2)
                        # SSreg = sum((yhat - ybar)**2)
                        # R2 = SSreg/SST  
                        y =0
                        for i in range (0,r+1):
                            y = y + z[i]*wl**(r-i) 
                            
                    else:
                        y = fracNoT14_mean*np.sqrt(2)
                    
                        
                    plt.plot(wl, y, '-', color='r', linewidth=2) 
                    plt.grid(True)
                    
                    cr = res_dict['input_spec']
                    cr_wl = res_dict['input_spec_wl']
                    
                    print (cr)
    
                    idx0 = np.argwhere ((np.array(cr_wl)>=wavlim[0])&(np.array(cr_wl)<=wavlim[1])).T[0]   
                    cr = cr[idx0]
                    cr_wl = cr_wl[idx0]
    
                    f = interpolate.interp1d(cr_wl,cr, bounds_error=False)
                    rand_spec = np.array(f(wl))
                
                
                    for ntransits in [1,10,100]:
                        plt.figure('sample spectrum for %s transit %s'%(ntransits, res_dict['time_tag']))  
                        plt.plot(cr_wl,cr, '-', color='r', linewidth=2, label='input spectrum')
                        for i in range(len(wl)):
                            rand_spec[i] = np.random.normal(rand_spec[i], y[i]/1e6/np.sqrt(ntransits))
                        plt.plot(wl, rand_spec, 'o-', color='b', label = 'randomized spectrum')
                        plt.errorbar(wl, rand_spec, y/1e6/np.sqrt(ntransits), ecolor='b')
                        plt.legend(loc='best')
                        plt.ylabel('Contrast ratio')
                        plt.xlabel('Wavelength ($\mu m$)')
                        plt.grid(True)
                    if DEBUG:plt.show()    #Sample spectrum with different transits
            if  not DEBUG or  (DEBUG and showBadPixels):
                plt.figure('bad pixels %s'%(res_dict['time_tag']))          
                plt.imshow(no_dict[key]['bad_map'], interpolation='none', aspect='auto')
                ticks = np.arange(no_dict[key]['bad_map'].shape[1])[0::int(no_dict[key]['bad_map'].shape[1]/10)]  
                ticklabels =  np.round(no_dict[key]['pixel wavelengths'], 2)[0::int(no_dict[key]['bad_map'].shape[1]/10)]  
                plt.xticks(ticks=ticks, labels = ticklabels)
                plt.ylabel('Spatial pixel')
                plt.xlabel('Wavelength ($\mu m$)')
            if  not DEBUG or  (DEBUG and showExposure):
                plt.figure('example integration image %s'%(res_dict['time_tag']))          
                plt.imshow(no_dict[key]['example_exposure_image'], interpolation='none', aspect='auto', vmin=0, vmax=no_dict[key]['example_exposure_image'].max(), cmap='jet')
                ticks = np.arange(no_dict[key]['example_exposure_image'].shape[1])[0::int(no_dict[key]['example_exposure_image'].shape[1]/10)]  
                ticklabels =  np.round(no_dict[key]['pixel wavelengths'], 2)[0::int(no_dict[key]['example_exposure_image'].shape[1]/10)]  
                plt.xticks(ticks=ticks, labels = ticklabels)
                plt.ylabel('Spatial pixel')
                plt.xlabel('Wavelength ($\mu m$)')           
                cbar = plt.colorbar() 
                cbar.set_label('Count (e$^-$)',size=12)
                if DEBUG:plt.show()    #SHOW : integration image and bad pixel

    elif res_dict['simulation_mode'] == 3 or res_dict['simulation_mode'] == 4:
        no_dict =  res_dict['noise_dic']  
  
        for key in no_dict.keys():
   
            idx = np.argwhere(no_list==key)[0].item()
            col = color[idx]
            
            idx = np.argwhere ((no_dict[key]['wl']>=wavlim[0])&(no_dict[key]['wl']<=wavlim[1])).T[0]          

            noise_type = key
            wl = no_dict[key]['wl'][idx]
            
            if res_dict['simulation_realisations'] == 1:          
                sig_stack = no_dict[key]['signal_mean_stack'][idx]
                no_stack = no_dict[key]['signal_std_stack'][idx]
                if 'fracNoT14_mean' in no_dict[key].keys():
                    fracNoT14_stack = no_dict[key]['fracNoT14_stack'][idx]
            else:
                sig_stack = no_dict[key]['signal_mean_stack'][:,idx]
                no_stack = no_dict[key]['signal_std_stack'][:,idx]
                if 'fracNoT14_mean' in no_dict[key].keys():
                    fracNoT14_stack = no_dict[key]['fracNoT14_stack'][:,idx]
            
            sig_mean = no_dict[key]['signal_mean_mean'][idx]
            no_mean = no_dict[key]['signal_std_mean'][idx]
            if 'fracNoT14_mean' in no_dict[key].keys():
                fracNoT14_mean = no_dict[key]['fracNoT14_mean'][idx]
    
            plt.figure('signal %s'%(res_dict['time_tag']))
            plt.plot(wl,sig_mean, 'o-', color = col, label = noise_type)
            if res_dict['simulation_realisations'] > 1:    
                for i in range(sig_stack.shape[0]):
                    plt.plot(wl,sig_stack[i], ':', color = col, alpha=0.5)           
            plt.legend(loc='best', ncol = 2, borderpad =0.3, fontsize=10)
            plt.ylabel('Signal (e$^-$)')
            plt.xlabel('Wavelength ($\mu m$)')
            plt.grid(True)
     
            plt.figure('noise %s'%(res_dict['time_tag']))
            plt.plot(wl,no_mean, 'o-', color = col, label = noise_type)
            if res_dict['simulation_realisations'] > 1:
                for i in range(no_stack.shape[0]):
                    plt.semilogy(wl,no_stack[i], '.', color = col, alpha=0.5)           
            plt.legend(loc='best', ncol = 3, borderpad =0.3, fontsize=10)
            plt.ylabel('Noise (e$^-$)')
            plt.xlabel('Wavelength ($\mu m$)')
            plt.grid(True)        
            
            # plt.figure('sigp')
            # N = 2367
            # sigp = np.sqrt(2)*(no_mean/sig_mean)/np.sqrt(N)
            # plt.plot(wl, sigp*1e6)
            
  
            
            if res_dict['simulation_realisations'] > 1:
                for i in range(no_stack.shape[0]):
                    plt.plot(wl,no_stack[i], '.', color = col, alpha=0.5)                 
            
            if 'fracNoT14_mean' in no_dict[key].keys():
                plt.figure('fractional noise %s'%(res_dict['time_tag']))
                plt.semilogy(wl,fracNoT14_mean, 'o-', color = col, label = noise_type)
                if res_dict['simulation_realisations'] > 1:
                    for i in range(fracNoT14_stack.shape[0]):
                        plt.plot(wl,fracNoT14_stack[i], '.', color = col, alpha=0.5)           
                plt.legend(loc='best', ncol = 3, borderpad =0.3, fontsize=10)
                plt.ylabel('Fractional noise at T14 (ppm)')
                plt.xlabel('Wavelength ($\mu m$)')
                plt.grid(True)
        

            if key == 'All noise':
              if not DEBUG or (DEBUG and showBadPixels): 
                plt.figure('bad pixels %s'%(res_dict['time_tag']))          
                plt.imshow(no_dict[key]['bad_map'], interpolation='none', aspect='auto')
                ticks = np.arange(no_dict[key]['bad_map'].shape[1])[0::int(no_dict[key]['bad_map'].shape[1]/10)]  
                ticklabels =  np.round(no_dict[key]['pixel wavelengths'], 2)[0::int(no_dict[key]['bad_map'].shape[1]/10)]  
                plt.xticks(ticks=ticks, labels = ticklabels)
                plt.ylabel('Spatial pixel')
                plt.xlabel('Wavelength ($\mu m$)')
                
                plt.figure('example integration image %s'%(res_dict['time_tag']))          
                plt.imshow(no_dict[key]['example_exposure_image'], interpolation='none', aspect='auto', vmin=0, vmax=no_dict[key]['example_exposure_image'].max(), cmap='jet')
                ticks = np.arange(no_dict[key]['example_exposure_image'].shape[1])[0::int(no_dict[key]['example_exposure_image'].shape[1]/10)]  
                ticklabels =  np.round(no_dict[key]['pixel wavelengths'], 2)[0::int(no_dict[key]['example_exposure_image'].shape[1]/10)]  
                plt.xticks(ticks=ticks, labels = ticklabels)
                plt.ylabel('Spatial pixel')
                plt.xlabel('Wavelength ($\mu m$)')           
                cbar = plt.colorbar() 
                cbar.set_label('Count (e$^-$)',size=12)
            
    plt.show()

if __name__ == "__main__":     

    run('OOT_SNR_AIRS CH1_GJ 1214 b_2021_03_05_1332_31.pickle')
