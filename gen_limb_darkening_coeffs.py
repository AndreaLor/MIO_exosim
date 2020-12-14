"""

"""
import numpy as np
import matplotlib.pyplot as plt
import os
import pickle
from exotethys import sail


class Anc():
    def __init___(self):
        pass
    
 
 
default_path =  '/Users/user1/Desktop/ExoSim-A/exosim_a'

exotethys_path = '/Users/user1/Downloads/ExoTETHyS-new'

dfp = '%s/data/LDC_3'%(default_path)

 
def run():

    T_s_list = np.arange(3000,7000,100)
#    T_s_list = np.arange(5000,7000,100)
    
#    above 4900 missing log 6.0

        
#    wl = np.arange(0.375, 10.0, 0.05)*1e4
    wl = np.arange(1.91, 4.0, 0.02)*1e4  # for ch0



    for i in range(len(T_s_list)):
        T_s = T_s_list[i]
        if T_s >= 4900:
            log_g_list = np.arange(3.0,6.0,0.5)
        else:
            log_g_list = np.arange(3.0,6.5,0.5)
            
            
        for j in range(len(log_g_list)):
            
            log_g = log_g_list[j]
            
            print T_s, log_g
            b = '%s/Wavelength_bins_files/Wavelength_bins.txt'%(exotethys_path)  
            open(b, 'w').close()      
            file = open(b,'w') 
            for q in range (len(wl)-1):
                file.write(str(wl[q])+'     '+str(wl[q+1])+'\n')
            file.close()  
                                  
            b = '%s/examples/input_ldc.txt'%(exotethys_path) 
            print b
            open(b, 'w').close() 
            file = open(b,'w')         
            file.write("##GENERIC \n")
            file.write("calculation_type			!grid		individual \n")
            file.write("stellar_models_grid			!Phoenix_2018	Phoenix_2012_13		!Atlas \n")
            file.write("limb_darkening_laws			claret4		!power2		!square_root	    quadratic	!linear		!gen_claret	!gen_poly \n")
            file.write("passbands				uniform_phoenix_2012_13.pass \n")
            file.write("wavelength_bins_files			Wavelength_bins.txt \n")
            file.write("user_output				basic		!complete \n")
            file.write("output_path                  . \n")
            file.write("target_names				Star \n")	 	
            file.write("star_effective_temperature		%s \n"%(T_s))	 
            file.write("star_log_gravity			%s \n"%(log_g))	 
            file.write("star_metallicity			0.0 \n")
            file.close()
                          
            sail.ldc_calculate(b)              
            a = '%s/Star_ldc.pickle'%(exotethys_path)                              
            a = pickle.load(open(a, 'rb'))
#            print a.keys()
#            print a['star_params']
#            print a['passbands'].keys()
#                                
#            for i in a['passbands'].keys():
##                                    print i, a['passbands'][i]['quadratic']['coefficients']
#                                    print i, a['passbands'][i]['claret4']['coefficients']
                    
                    #xxxx  
                    
            wl_list1 = []
            wl_list2 = []
            u1_list = []
            u2_list = []
            u3_list = []
            u4_list = []            
            
            for k in a['passbands'].keys():
                    
                #    print i, a['passbands'][i]['quadratic']['coefficients']
                    u1 =  a['passbands'][k]['quadratic']['coefficients'][0]
                    u2 =  a['passbands'][k]['quadratic']['coefficients'][1]
                    
                    u1 =  a['passbands'][k]['claret4']['coefficients'][0]
                    u2 =  a['passbands'][k]['claret4']['coefficients'][1]                    
                    u3 =  a['passbands'][k]['claret4']['coefficients'][2]  
                    u4 =  a['passbands'][k]['claret4']['coefficients'][3]  
                   
                    r = k[24:]
                    
                    idx = r.find('_')
                    r1 = np.float(r[:idx])
                    r2 = np.float(r[idx+1:])
                    
                    wl_list1.append(r1)
                    wl_list2.append(r2) 
                  
                    u1_list.append(u1)
                    u2_list.append(u2)
                    u3_list.append(u3)
                    u4_list.append(u4)                    
                    
                    
                #                print r1, r2, u1
            ll = zip(wl_list1,wl_list2, u1_list,u2_list)
            ll = zip(wl_list1,wl_list2, u1_list,u2_list, u3_list,u4_list)
            
            
            ll.sort() 
            wl_list1 = [x[0] for x in ll]
            wl_list2 = [x[1] for x in ll]   
            u1_list = [x[2] for x in ll]   
            u2_list = [x[3] for x in ll]  
            u3_list = [x[4] for x in ll]   
            u4_list = [x[5] for x in ll]              
            
            plt.figure('ldcs')         
            plt.plot(wl_list1, u1_list, 'ro')
            plt.plot(wl_list1, u2_list, 'bo') 
            plt.plot(wl_list1, u3_list, 'go')
            plt.plot(wl_list1, u4_list, 'mo')         
#            print wl_list1
#            print wl_list2
            
            wl_middle = (np.array(wl_list1) + np.array(wl_list2))/2.0
                      
            b = '%s/T=%s_logg=%s_z=0.txt'%(dfp,T_s,log_g)  
            aa = np.vstack((np.array(wl_middle)/1e4, u1_list,u2_list, u3_list, u4_list)).T
            np.savetxt(b,aa)

    b = '%s/LDC_wavelength.txt'%(dfp)
    np.savetxt(b,wl/1e4)
     
            

run() 
