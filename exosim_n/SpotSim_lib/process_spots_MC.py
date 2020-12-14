# -*- coding: utf-8 -*-
"""
Created on Thu May 24 06:21:29 2018

@author: user1

Process monte carlo results from spots

"""

import numpy as np
import matplotlib.pyplot as plt

#plt.clf()
res_file = '/Users/user1/Desktop/spot_results'
res_file = '/Users/user1/Desktop/spot_results_HD 209458 b'

spatDist = 'Uniform'; col='b'; shape = 'o'
#spatDist = 'Longitudinal'; col='r'; shape = 's'
#spatDist = 'Equatorial'; col = 'g'; shape = 'D'
##




b = 0.37
b = 0.53

#b = 0.0

filling = 1.0
ch_list= ['AIRS CH1'] 
#ch_list = ['NIR Spec']
#ch_list= ['FGS Red'] 
ch_list= ['AIRS CH0', 'AIRS CH1', 'NIR Spec', 'FGS Red', 'FGS Prime', 'NIR Phot'] 
#ch_list= ['AIRS CH1','NIR Spec', 'FGS Red', 'FGS Prime', 'NIR Phot'] 


spatDist_list = ['Uniform', 'Equatorial', 'Longitudinal']


filling_list = [0.1, 1.0, 5.0, 10.0, 20.0, 30.0]
b_list = [0.37,0.0]
b_list = [0.53, 0.0]

input_cr = 0.0135
input_cr = 0.01464


col_list = ['b', 'r', 'g', 'm', 'y', 'c']
shape_list = ['o', 's', 'D', 'v', '^', 'p']

zip_list = zip(filling_list, col_list, shape_list)

 
n_type = 'SN'
variation =0
tag = 'spot_facs'

#n_type = 'PN'
#filling_list = [0.1]
#b_list = [0.0]
#spatDist_list = ['Uniform']
#col_list = ['b']
#shape_list = ['o']
#zip_list = zip(filling_list, col_list, shape_list)

f_p_mean_list = []
f_p_err_list = []
wl_list=[]

fig = plt.figure()
xxx= 0
for b in b_list:
    for spatDist in spatDist_list:
    
    
#        plt.figure('%s %s %s variation: %s, b =%s '%(spatDist, tag, n_type, variation, b))  
        xxx = xxx+1
        plot = 320 + xxx
        ax = fig.add_subplot(plot)
        for filling, col, shape in zip_list:
            for ch in ch_list:
                f = np.load('%s/%s/b=%s/%s_%s_fill_%s_var_%s_%s.npy'%(res_file, spatDist, b, ch, n_type, filling,  variation , tag ) )    
                f_p = np.load('%s/Uniform/b=0.0/%s_PN_fill_0.1_var_0_spot_facs.npy'%(res_file, ch) )    
#                f_p=f
                wl = np.load('%s/%s/b=%s/%s_wavelength.npy'%(res_file, spatDist, b, ch))
                if ch == 'AIRS CH0':
                    wl = wl[8:-1]
                    f = f[:,8:-1]
                    f_p = f_p[:,8:-1]
                elif ch == 'AIRS CH1':
                    wl = wl[1:]
                    f = f[:,1:] 
                    f_p=f_p[:,1:] 
                elif ch == 'NIR Spec':
                    wl = wl[2:-4]                    
                    f = f[:,2:-4]  
                    f_p=f_p[:,2:-4]  
                
                if ch == ch_list[0]:
#                    if filling == zip_list[0][0]:
#                        plt.plot(wl, f_p.mean(axis=0), 'kd', label = 'Photon noise')                    
                    plt.plot(wl, f.mean(axis=0), '%s%s'%(col,shape), label = filling)

            
                else:
#                    if filling == zip_list[0][0]:                    
#                        plt.plot(wl, f_p.mean(axis=0), 'kd')                       
                    plt.plot(wl, f.mean(axis=0), '%s%s'%(col,shape))
    
#                if filling == zip_list[0][0]:                    
#                   plt.errorbar(wl, f_p.mean(axis=0),  f_p.std(axis=0), ecolor = 'k', fmt=None)
#                   plt.fill_between(wl, f_p.mean(axis=0)-f_p.std(axis=0), f_p.mean(axis=0)+ 0.1)
            
                plt.errorbar(wl, f.mean(axis=0),  f.std(axis=0), ecolor = col, fmt=None)
                
                if filling == zip_list[0][0] and b == b_list[0] and spatDist == spatDist_list[0]:  
                    f_p_mean_list= np.hstack((f_p_mean_list, f_p.mean(axis=0)))
                    f_p_err_list = np.hstack((f_p_err_list, f_p.std(axis=0)))
                    wl_list = np.hstack((wl_list, wl))
 
         
 
#        plt.fill_between(wl_list, f_p_mean_list-f_p_err_list, f_p_mean_list+f_p_err_list,zorder=0)   
#        plt.fill_between(wl_list, f_p_mean_list-0.001, f_p_mean_list+0.001,zorder=0)      
#        plt.errorbar(wl_list, f_p_mean_list, f_p_err_list, ecolor = 'k', fmt = None,zorder=0)   
 
        idx = wl_list.argsort()
        wl_list = wl_list[idx[::-1]]
        f_p_mean_list = f_p_mean_list[idx[::-1]]
        f_p_err_list = f_p_err_list[idx[::-1]]
        
        plt.fill_between(wl_list, f_p_mean_list+f_p_err_list, f_p_mean_list-f_p_err_list, facecolor='black', alpha=0.5, zorder=0)   

        
        plt.title('%s %s %s variation: %s, b =%s '%(spatDist, tag, n_type, variation, b))
        plt.plot([0.5,8.0], [input_cr,input_cr], 'k--', linewidth=2, label = 'input contrast ratio')
        
        
        plt.grid(True)
        #plt.ylim(0.0133,0.0137)
        plt.xlim(0.5, 8.0)
        ax = plt.gca()
        legend = ax.legend(loc='lower right', shadow=True, labelspacing =0.2)
        frame = legend.get_frame()
        frame.set_facecolor('0.90')
        for label in legend.get_texts():
        #    label.set_fontsize('medium')
            label.set_fontsize(12)
        for label in legend.get_lines():
            label.set_linewidth(1.5)  # the legend line width
        for item in ([ax.title, ax.xaxis.label, ax.yaxis.label] +
                     ax.get_xticklabels() + ax.get_yticklabels()):
            item.set_fontsize(16)