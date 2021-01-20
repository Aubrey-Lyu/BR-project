#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 17 21:09:54 2020

@author: dian
"""
import matplotlib.pyplot as plt
import rpy2.robjects as ro
from rpy2.robjects.conversion import localconverter
from rpy2.robjects import pandas2ri
import numpy as np
import os
from mne.viz import circular_layout, plot_connectivity_circle

os.chdir('/media/dian/D/data/Binocular_Rivalry/HMM_analyses/output_test/embedded_lag/DMN_V1_1000_100_K4_pca16/run1_analyses')

for k in range(1,5):
    for freq in ['Theta', 'Alpha-1', 'Beta-1', 'Beta-2']:
        
        print('Plotting for state %d (%s)...' %(k, freq) )
        
        fn_keys = 'PCh_conn_K' + str(k) + '_' + freq
        fn = fn_keys + '_Tvalue.Rdata'
        
        ro.r['load'](fn)
        
        From = np.array(ro.r['from'])
        To = np.array(ro.r['to'])
        
        if ro.r['middle_point'] == ro.rinterface.NULL:
            middle_point= np.NaN
        else:
            middle_point = np.array(ro.r['middle_point'])
        
        if ro.r['conn_values'] == ro.rinterface.NULL: # if no significant connectivity
            print('Skip!')
            continue # then skip this iteration
       
        conn = np.array(ro.r['conn_values'])

        # To pandas DataFrame
        with localconverter(ro.default_converter + pandas2ri.converter):
            vertices = ro.conversion.rpy2py(ro.r['vertices'])
        
        node_names = list(vertices.name[3:])
        node_order = ['V1', 'PCU', 'PCC', 'lIPL', 'rIPL','lHP','rHP','ACC']
        
        # find V1
        v1_pos = [i for i,x in enumerate(node_order) if (x=='V1')]
        
        if v1_pos[0] == 7:
            v1_pos.append(0)
        else:
            v1_pos.append(v1_pos[0]+1)
        
        v1_pos.sort()
        
        node_angles = circular_layout(node_names, node_order, start_pos=90,
                                      group_boundaries=v1_pos)
        
        indices = (From-4, To-4) 
        
        # undirected connectivity
        con, ind = np.unique(conn, return_index=True)
        indices = (From[ind]-4, To[ind]-4)
        
        #color rgb
        n_colors = [None]*len(node_order)
        n = 0 # counter
        for i in range(len(node_order)):
            if node_order[i]=='V1':
                n_colors[i] = (0, 0, 0)
            else:
                n += 1
                high = 240
                low = 100
                c = (n*(high-low)/(len(node_order)-1)+low)/256
                dc = 15/256
                n_colors[i]=(c+dc*0.8,c+dc*0.4,c-dc*0.5)
                
        matches = [node_order.index(i) for i in node_names]
        node_colors = [n_colors[i] for i in matches]
        
        # colorbar limits
        #vmax = max(abs(con))
        #vmin = -vmax
        mdp = middle_point
        #if max(con)==mdp: # all decreased conn
         #   vmin = min(con)
         #   vmax = 2*mdp+vmin
       # else:
         #   if min(con)==mdp:
          #          vmax = max(con)
          #          vmin = vmax- 2*mdp
           # else:
            #    vmax = max(abs(con-mdp))*0.4 + mdp
             #   vmin = mdp-max(abs(con-mdp))*0.4
        
        
        if ro.r['note'] != ro.rinterface.NULL:
            note = ro.r['note'][0]
            if note =='all increased connectivity':
                vmin = 0
                vmax = max(con)*0.7
            elif note=='all decreased connectivity':
                vmin = min(con)
                vmax = vmin+mdp*2
        else:
             vmax = max(abs(con-mdp))*0.9 + mdp
             vmin = mdp-max(abs(con-mdp))*0.9   
                
  #----------------------------------------------------------      
# deal with extreme value in a special case: CONN between PCu and V1 is highest among all connections, but decreased across states
# the problem only exists in the 'pca=4' file
        if k==2 and freq=='Theta':
            con[np.where(con==np.max(con))]=mdp*0.5
            mdp = (max(con)-min(con))/2 + min(con)
            vmax = max(abs(con-mdp))*3 + mdp
            vmin = mdp-max(abs(con-mdp))*3
  #----------------------------------------------------------           
       # if not np.isnan(mdp):
         #   vmax = max(abs(con-mdp)) + mdp
         #   vmin = mdp-max(abs(con-mdp))
       # else:
        #    vmax = max(con) + 0.1*np.mean(con)
         #   vmin = min(con)-0.1*np.mean(con)

        
        # create rescale function
        def rescale(x, range2):
            range1 = [min(x), max(x)]
            delta1 = range1[1] - range1[0]
            delta2 = range2[1] - range2[0]
            return (delta2 * (x - range1[0]) / delta1) + range2[0]
        
        # title
        tlt = 'Phase-Coherence Connectivity among DMN regions and V1 \n for \n State %d (%s)' %(k, freq)
        
        fig = plt.figure(num=None, figsize=(8, 8))
        plot_connectivity_circle(con, node_names, indices=indices, node_angles = node_angles,
                                         facecolor='white', textcolor = 'black', node_colors = node_colors,
                                         linewidth = 3,
                                         node_edgecolor='black', colormap='RdBu_r',
                                         fontsize_names=12,
                                         fontsize_colorbar=12,
                                         vmin=vmin, vmax=vmax,
                                         colorbar_pos=(-0.7, 0.25),
                                         title = tlt, fig=fig)

        # save plot
        fig.savefig('fig/CirclePlt_'+fn_keys+'_Tvalue.png', dpi=300, transparent=True)
