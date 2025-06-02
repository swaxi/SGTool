"""
Plot functions 

A Python program to plot the total-field anomaly and the estimates
on classic plot. 

This code plot the figures 2d, 4 and 7 in the folder 'figures'.

This code is released from the paper: 
Reliable Euler deconvolution estimates throughout the
vertical derivatives of the total-field anomaly

The program is under the conditions terms in the file README.txt

authors: Felipe F. Melo and Valeria C.F. Barbosa, 2019
email: felipe146@hotmail.com, valcris@on.br
"""

import numpy as np
import matplotlib.pylab as plt
import matplotlib.patches as patches

#######################################################################

def plot_input_data(data,xi,yi,zi,shape):

    '''
    Plot the input data - Figure 2d
    '''
    
    fig=plt.figure(figsize=(5, 4))
    
    rect0 = patches.Rectangle((26,0),2,24,linewidth=1,edgecolor='black',
                              facecolor='none',linestyle='-',zorder=2)
    rect1 = patches.Rectangle((16.850,10),0.3,24,linewidth=1,
                              edgecolor='black',facecolor='none',
                              linestyle='-',zorder=2)
    rect2 = patches.Rectangle((6.850,15.850),0.3,0.3,linewidth=1,
                              edgecolor='black',facecolor='none',
                              linestyle='-',zorder=2)
    rect3 = patches.Rectangle((6.800,6.800),0.4,0.4,linewidth=1,
                              edgecolor='black',facecolor='none',
                              linestyle='-',zorder=2)

    im=plt.contourf(yi.reshape(shape)/1000.,xi.reshape(shape)/1000.,
                   data.reshape(shape), 30, cmap='jet')
    ax = plt.gca()
    ax.set_ylabel('Northing (km)', fontsize = 14)
    ax.set_xlabel('Easting (km)', fontsize = 14)
    ax.tick_params(labelsize=13)
    cbar=fig.colorbar(im,pad=0.01,shrink=1)
    cbar.set_label('nT',labelpad=-21,y=-0.03, rotation=0,fontsize=13)
    ax.add_patch(rect0)
    ax.add_patch(rect1)
    ax.add_patch(rect2)
    ax.add_patch(rect3)
    cbar.ax.tick_params(labelsize=13)
    ax.set_xlim(np.min(yi/1000.),np.max(yi/1000.))
    ax.set_ylim(np.min(xi/1000.),np.max(xi/1000.))
    plt.text(22.5,5,'P0',color='k', size='large')
    plt.text(13.5,14,'P1',color='k', size='large')
    plt.text(3,16,'P2',color='k', size='large')
    plt.text(3,5,'P3',color='k', size='large')
    
    plt.savefig('figures/FIG2d.png', bbox_inches='tight', dpi = 600)
    plt.close('all')
    
    return


def plot_classic(data,est_classic,xi,yi,zi,shape):
    '''
    Classic plot of the depth and base level estimates for all SIs
    '''
    
    '''
    Figure 4 source-position (depth) estimates
    '''
    vet_title=["(a)","(b)","(c)","(d)"]
    
    
    minz = np.min(0.)
    maxz = np.max(2.)
    levelsz = np.linspace(minz,maxz,11)   
    
    #depth plots
    plt.figure(figsize=(12, 8.5))
    for i in range (4):
     
        plt.subplot(2,2,i+1)  
        plt.title(vet_title[i], fontsize = 14, loc='center',y=-0.27)
        plt.contourf(yi.reshape(shape)/1000.,xi.reshape(shape)/1000.,
                       data.reshape(shape), 30, cmap='gray')
        ax = plt.gca()
        ax.set_ylabel('Northing (km)', fontsize = 14)
        ax.set_xlabel('Easting (km)', fontsize = 14)
        ax.tick_params(labelsize=13)
        scat=plt.scatter(est_classic[i][:,1]/1000., 
                         est_classic[i][:,0]/1000.,s=40,
                         c=(est_classic[i][:,2]/1000.),
                         cmap='terrain_r',vmin=minz,vmax=maxz,
                         edgecolors='k')             
        cbar=plt.colorbar(scat,ticks=levelsz,pad=0.01,shrink=1,
                          format='%0.1f')
        cbar.set_label('$\^z_o$ (km)',labelpad=-18,y=-0.03, rotation=0,
                       fontsize=13)
        cbar.ax.tick_params(labelsize=13)
        ax.set_xlim(np.min(yi/1000.),np.max(yi/1000.))
        ax.set_ylim(np.min(xi/1000.),np.max(xi/1000.))
        ax.set_xticks([0,5,10,15,20,25])
        plt.text(22.5,5,'P0',color='w', size='large')
        plt.text(13.5,14,'P1',color='w', size='large')
        plt.text(3,16,'P2',color='w', size='large')
        plt.text(3,5,'P3',color='w', size='large')
        
    plt.subplots_adjust(wspace=0.15,hspace=0.32)
        
    plt.savefig('figures/FIG4.png',bbox_inches='tight', dpi = 600)
    plt.close('all')
             
    '''
    Figure 7 base level-position (base-level) estimates
    '''
    plt.figure(figsize=(12, 8.5))
    for i in range (4):
     
        if i == 0:
            #base level estimates for SI = 0 have higher amplitude
            minb = np.min(-70)
            maxb = np.max(20)
            levelsb = np.linspace(minb,maxb,7)
            
            plt.subplot(2,2,i+1) 
            plt.title(vet_title[i], fontsize = 14, loc='center',
                      y=-0.27)
            plt.contourf(yi.reshape(shape)/1000.,xi.reshape(shape)/1000.,
                           data.reshape(shape), 30, cmap='gray')
            ax=plt.gca()
            ax.set_ylabel('Northing (km)', fontsize = 14)
            ax.set_xlabel('Easting (km)', fontsize = 14)
            ax.tick_params(labelsize=13)
            scat=plt.scatter(est_classic[i][:,1]/1000.,
                             est_classic[i][:,0]/1000.,s=40,
                             c=(est_classic[i][:,3]/1000.), 
                             cmap='jet',vmin=minb,
                             vmax=maxb,edgecolors='k')             
            cbar=plt.colorbar(scat,ticks=levelsb,pad=0.01,shrink=1,
                              format='%d')
            cbar.set_label('$\^b$ (nT)',labelpad=-18,y=-0.03, rotation=0,
                           fontsize=13)
            cbar.ax.tick_params(labelsize=13)
            ax.set_xlim(np.min(yi/1000.),np.max(yi/1000.))
            ax.set_ylim(np.min(xi/1000.),np.max(xi/1000.))
            ax.set_xticks([0,5,10,15,20,25])
            plt.text(22.5,5,'P0',color='w', size='large')
            plt.text(13.5,14,'P1',color='w', size='large')
            plt.text(3,16,'P2',color='w', size='large')
            plt.text(3,5,'P3',color='w', size='large') 
            plt.text(25,25,'x10$^{3}$',color='k',size='medium')
            
        else:
            minb = np.min(-30)
            maxb = np.max(210)
            levelsb = np.linspace(minb,maxb,7)
      
            plt.subplot(2,2,i+1) 
            plt.title(vet_title[i], fontsize = 14, loc='center',
                      y=-0.27)
            plt.contourf(yi.reshape(shape)/1000.,xi.reshape(shape)/1000.,
                           data.reshape(shape), 30, cmap='gray')
            ax=plt.gca()
            ax.set_ylabel('Northing (km)', fontsize = 14)
            ax.set_xlabel('Easting (km)', fontsize = 14)
            ax.tick_params(labelsize=13)
            scat=plt.scatter(est_classic[i][:,1]/1000., 
                             est_classic[i][:,0]/1000.,s=40,
                             c=(est_classic[i][:,3]),cmap='jet',
                             vmin=minb, vmax=maxb,edgecolors='k')        
            cbar=plt.colorbar(scat,ticks=levelsb,pad=0.01,shrink=1,
                              format='%d')
            cbar.set_label('$\^b$ (nT)',labelpad=-18,y=-0.03, rotation=0,
                           fontsize=13)
            cbar.ax.tick_params(labelsize=13)
            ax.set_xlim(np.min(yi/1000.),np.max(yi/1000.))
            ax.set_ylim(np.min(xi/1000.),np.max(xi/1000.))
            ax.set_xticks([0,5,10,15,20,25])
            plt.text(22.5,5,'P0',color='w', size='large')
            plt.text(13.5,14,'P1',color='w', size='large')
            plt.text(3,16,'P2',color='w', size='large')
            plt.text(3,5,'P3',color='w', size='large')   
            
    plt.subplots_adjust(wspace=0.15,hspace=0.32)
        
    plt.savefig('figures/FIG7.png',bbox_inches='tight', dpi = 600)
    plt.close('all')
    
    return