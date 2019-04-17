# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 15:11:34 2019

@author: haessler, adapted from Matlab to Python by Marie
"""
#comments :
#Much shorter than the original Matlab script (which remains here as comments) but it seems that only plots parameters have been deleted (?) 
#Anyways it does return the same B and C values 


#import libraries
import numpy as np
import scipy as sc
from matplotlib import pyplot as plt


#define the function
def fit_calibration_stefan(nn,xn,deltaE,l,inc_angle,pxsize,gratinglines):

    hbar = 1.0545718e-34  #Planck constant 
    c = 299792458  # speed of light
    e = 1.6021766208e-19  #J = 1eV
    
    # Set up figure to receive data sets and fits
    plt.figure()
#    set(f_,'Units','Pixels','Position',[445 129 688 485]);
#    # Line handles and text for the legend.
#    legh_ = [];
#    legt_ = {};
#    # Limits of the x-axis.
#    xlim_ = [Inf -Inf];
#    # Axes for the plot.
#    ax_ = axes;
#    set(ax_,'Units','normalized','OuterPosition',[0 0 1 1]);
#    set(ax_,'Box','on');
#    axes(ax_);
#    hold on;
#    
#    # --- Plot data that was originally in data set "xn vs. nn"
    plt.plot(nn,xn, 'o', label = 'xn (position) vs nn (offseted harmonic number)')
#    nn = nn(:);
#    xn = xn(:);
#    h_ = line(nn,xn,'Parent',ax_,'Color',[0.333333 0 0.666667],...
#        'LineStyle','none', 'LineWidth',1,...
#        'Marker','.', 'MarkerSize',12);
#    xlim_(1) = min(xlim_(1),min(nn));
#    xlim_(2) = max(xlim_(2),max(nn));
#    legh_(end+1) = h_;
#    legt_{end+1} = 'xn vs. nn';
#    
#    # Nudge axis limits beyond data limits
#    if all(isfinite(xlim_))
#        xlim_ = xlim_ + [-1 1] * 0.01 * diff(xlim_);
#        set(ax_,'XLim',xlim_)
#    else
#        set(ax_, 'XLim',[-0.41000000000000003, 41.409999999999997]);
#    end
#    
    # --- Fit nn, xn to get C and B :
    
    const1 = hbar /e *2*np.pi*c *gratinglines    # hbar/e * 2*pi*c /d
    def expr (n, B, C):
        return -l/pxsize / np.tan( np.arcsin( const1/((n+C)*deltaE) - np.sin(inc_angle/180*np.pi) )  ) - B/pxsize
    [cf_,gof] = sc.optimize.curve_fit(expr, nn, xn)  #cf_ = [B, C]

    plt.plot( np.arange(nn[-1], nn[0], 0.01) , expr( np.arange(nn[-1], nn[0], 0.01), cf_[0] ,cf_[1] ) , color= 'red' , label = 'fit')
    plt.legend()
#    # Plot this fit
#    h_ = plot(cf_,'fit',0.95);
#    set(h_(1),'Color',[1 0 0],...
#        'LineStyle','-', 'LineWidth',2,...
#        'Marker','none', 'MarkerSize',6);
#    # Turn off legend created by plot method.
#    legend off;
#    # Store line handle and fit name for legend.
#    legh_(end+1) = h_(1);
#    legt_{end+1} = 'fit 1';
#    
#    # --- Finished fitting and plotting data. Clean up.
#    hold off;
#    # Display legend
#    leginfo_ = {'Orientation', 'vertical', 'Location', 'NorthEast'};
#    h_ = legend(ax_,legh_,legt_,leginfo_{:});
#    set(h_,'Interpreter','none');
#    # Remove labels from x- and y-axes.


    return [cf_, gof]