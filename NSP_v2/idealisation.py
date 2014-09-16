# -*- coding: utf-8 -*-
"""
Created on Wed May 28 15:23:24 2014

@author: chiaracasotto
"""
import numpy as np
import matplotlib.pyplot as plt

def bilinear(droof,Vb,flag,linew,fontsize,units):
#   FEMA method
    droof = np.array(droof)
    Fy = np.max(Vb)
    du = np.max(droof)
    dmax = droof[Vb == Fy]
    Ay = 0.6*Fy
    Ax = np.interp(Ay, Vb[droof<=dmax],droof[droof<=dmax])
    slp = Ay/Ax
    dy = Ax+ (Fy-Ay)/slp
    
    if flag:
        # Plot pushover curve and bilinear curve
        plt.plot(droof,Vb,color='b',label='pushover input',linewidth=linew)
        x = np.array([0, dy,du])
        y = np.array([0, Fy, Fy])
        plt.plot(x,y,color='r',marker = 'o',linewidth=linew,label='bilinear idealisation')
        plt.xlabel('roof displacement, droof '+units[0],fontsize = fontsize)
        plt.ylabel('base shear, Vb '+units[1],fontsize = fontsize)
        plt.suptitle('Pushover curve',fontsize = fontsize)
        plt.legend(loc='lower right',frameon = False)
        plt.show()
        
    return [dy,du,Fy]
    
def quadrilinear(droof,Vb,flag,linew,fontsize,units):
    droof = np.array(droof)
    Fmax = np.max(Vb)
    fmax = [i for i in range(0,len(Vb)) if Vb[i] == Fmax]
    fmax = fmax[0]
    dmax = droof[fmax]
 
#   Yielding point:
#   FEMA method    
#    By = np.linspace(0.6*Fmax,Fmax,100)
#    Ay = 0.6*By
#    Ax = np.interp(Ay, Vb[Vb<=Fmax],droof[Vb<=Fmax])
#    slope = Ay/Ax
#    Bx = By/slope
#    error = np.zeros_like(By)
#    for i in range(0,len(error)):
#        r1 = np.interp(droof[droof<=Bx[i]], [0, Bx[i]], [0, By[i]])
#        droofBx = droof[droof>Bx[i]]
#        r2 = np.interp(droof[droofBx<dmax], [Bx[i], dmax], [By[i], Fmax])
#        difference = np.concatenate((r1-Vb[droof<=Bx[i]], r2-Vb[droofBx<dmax]))
#        error[i] = np.abs(np.sum(difference))
#    
#    Fy = Vb[error == np.min(error)][0]
#    dy = droof[error == np.min(error)][0]
 
#   Yielding point:
#   Vulnerability guidelines method
#   Find yielding displacement with equal energy principle n the interval from 0 to Dmax
    Areas = np.array([(Vb[i+1]+Vb[i])/2 for i in range(0,fmax)])
    dd = np.array([droof[i+1]-droof[i] for i in range(0,fmax)])    
    Edmax = np.sum(dd*Areas) #Area under the pushover curve in the interval from 0 to Dmax         
    dy = 2*(dmax-Edmax/Fmax)
    Fy = Fmax
    
#   Onset of plateu
#   Find were slope of pushover curve before decreasing in the plateu
    Vb_norm = Vb/Fy
    d_norm = droof/dy
    slp = [(Vb_norm[i]-Vb_norm[i-1])/(d_norm[i]-d_norm[i-1]) for i in xrange(1,len(Vb))]   
    indy_soft = np.nonzero(abs(np.array(slp))>0.2)
    fmin = indy_soft[0][-1]
    Fmin = Vb[fmin]
    dmin = droof[fmin]
    
#   Onset of softening
#   Find yielding displacement with equal energy principle n the interval from Dmax to Dmin (onset of plateu)
    Areas = np.array([(Vb[i+1]+Vb[i])/2 for i in range(fmax,fmin)])
    dd = np.array([droof[i+1]-droof[i] for i in range(fmax,fmin)])
    Edmin = np.sum(dd*Areas)
    ds = 2/(Fmax-Fmin)*(Edmin - (dmin-dmax)*Fmax + 0.5*dmin*(Fmax-Fmin))
    du = np.max(droof)
    
#   Residual Plateu     
    Areas= np.array([(Vb[i+1]+Vb[i])/2 for i in range(fmin,len(Vb)-1)])
    dd = np.array([droof[i+1]-droof[i] for i in range(fmin,len(Vb)-1)])
    Edplat = np.sum(dd*Areas)
    Fres = Edplat/(droof[-1]-dmin)
    slp_soft = abs((Fmax-Fmin)/(ds-dmin))
    dmin = dmin+(Fmin-Fres)/slp_soft
    Fmin = Fres
    
    if flag:
        # Plot pushover curve and bilinear curve
        plt.plot(droof,Vb,color='b',linewidth=linew,label='pushover input')
        x = np.array([0, dy, ds, dmin, du])
        y = np.array([0, Fy, Fy, Fmin, Fmin])
        plt.plot(x,y,color='r',marker = 'o', linewidth=linew, label='quadrilinear idealisation')
        plt.xlabel('roof displacement, droof '+units[0],fontsize = fontsize)
        plt.ylabel('base shear, Vb '+units[1],fontsize = fontsize)
        plt.suptitle('Pushover curve',fontsize = fontsize)
        plt.legend(loc='lower right',frameon = False)
        plt.show()
    return [dy,ds,dmin,du,Fy,Fmax,Fmin]