# Module: Postprocessing and Plotter

import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import ScalarFormatter

# dr.Patočka:
#vpl=np.loadtxt("spada.dat")[:,11]
#vpt=np.loadtxt("spada.dat")[:,0]

# LOAD DATA (1: radius 2: density 3: mass 4: gravity 5: material parameters)
(rCenters,rhoCenters,mCenters,gCenters)=np.loadtxt('run/gridcenters.txt',dtype=float,unpack=True)
(rFaces,rhoFaces,mFaces,gFaces,muFaces,etaFaces,lmdFaces)=np.loadtxt('run/gridfaces.txt',dtype=float,unpack=True)
(time,dt,urtop,urtopIm,urbot,urbotIm,k2f,h2f,l2f,drpsg)=np.loadtxt('run/output_evol.dat',dtype=float,unpack=True)
KFaces=muFaces-(2./3.)*muFaces

# PLOT THE MODEL AND PARAMETER PROFILES:
fig,axs=plt.subplots(nrows=2,ncols=2,figsize=(10,8))
title="Parameters"
fig.subplots_adjust(hspace=0.40,top=0.90)
fig.suptitle(title,fontsize=16)
axs[0,0].scatter(rhoCenters[:-1],rCenters[:-1],s=1,marker='x',color='purple',label='Mantle')
axs[0,0].scatter(rhoCenters[-1],rCenters[-1],s=1,marker='x',color='red',label='Core')
axs[0,0].set_ylabel("Radius [m]",fontsize=12)    
axs[0,0].set_xlabel("Density [m]",fontsize=12)    
axs[0,0].legend(fontsize=12)
axs[0,1].scatter(muFaces,rFaces,s=1,marker='x',color='purple',label='Mantle')
axs[0,1].set_ylabel("Radius [m]",fontsize=12)    
axs[0,1].set_xlabel("Shear modulus [Pa]",fontsize=12)    
axs[0,1].legend(fontsize=12)
axs[1,0].scatter(KFaces,rFaces,s=1,marker='x',color='purple',label='Mantle')
axs[1,0].set_ylabel("Radius [m]",fontsize=12)    
axs[1,0].set_xlabel("Bulk Modulus [Pa]",fontsize=12)    
axs[1,0].legend(fontsize=12)
axs[1,1].scatter(gCenters,rCenters,s=1,marker='x',color='purple',label='Mantle')
axs[1,1].set_ylabel("Radius [m]",fontsize=12)    
axs[1,1].set_xlabel("Acceleration [ms$^{-2}$]",fontsize=12)    
axs[1,1].legend(fontsize=12)
fig.savefig('run/profiles.pdf',dpi=50)

fig,axs=plt.subplots(nrows=2,ncols=2,figsize=(16,10))
#ax.yaxis.set_major_formatter(ScalarFormatter())
#ax.minorticks_off()
label_k2e='$k_{2e}^T$='+str(k2f[0])
label_k2f="$k_{2f}^T$="+str(k2f[-1])
label_h2e='$h_{2e}^T$='+str(h2f[0])
label_h2f="$h_{2f}^T$="+str(h2f[-1])
label_l2e='$l_{2e}^T$='+str(l2f[0])
label_l2f="$l_{2f}^T$="+str(l2f[-1])
title="Time Evolution"
fig.subplots_adjust(hspace=0.40,top=0.90)
fig.suptitle(title,fontsize=18)
axs[0,0].plot(np.real(time[1:]),(-1)*np.real(urbot[1:]),color="black",linestyle="-",label='urbot')
axs[0,0].plot(np.real(time[1:]),(-1)*np.real(urtop[1:]),color="red",linestyle="-",label='urtop')
axs[0,0].plot(np.real(time[1:]),np.real(drpsg[1:]),color="orange",linestyle="-",label='(driving+sg)/$g_0$')
axs[0,0].set_xscale('log')
axs[0,0].set_title("Displacement",color="black",fontsize=14)
axs[0,0].set_ylabel("Dislacement [m]",fontsize=14)
axs[0,0].set_xlabel("Time [yr]",fontsize=14)
axs[0,0].legend(fontsize=14)
axs[0,0].margins(y=0.1,x=0)
##########################################################
axs[0,1].plot(np.real(time[1:]),np.real(k2f[1:]),color="black",label=label_k2f)
axs[0,1].set_title("Tidal Love Number $k_2$",color="black",fontsize=14)
#axs[0,1].plot(np.real(vpt[1:]),np.real(vpl[1:]),color="blue",label="dr.Patočka")
axs[0,1].axhline(y=k2f[0],color="black",linestyle=":",label=label_k2e)
axs[0,1].set_xscale('log')
axs[0,1].set_ylabel("Love number $k_2$",fontsize=14)
axs[0,1].set_xlabel("Time [yr]",fontsize=14)
axs[0,1].legend(loc="best",fontsize=14)
axs[0,1].margins(y=0.1,x=0)
##########################################################
axs[1,0].plot(np.real(time[1:]),np.real(h2f[1:]),color="black",label=label_h2f)
axs[1,0].set_title("Tidal Love Number $h_2$",color="black",fontsize=14)
#ax.plot(np.real(vpt[1:]),np.real(vpl[1:]),color="blue",label="dr.Patočka")
axs[1,0].axhline(y=h2f[0],color="black",linestyle=":",label=label_h2e)
axs[1,0].set_xscale('log')
axs[1,0].set_ylabel("Love number $h_2$",fontsize=14)
axs[1,0].set_xlabel("Time [yr]",fontsize=14)
axs[1,0].legend(loc="best",fontsize=14)
axs[1,0].margins(y=0.1,x=0)
##########################################################
axs[1,1].plot(np.real(time[1:]),np.real(l2f[1:]),color="black",label=label_l2f)
axs[1,1].set_title("Tidal Love Number $l_2$",color="black",fontsize=14)
#ax.plot(np.real(vpt[1:]),np.real(vpl[1:]),color="blue",label="dr.Patočka")
axs[1,1].axhline(y=l2f[0],color="black",linestyle=":",label=label_l2e)
axs[1,1].set_xscale('log')
axs[1,1].set_ylabel("Love number $l_2$",fontsize=14)
axs[1,1].set_xlabel("Time [yr]",fontsize=14)
axs[1,1].legend(loc="best",fontsize=14)
axs[1,1].margins(y=0.1,x=0)
plt.savefig('run/love2f.pdf',dpi=50)    
##########################################################
# end of the module!