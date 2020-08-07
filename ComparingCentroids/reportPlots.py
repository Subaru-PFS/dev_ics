"""

Collection of routines which were used to make plots for the report. 

The routines mostly read in the data cubes created earlier, or the
output of the differences routine, does some basic math on them, and
plots the results in a pretty fashion.

Figure numbers are from version 4 of the report. 

For reference, the tags are as follows
  w: result of windowed centroiding (sextrator)
  p: results of psf fitting (sextractor)
  l: current implmentation of MCS code
  t: simple thresholding
  c: simulated data

The two sets of data sets chosen
    5854,15864,15874,15884
    762200,1762300,1762400,1762500
were selected as typical, unremarkable results of rotation and exposure time, respectively. 
The offset variable selects a specific frame in the sets.


"""


import numpy as np
import matplotlib.pyplot as plt
import analysisRoutinesNew as ar
import numpy.ma as ma
import matplotlib as mlab

def offsets():

    """ 

    Plot offsets between different methods. 
    Figures 42, 43

"""

    #read in the output produced by differences
    
    ff1=np.loadtxt("w_p_diffs1.dat")
    ff2=np.loadtxt("w_l_diffs1.dat")
    ff3=np.loadtxt("l_p_diffs1.dat")
    
    fig,ax=plt.subplots(3,2,figsize=(14,10))

    #pick rotation=1, to avoid axis flipping problems
    ind=np.where(ff1[:,3]==0)

    #plot histograms of the differences. The - signs keep things oriented the same way 
    
    ax[0,0].hist(ff1[ind,4].ravel(),10)
    ax[0,0].set_title("Mean dx (Windowed - PSF Fitting)")
    ax[1,0].hist(-ff2[ind,4].ravel(),10)
    ax[1,0].set_title("Mean dx (MCS-Windowed)")
    ax[2,0].hist(ff3[ind,4].ravel(),10)
    ax[2,0].set_title("Mean dx (MCS - PSF Fitting)")
    ax[0,1].hist(ff1[ind,5].ravel(),10)
    ax[0,1].set_title("Mean dy (Windowed - PSF Fitting)")
    ax[1,1].hist(-ff2[ind,5].ravel(),10)
    ax[1,1].set_title("Mean dy (MCS-Windowed)")
    ax[2,1].hist(ff3[ind,5].ravel(),10)
    ax[2,1].set_title("Mean dy (MCS - PSF Fitting)")

    #print(np.median(ff1[ind,4]),np.median(ff2[ind,4]),np.median(ff3[ind,4]))
    #print(np.median(ff1[ind,5]),np.median(ff2[ind,5]),np.median(ff3[ind,5]))

    #set a vertical line at zero
    
    ax[0,0].axvline(x=0,color="k")
    ax[1,0].axvline(x=0,color="k")
    ax[2,0].axvline(x=0,color="k")
    ax[0,1].axvline(x=0,color="k")
    ax[1,1].axvline(x=0,color="k")
    ax[2,1].axvline(x=0,color="k")

    #some labels
    ax[2,0].set_xlabel("dx (microns)")
    ax[2,1].set_xlabel("dy (microns)")

    ax[0,0].set_ylabel("N")
    ax[1,0].set_ylabel("N")
    ax[2,0].set_ylabel("N")
    
    plt.savefig("diffHist.png")


    #plot histograms of the rms. The - signs keep things oriented the same way 
    
    fig,ax=plt.subplots(3,2,figsize=(14,10))

    ind=np.where(ff1[:,3]==0)
    ax[0,0].hist(ff1[ind,10].ravel(),10)
    ax[0,0].set_title("RMS dx (Windowed - PSF Fitting)")
    ax[1,0].hist(ff2[ind,10].ravel(),10)
    ax[1,0].set_title("RMS dx (MCS-Windowed)")
    ax[2,0].hist(ff3[ind,10].ravel(),10)
    ax[2,0].set_title("RMS dx (MCS - PSF Fitting)")
    ax[0,1].hist(ff1[ind,11].ravel(),10)
    ax[0,1].set_title("RMS dy (Windowed - PSF Fitting)")
    ax[1,1].hist(ff2[ind,11].ravel(),10)
    ax[1,1].set_title("RMS dy (MCS-Windowed)")
    ax[2,1].hist(ff3[ind,11].ravel(),10)
    ax[2,1].set_title("RMS dy (MCS - PSF Fitting)")

    #print(np.median(ff1[ind,4]),np.median(ff2[ind,4]),np.median(ff3[ind,4]))
    #print(np.median(ff1[ind,5]),np.median(ff2[ind,5]),np.median(ff3[ind,5]))
    
    ax[2,0].set_xlabel("dx (microns)")
    ax[2,1].set_xlabel("dy (microns)")

    ax[0,0].set_ylabel("N")
    ax[1,0].set_ylabel("N")
    ax[2,0].set_ylabel("N")
    
    plt.savefig("rmsHist.png")

    #a scatter plot of results to check the averages for different orientations and elevations
    
    fig,ax=plt.subplots()
    for i in range(len(ff3[:,0].ravel())):
        if(ff3[i,1]<=0.5):
            col="red"
        elif(ff3[i,1]==1):
            col="orange"
        elif(ff3[i,1]==2):
            col="yellow"
        elif(ff3[i,1]>=5):
            col="green"
        if(ff3[i,2]==60):
            sym="x"
        elif(ff3[i,2]==45):
            sym="s"
        elif(ff3[i,2]==30):
            sym="d"
        elif(ff3[i,2]==15):
            sym="^"
        elif(ff3[i,2]==0):
            sym="o"
        ax.plot([ff3[i,7]],ff3[i,4],color=col,marker=sym)
    plt.savefig("test.png")

    fig,ax=plt.subplots()
    for i in range(len(ff3[:,0].ravel())):
        if(ff3[i,3]==-180):
            col="red"
        if(ff3[i,3]==-90):
            col="orange"
        if(ff3[i,3]==0):
            col="yellow"
        if(ff3[i,3]==90):
            col="green"

        ax.plot(ff3[i,7],ff3[i,4],color=col,marker="d")


    plt.savefig("test1.png")

def simVals():

    """ 
    Compare mean RMS of different centroiding methods as a function of exposure time
    Figure 49
    """

    #load the saved data
    ff=np.loadtxt("rmsMic.dat")

    fig,ax=plt.subplots()

    ind1=np.where(ff[:,0]==1)
    ind2=np.where(ff[:,0]==2)
    ind3=np.where(ff[:,0]==3)
    ind4=np.where(ff[:,0]==4)
    ind5=np.where(ff[:,0]==5)
    ind6=np.where(ff[:,0]==6)

    ax.plot(ff[ind5,1].ravel(),ff[ind5,5].ravel(),marker="d",label="PSF")
    ax.plot(ff[ind4,1].ravel(),ff[ind4,5].ravel(),marker="d",label="Win")
    ax.plot(ff[ind6,1].ravel(),ff[ind6,5].ravel(),marker="d",label="MCS")

    ax.set_xlabel("Exposure time (s)")
    ax.set_ylabel("Mean RMS (microns)")
    plt.title("Comparison of RMS of different algorithms")
    plt.legend()
    plt.savefig("rmsVar.png")

def simStats(frameID):

    """ 
    calculate differences between methods for a single frame. 
    This is just loading the cubes and subtracting 
    Outputs numbers
    """
    
    loadPref="simSet/"
    source="1"

    cubeX1,cubeY1,cubeFX1,cubeFY1,cubeP1,cubeFXs1,cubeFYs1,cubePs1,cubeXD1,cubeYD1,cubeXDs1,cubeYDs1,xAv1,yAv1,fxAv1,fxAv1,fyAv1,pAv1,fxAvs1,fyAvs1,pAvs1,rmsVal1,rmsX1,rmsY1,rmsVals1,rmsXs1,rmsYs1=ar.loadCubes(loadPref,source,frameID)

    source="2"

    cubeX2,cubeY2,cubeFX2,cubeFY2,cubeP2,cubeFXs2,cubeFYs2,cubePs2,cubeXD2,cubeYD2,cubeXDs2,cubeYDs2,xAv2,yAv2,fxAv2,fxAv2,fyAv2,pAv2,fxAvs2,fyAvs2,pAvs2,rmsVal2,rmsX2,rmsY2,rmsVals2,rmsXs2,rmsYs2=ar.loadCubes(loadPref,source,frameID)

    source="3"

    cubeX3,cubeY3,cubeFX3,cubeFY3,cubeP3,cubeFXs3,cubeFYs3,cubePs3,cubeXD3,cubeYD3,cubeXDs3,cubeYDs3,xAv3,yAv3,fxAv3,fxAv3,fyAv3,pAv3,fxAvs3,fyAvs3,pAvs3,rmsVal3,rmsX3,rmsY3,rmsVals3,rmsXs3,rmsYs3=ar.loadCubes(loadPref,source,frameID)


    source="4"

    cubeX4,cubeY4,cubeFX4,cubeFY4,cubeP4,cubeFXs4,cubeFYs4,cubePs4,cubeXD4,cubeYD4,cubeXDs4,cubeYDs4,xAv4,yAv4,fxAv4,fxAv4,fyAv4,pAv4,fxAvs4,fyAvs4,pAvs4,rmsVal4,rmsX4,rmsY4,rmsVals4,rmsXs4,rmsYs4=ar.loadCubes(loadPref,source,frameID)


    source="5"

    cubeX5,cubeY5,cubeFX5,cubeFY5,cubeP5,cubeFXs5,cubeFYs5,cubePs5,cubeXD5,cubeYD5,cubeXDs5,cubeYDs5,xAv5,yAv5,fxAv5,fxAv5,fyAv5,pAv5,fxAvs5,fyAvs5,pAvs5,rmsVal5,rmsX5,rmsY5,rmsVals5,rmsXs5,rmsYs5=ar.loadCubes(loadPref,source,frameID)

    source="6"

    cubeX6,cubeY6,cubeFX6,cubeFY6,cubeP6,cubeFXs6,cubeFYs6,cubePs6,cubeXD6,cubeYD6,cubeXDs6,cubeYDs6,xAv6,yAv6,fxAv6,fxAv6,fyAv6,pAv6,fxAvs6,fyAvs6,pAvs6,rmsVal6,rmsX6,rmsY6,rmsVals6,rmsXs6,rmsYs6=ar.loadCubes(loadPref,source,frameID)

    print("%8d -1 %6.3f %6.3f" % (frameID,1000*(xAv1-xAv2).mean(),1000*(yAv1-yAv2).mean()))
    print("%8d -2 %6.3f %6.3f" % (frameID,1000*(xAv1-xAv3).mean(),1000*(yAv1-yAv3).mean()))
    print("%8d -3 %6.3f %6.3f" % (frameID,1000*(xAv2-xAv3).mean(),1000*(yAv2-yAv3).mean()))
    print("%8d -4 %6.3f %6.3f" % (frameID,1000*(xAv4-xAv5).mean(),1000*(yAv4-yAv5).mean()))
    print("%8d -5 %6.3f %6.3f" % (frameID,1000*(xAv4-xAv6).mean(),1000*(yAv4-yAv6).mean()))
    print("%8d -6 %6.3f %6.3f" % (frameID,1000*(xAv5-xAv6).mean(),1000*(yAv5-yAv6).mean()))

    
def plotSimStats():

    """  
    Compare systematic differences between methods using simulations. 
    Loads the output of simSats
    Figure 50
    """
    

    ff=np.loadtxt("deltaMic.dat")

    ind1=np.where(ff[:,2]==-1)
    ind2=np.where(ff[:,2]==-2)
    ind3=np.where(ff[:,2]==-3)
    ind4=np.where(ff[:,2]==-4)
    ind5=np.where(ff[:,2]==-5)
    ind6=np.where(ff[:,2]==-6)
    fig,ax=plt.subplots(2,1)

    ax[0].plot(ff[ind1,1].ravel(),ff[ind1,3].ravel(),marker="d",linestyle="-",label="Win-PSF (x5)")
    ax[0].plot(ff[ind2,1].ravel(),ff[ind2,3].ravel(),marker="d",linestyle="-",label="Win-MCS (x5)")
    ax[0].plot(ff[ind3,1].ravel(),ff[ind3,3].ravel(),marker="d",linestyle="-",label="PSF-MCS (x5)")
    ax[0].plot(ff[ind4,1].ravel(),ff[ind4,3].ravel(),marker="d",linestyle="-",label="Win-PSF")
    ax[0].plot(ff[ind5,1].ravel(),ff[ind5,3].ravel(),marker="d",linestyle="-",label="Win-MCS")
    ax[0].plot(ff[ind6,1].ravel(),ff[ind6,3].ravel(),marker="d",linestyle="-",label="PSF-MCS")
    ax[0].set_title("X Difference between methods")
    ax[0].set_xlabel("t")
    ax[0].set_ylabel("dx (microns)")

    ax[0].legend(loc=(0,.3),prop={'size': 8})

    ax[1].plot(ff[ind1,1].ravel(),ff[ind1,4].ravel(),marker="d",linestyle="-",label="Win-PSF (x5)")
    ax[1].plot(ff[ind2,1].ravel(),ff[ind2,4].ravel(),marker="d",linestyle="-",label="Win-MCS (x5)")
    ax[1].plot(ff[ind3,1].ravel(),ff[ind3,4].ravel(),marker="d",linestyle="-",label="PSF-MCS (x5)")
    ax[1].plot(ff[ind4,1].ravel(),ff[ind4,4].ravel(),marker="d",linestyle="-",label="Win-PSF")
    ax[1].plot(ff[ind5,1].ravel(),ff[ind5,4].ravel(),marker="d",linestyle="-",label="Win-MCS")
    ax[1].plot(ff[ind6,1].ravel(),ff[ind6,4].ravel(),marker="d",linestyle="-",label="PSF-MCS")
    ax[1].set_title("Y Difference between methods")
    ax[1].set_xlabel("t")
    ax[1].set_ylabel("dy (microns)")

    plt.savefig("deltaSim.png")

    #plt.legend()
    plt.savefig("dY.png")
        
    
def modelMoire():


    """  
    Look for moire pattern in simlated data
    Figure 48
    """
    
    source="c"
    
    loadPref="dbSet/"
    frameID=1762500
    
    fig,ax=plt.subplots(2,2,figsize=(14,10))
    st=ar.loadPhys("Phys/nphys_",frameID)
    cubeX,cubeY,cubeFX,cubeFY,cubeP,cubeFXs,cubeFYs,cubePs,cubeXD,cubeYD,cubeXDs,cubeYDs,xAv,yAv,fxAv,fxAv,fyAv,pAv,fxAvs,fyAvs,pAvs,rmsVal,rmsX,rmsY,rmsVals,rmsXs,rmsYs=ar.ar.loadCubes(loadPref,source,frameID)

    fig,ax=plt.subplots(2,2,figsize=(14,10))

    sc1=ax[0,0].scatter(xAv,yAv,c=fxAvs)
    ax[0,0].set_title("average spot size (x)")
    fig.colorbar(sc1,ax=ax[0,0])

    sc1=ax[0,1].scatter(xAv,yAv,c=fyAvs)
    ax[0,1].set_title("average spot size (y)")
    fig.colorbar(sc1,ax=ax[0,1])

    sc1=ax[1,0].scatter(xAv,yAv,c=cubeXDs[:,0])
    ax[1,0].set_title("dx")
    fig.colorbar(sc1,ax=ax[1,0])

    sc1=ax[1,1].scatter(xAv,yAv,c=pAv)
    ax[1,1].set_title("peak value")
    fig.colorbar(sc1,ax=ax[1,1])

    plt.suptitle("Patterns in Modelled Data ["+str(frameID)+"]}")
    plt.savefig("modelMoire.png")

def rotAvParm():

    """

    Plot averages of parameters (brightness, size)  as a function of rotation 
    for the specific data sets. 

    Figures 1-3

    """
    
    frames=np.array([15854,15864,15874,15884]).astype('int')
    offset=[4,1,5,0]

    loadPref="dbSet"
    source="l"

    fig,ax=plt.subplots(2,2,figsize=(14,10))
    frameID=frames[0]
    st=ar.loadPhys("Phys/nphys_",frameID)
    cubeX,cubeY,cubeFX,cubeFY,cubeP,cubeFXs,cubeFYs,cubePs,cubeXD,cubeYD,cubeXDs,cubeYDs,xAv,yAv,fxAv,fxAv,fyAv,pAv,fxAvs,fyAvs,pAvs,rmsVal,rmsX,rmsY,rmsVals,rmsXs,rmsYs=ar.loadCubes(loadPref,source,frameID)
    
    vmi,vma=ar.bounds(pAv,2)

    sc1=ax[0,0].scatter(xAv,yAv,c=pAv,vmin=vmi,vmax=vma)
    ax[0,0].set_title(str(int(frameID))+" ["+st+"]")
    fig.colorbar(sc1,ax=ax[0,0])

    frameID=frames[1]
    st=ar.loadPhys("Phys/nphys_",frameID)
    cubeX,cubeY,cubeFX,cubeFY,cubeP,cubeFXs,cubeFYs,cubePs,cubeXD,cubeYD,cubeXDs,cubeYDs,xAv,yAv,fxAv,fxAv,fyAv,pAv,fxAvs,fyAvs,pAvs,rmsVal,rmsX,rmsY,rmsVals,rmsXs,rmsYs=ar.loadCubes(loadPref,source,frameID)
    vmi,vma=ar.bounds(pAv,2)
    sc1=ax[0,1].scatter(xAv,yAv,c=pAv,vmin=vmi,vmax=vma)
    ax[0,1].set_title(str(int(frameID))+" ["+st+"]")
    fig.colorbar(sc1,ax=ax[0,1])

    frameID=frames[2]
    st=ar.loadPhys("Phys/nphys_",frameID)
    cubeX,cubeY,cubeFX,cubeFY,cubeP,cubeFXs,cubeFYs,cubePs,cubeXD,cubeYD,cubeXDs,cubeYDs,xAv,yAv,fxAv,fxAv,fyAv,pAv,fxAvs,fyAvs,pAvs,rmsVal,rmsX,rmsY,rmsVals,rmsXs,rmsYs=ar.loadCubes(loadPref,source,frameID)
    vmi,vma=ar.bounds(pAv,2)
    sc1=ax[1,0].scatter(xAv,yAv,c=pAv,vmin=vmi,vmax=vma)
    ax[1,0].set_title(str(int(frameID))+" ["+st+"]")
    fig.colorbar(sc1,ax=ax[1,0])

    frameID=frames[3]
    st=ar.loadPhys("Phys/nphys_",frameID)
    cubeX,cubeY,cubeFX,cubeFY,cubeP,cubeFXs,cubeFYs,cubePs,cubeXD,cubeYD,cubeXDs,cubeYDs,xAv,yAv,fxAv,fxAv,fyAv,pAv,fxAvs,fyAvs,pAvs,rmsVal,rmsX,rmsY,rmsVals,rmsXs,rmsYs=ar.loadCubes(loadPref,source,frameID)
    vmi,vma=ar.bounds(pAv,2)
    sc1=ax[1,1].scatter(xAv,yAv,c=pAv,vmin=vmi,vmax=vma)
    ax[1,1].set_title(str(int(frameID))+" ["+st+"]")
    fig.colorbar(sc1,ax=ax[1,1])

    plt.suptitle("Average Brightness with Rotation")
    plt.savefig("peakAvRot.png")

    #----------------------------------------------------------------------
    source="w"
    loadPref="dbSet/"
    
    fig,ax=plt.subplots(2,2,figsize=(14,10))
    frameID=frames[0]
    st=ar.loadPhys("Phys/nphys_",frameID)
    cubeX,cubeY,cubeFX,cubeFY,cubeP,cubeFXs,cubeFYs,cubePs,cubeXD,cubeYD,cubeXDs,cubeYDs,xAv,yAv,fxAv,fxAv,fyAv,pAv,fxAvs,fyAvs,pAvs,rmsVal,rmsX,rmsY,rmsVals,rmsXs,rmsYs=ar.loadCubes(loadPref,source,frameID)
    vmi,vma=ar.bounds(fxAv,2)
    sc1=ax[0,0].scatter(xAv,yAv,c=fxAv,vmin=vmi,vmax=vma)
    ax[0,0].set_title(str(int(frameID))+" ["+st+"]")
    fig.colorbar(sc1,ax=ax[0,0])

    frameID=frames[1]
    st=ar.loadPhys("Phys/nphys_",frameID)
    cubeX,cubeY,cubeFX,cubeFY,cubeP,cubeFXs,cubeFYs,cubePs,cubeXD,cubeYD,cubeXDs,cubeYDs,xAv,yAv,fxAv,fxAv,fyAv,pAv,fxAvs,fyAvs,pAvs,rmsVal,rmsX,rmsY,rmsVals,rmsXs,rmsYs=ar.loadCubes(loadPref,source,frameID)
    vmi,vma=ar.bounds(fxAv,2)
    sc1=ax[0,1].scatter(xAv,yAv,c=fxAv,vmin=vmi,vmax=vma)
    ax[0,1].set_title(str(int(frameID))+" ["+st+"]")
    fig.colorbar(sc1,ax=ax[0,1])

    frameID=frames[2]
    st=ar.loadPhys("Phys/nphys_",frameID)
    cubeX,cubeY,cubeFX,cubeFY,cubeP,cubeFXs,cubeFYs,cubePs,cubeXD,cubeYD,cubeXDs,cubeYDs,xAv,yAv,fxAv,fxAv,fyAv,pAv,fxAvs,fyAvs,pAvs,rmsVal,rmsX,rmsY,rmsVals,rmsXs,rmsYs=ar.loadCubes(loadPref,source,frameID)
    vmi,vma=ar.bounds(fxAv,2)
    sc1=ax[1,0].scatter(xAv,yAv,c=fxAv,vmin=vmi,vmax=vma)
    ax[1,0].set_title(str(int(frameID))+" ["+st+"]")
    fig.colorbar(sc1,ax=ax[1,0])

    frameID=frames[3]
    st=ar.loadPhys("Phys/nphys_",frameID)
    cubeX,cubeY,cubeFX,cubeFY,cubeP,cubeFXs,cubeFYs,cubePs,cubeXD,cubeYD,cubeXDs,cubeYDs,xAv,yAv,fxAv,fxAv,fyAv,pAv,fxAvs,fyAvs,pAvs,rmsVal,rmsX,rmsY,rmsVals,rmsXs,rmsYs=ar.loadCubes(loadPref,source,frameID)
    vmi,vma=ar.bounds(fxAv,2)
    sc1=ax[1,1].scatter(xAv,yAv,c=fxAv,vmin=vmi,vmax=vma)
    ax[1,1].set_title(str(int(frameID))+" ["+st+"]")
    fig.colorbar(sc1,ax=ax[1,1])

    plt.suptitle("Average Spot Size (x) with Rotation")
    plt.savefig("fxAvRot.png")

    #----------------------------------------------------------------------

    fig,ax=plt.subplots(2,2,figsize=(14,10))
    frameID=frames[0]
    st=ar.loadPhys("Phys/nphys_",frameID)
    cubeX,cubeY,cubeFX,cubeFY,cubeP,cubeFXs,cubeFYs,cubePs,cubeXD,cubeYD,cubeXDs,cubeYDs,xAv,yAv,fxAv,fxAv,fyAv,pAv,fxAvs,fyAvs,pAvs,rmsVal,rmsX,rmsY,rmsVals,rmsXs,rmsYs=ar.loadCubes(loadPref,source,frameID)
    vmi,vma=ar.bounds(fyAv,2)
    sc1=ax[0,0].scatter(xAv,yAv,c=fyAv,vmin=vmi,vmax=vma)
    ax[0,0].set_title(str(int(frameID))+" ["+st+"]")
    fig.colorbar(sc1,ax=ax[0,0])

    frameID=frames[1]
    st=ar.loadPhys("Phys/nphys_",frameID)
    cubeX,cubeY,cubeFX,cubeFY,cubeP,cubeFXs,cubeFYs,cubePs,cubeXD,cubeYD,cubeXDs,cubeYDs,xAv,yAv,fxAv,fxAv,fyAv,pAv,fxAvs,fyAvs,pAvs,rmsVal,rmsX,rmsY,rmsVals,rmsXs,rmsYs=ar.loadCubes(loadPref,source,frameID)
    vmi,vma=ar.bounds(fyAv,2)
    sc1=ax[0,1].scatter(xAv,yAv,c=fyAv,vmin=vmi,vmax=vma)
    ax[0,1].set_title(str(int(frameID))+" ["+st+"]")
    fig.colorbar(sc1,ax=ax[0,1])

    frameID=frames[2]
    st=ar.loadPhys("Phys/nphys_",frameID)
    cubeX,cubeY,cubeFX,cubeFY,cubeP,cubeFXs,cubeFYs,cubePs,cubeXD,cubeYD,cubeXDs,cubeYDs,xAv,yAv,fxAv,fxAv,fyAv,pAv,fxAvs,fyAvs,pAvs,rmsVal,rmsX,rmsY,rmsVals,rmsXs,rmsYs=ar.loadCubes(loadPref,source,frameID)
    vmi,vma=ar.bounds(fyAv,2)
    sc1=ax[1,0].scatter(xAv,yAv,c=fyAv,vmin=vmi,vmax=vma)
    ax[1,0].set_title(str(int(frameID))+" ["+st+"]")
    fig.colorbar(sc1,ax=ax[1,0])

    frameID=frames[3]
    st=ar.loadPhys("Phys/nphys_",frameID)
    cubeX,cubeY,cubeFX,cubeFY,cubeP,cubeFXs,cubeFYs,cubePs,cubeXD,cubeYD,cubeXDs,cubeYDs,xAv,yAv,fxAv,fxAv,fyAv,pAv,fxAvs,fyAvs,pAvs,rmsVal,rmsX,rmsY,rmsVals,rmsXs,rmsYs=ar.loadCubes(loadPref,source,frameID)
    vmi,vma=ar.bounds(fyAv,2)
    sc1=ax[1,1].scatter(xAv,yAv,c=fyAv,vmin=vmi,vmax=vma)
    ax[1,1].set_title(str(int(frameID))+" ["+st+"]")
    fig.colorbar(sc1,ax=ax[1,1])

    plt.suptitle("Average Spot Size (y) with Rotation")
    plt.savefig("fyAvRot.png")
   
    #----------------------------------------------------------------------

def rotDeltaS():

    """
    
    Plot rotation examples - detla parameters, with the smooth component subtracted.  
    dx,dy,fxs,fys

    Figures  24, 25, 27, 26, 5


    """

    
    frames=np.array([15854,15864,15874,15884]).astype('int')
    offset=[4,1,5,0]
    
    loadPref="dbSet"
    source="l"

    fig,ax=plt.subplots(2,2,figsize=(14,10))
    frameID=frames[0]
    off=offset[0]
    st=ar.loadPhys("Phys/nphys_",frameID)

    cubeX,cubeY,cubeFX,cubeFY,cubeP,cubeFXs,cubeFYs,cubePs,cubeXD,cubeYD,cubeXDs,cubeYDs,xAv,yAv,fxAv,fxAv,fyAv,pAv,fxAvs,fyAvs,pAvs,rmsVal,rmsX,rmsY,rmsVals,rmsXs,rmsYs=ar.loadCubes(loadPref,source,frameID)
    xmi,xma=ar.bounds(cubeXDs[:,off],2)
    sc1=ax[0,0].scatter(xAv,yAv,c=cubeXDs[:,off],vmin=xmi,vmax=xma)
    ax[0,0].set_title(str(int(frameID))+" ["+st+"]")
    fig.colorbar(sc1,ax=ax[0,0])

    frameID=frames[1]
    off=offset[1]
    st=ar.loadPhys("Phys/nphys_",frameID)
    cubeX,cubeY,cubeFX,cubeFY,cubeP,cubeFXs,cubeFYs,cubePs,cubeXD,cubeYD,cubeXDs,cubeYDs,xAv,yAv,fxAv,fxAv,fyAv,pAv,fxAvs,fyAvs,pAvs,rmsVal,rmsX,rmsY,rmsVals,rmsXs,rmsYs=ar.loadCubes(loadPref,source,frameID)
    xmi,xma=ar.bounds(cubeXDs[:,off],2)
    sc1=ax[0,1].scatter(xAv,yAv,c=cubeXDs[:,off],vmin=xmi,vmax=xma)
    ax[0,1].set_title(str(int(frameID))+" ["+st+"]")
    fig.colorbar(sc1,ax=ax[0,1])

    frameID=frames[2]
    off=offset[2]
    st=ar.loadPhys("Phys/nphys_",frameID)
    cubeX,cubeY,cubeFX,cubeFY,cubeP,cubeFXs,cubeFYs,cubePs,cubeXD,cubeYD,cubeXDs,cubeYDs,xAv,yAv,fxAv,fxAv,fyAv,pAv,fxAvs,fyAvs,pAvs,rmsVal,rmsX,rmsY,rmsVals,rmsXs,rmsYs=ar.loadCubes(loadPref,source,frameID)
    xmi,xma=ar.bounds(cubeXDs[:,off],2)
    sc1=ax[1,0].scatter(xAv,yAv,c=cubeXDs[:,off],vmin=xmi,vmax=xma)
    ax[1,0].set_title(str(int(frameID))+" ["+st+"]")
    fig.colorbar(sc1,ax=ax[1,0])

    frameID=frames[3]
    off=offset[3]
    st=ar.loadPhys("Phys/nphys_",frameID)
    cubeX,cubeY,cubeFX,cubeFY,cubeP,cubeFXs,cubeFYs,cubePs,cubeXD,cubeYD,cubeXDs,cubeYDs,xAv,yAv,fxAv,fxAv,fyAv,pAv,fxAvs,fyAvs,pAvs,rmsVal,rmsX,rmsY,rmsVals,rmsXs,rmsYs=ar.loadCubes(loadPref,source,frameID)
    xmi,xma=ar.bounds(cubeXDs[:,off],2)
    sc1=ax[1,1].scatter(xAv,yAv,c=cubeXDs[:,off],vmin=xmi,vmax=xma)
    ax[1,1].set_title(str(int(frameID))+" ["+st+"]")
    fig.colorbar(sc1,ax=ax[1,1])

    plt.suptitle("dx with Rotation (smooth component subtracted)")
    plt.savefig("dxRotS.png")

    #----------------------------------------------------------------------

    frames=np.array([15854,15864,15874,15884]).astype('int')
    offset=[4,1,5,0]

    loadPref="dbSet"
    source="l"

    fig,ax=plt.subplots(2,2,figsize=(14,10))
    frameID=frames[0]
    off=offset[0]
    st=ar.loadPhys("Phys/nphys_",frameID)

    cubeX,cubeY,cubeFX,cubeFY,cubeP,cubeFXs,cubeFYs,cubePs,cubeXD,cubeYD,cubeXDs,cubeYDs,xAv,yAv,fxAv,fxAv,fyAv,pAv,fxAvs,fyAvs,pAvs,rmsVal,rmsX,rmsY,rmsVals,rmsXs,rmsYs=ar.loadCubes(loadPref,source,frameID)
    ymi,yma=ar.bounds(cubeYDs[:,off],2)
    sc1=ax[0,0].scatter(xAv,yAv,c=cubeYDs[:,off],vmin=ymi,vmax=yma)
    ax[0,0].set_title(str(int(frameID))+" ["+st+"]")
    fig.colorbar(sc1,ax=ax[0,0])

    frameID=frames[1]
    off=offset[1]
    st=ar.loadPhys("Phys/nphys_",frameID)
    cubeX,cubeY,cubeFX,cubeFY,cubeP,cubeFXs,cubeFYs,cubePs,cubeXD,cubeYD,cubeXDs,cubeYDs,xAv,yAv,fxAv,fxAv,fyAv,pAv,fxAvs,fyAvs,pAvs,rmsVal,rmsX,rmsY,rmsVals,rmsXs,rmsYs=ar.loadCubes(loadPref,source,frameID)
    ymi,yma=ar.bounds(cubeYDs[:,off],2)
    sc1=ax[0,1].scatter(xAv,yAv,c=cubeYDs[:,off],vmin=ymi,vmax=yma)
    ax[0,1].set_title(str(int(frameID))+" ["+st+"]")
    fig.colorbar(sc1,ax=ax[0,1])

    frameID=frames[2]
    off=offset[2]
    st=ar.loadPhys("Phys/nphys_",frameID)
    cubeX,cubeY,cubeFX,cubeFY,cubeP,cubeFXs,cubeFYs,cubePs,cubeXD,cubeYD,cubeXDs,cubeYDs,xAv,yAv,fxAv,fxAv,fyAv,pAv,fxAvs,fyAvs,pAvs,rmsVal,rmsX,rmsY,rmsVals,rmsXs,rmsYs=ar.loadCubes(loadPref,source,frameID)
    ymi,yma=ar.bounds(cubeYDs[:,off],2)
    sc1=ax[1,0].scatter(xAv,yAv,c=cubeYDs[:,off],vmin=ymi,vmax=yma)
    ax[1,0].set_title(str(int(frameID))+" ["+st+"]")
    fig.colorbar(sc1,ax=ax[1,0])

    frameID=frames[3]
    off=offset[3]
    st=ar.loadPhys("Phys/nphys_",frameID)
    cubeX,cubeY,cubeFX,cubeFY,cubeP,cubeFXs,cubeFYs,cubePs,cubeXD,cubeYD,cubeXDs,cubeYDs,xAv,yAv,fxAv,fxAv,fyAv,pAv,fxAvs,fyAvs,pAvs,rmsVal,rmsX,rmsY,rmsVals,rmsXs,rmsYs=ar.loadCubes(loadPref,source,frameID)
    ymi,yma=ar.bounds(cubeYDs[:,off],2)
    sc1=ax[1,1].scatter(xAv,yAv,c=cubeYDs[:,off],vmin=ymi,vmax=yma)
    ax[1,1].set_title(str(int(frameID))+" ["+st+"]")
    fig.colorbar(sc1,ax=ax[1,1])

    plt.suptitle("dy with Rotation (smooth component subtracted)")
    plt.savefig("dyRotS.png")

    #----------------------------------------------------------------------

    frames=np.array([15854,15864,15874,15884]).astype('int')
    offset=[4,1,5,0]

    loadPref="dbSet"
    source="w"

    fig,ax=plt.subplots(2,2,figsize=(14,10))
    frameID=frames[0]
    off=offset[0]
    st=ar.loadPhys("Phys/nphys_",frameID)

    cubeX,cubeY,cubeFX,cubeFY,cubeP,cubeFXs,cubeFYs,cubePs,cubeXD,cubeYD,cubeXDs,cubeYDs,xAv,yAv,fxAv,fxAv,fyAv,pAv,fxAvs,fyAvs,pAvs,rmsVal,rmsX,rmsY,rmsVals,rmsXs,rmsYs=ar.loadCubes(loadPref,source,frameID)
    fxmi,fxma=ar.bounds(cubeFXs.mean(axis=1),2)
    sc1=ax[0,0].scatter(xAv,yAv,c=cubeFXs.mean(axis=1),vmin=fxmi,vmax=fxma)
    ax[0,0].set_title(str(int(frameID))+" ["+st+"]")
    fig.colorbar(sc1,ax=ax[0,0])

    frameID=frames[1]
    off=offset[1]
    st=ar.loadPhys("Phys/nphys_",frameID)
    cubeX,cubeY,cubeFX,cubeFY,cubeP,cubeFXs,cubeFYs,cubePs,cubeXD,cubeYD,cubeXDs,cubeYDs,xAv,yAv,fxAv,fxAv,fyAv,pAv,fxAvs,fyAvs,pAvs,rmsVal,rmsX,rmsY,rmsVals,rmsXs,rmsYs=ar.loadCubes(loadPref,source,frameID)
    fxmi,fxma=ar.bounds(cubeFXs.mean(axis=1),2)
    sc1=ax[0,1].scatter(xAv,yAv,c=cubeFXs.mean(axis=1),vmin=fxmi,vmax=fxma)
    ax[0,1].set_title(str(int(frameID))+" ["+st+"]")
    fig.colorbar(sc1,ax=ax[0,1])

    frameID=frames[2]
    off=offset[2]
    st=ar.loadPhys("Phys/nphys_",frameID)
    cubeX,cubeY,cubeFX,cubeFY,cubeP,cubeFXs,cubeFYs,cubePs,cubeXD,cubeYD,cubeXDs,cubeYDs,xAv,yAv,fxAv,fxAv,fyAv,pAv,fxAvs,fyAvs,pAvs,rmsVal,rmsX,rmsY,rmsVals,rmsXs,rmsYs=ar.loadCubes(loadPref,source,frameID)
    fxmi,fxma=ar.bounds(cubeFXs.mean(axis=1),2)
    sc1=ax[1,0].scatter(xAv,yAv,c=cubeFXs.mean(axis=1),vmin=fxmi,vmax=fxma)
    ax[1,0].set_title(str(int(frameID))+" ["+st+"]")
    fig.colorbar(sc1,ax=ax[1,0])

    frameID=frames[3]
    off=offset[3]
    st=ar.loadPhys("Phys/nphys_",frameID)
    cubeX,cubeY,cubeFX,cubeFY,cubeP,cubeFXs,cubeFYs,cubePs,cubeXD,cubeYD,cubeXDs,cubeYDs,xAv,yAv,fxAv,fxAv,fyAv,pAv,fxAvs,fyAvs,pAvs,rmsVal,rmsX,rmsY,rmsVals,rmsXs,rmsYs=ar.loadCubes(loadPref,source,frameID)
    fxmi,fxma=ar.bounds(cubeFXs.mean(axis=1),2)
    sc1=ax[1,1].scatter(xAv,yAv,c=cubeFXs.mean(axis=1),vmin=fxmi,vmax=fxma)
    ax[1,1].set_title(str(int(frameID))+" ["+st+"]")
    fig.colorbar(sc1,ax=ax[1,1])

    plt.suptitle("Average Spot Size (x) with Rotation (smooth component subtracted)")
    plt.savefig("fxRotS.png")

    #----------------------------------------------------------------------

    frames=np.array([15854,15864,15874,15884]).astype('int')
    offset=[4,1,5,0]

    loadPref="dbSet"
    source="w"

    fig,ax=plt.subplots(2,2,figsize=(14,10))
    frameID=frames[0]
    off=offset[0]
    st=ar.loadPhys("Phys/nphys_",frameID)

    cubeX,cubeY,cubeFX,cubeFY,cubeP,cubeFXs,cubeFYs,cubePs,cubeXD,cubeYD,cubeXDs,cubeYDs,xAv,yAv,fxAv,fxAv,fyAv,pAv,fxAvs,fyAvs,pAvs,rmsVal,rmsX,rmsY,rmsVals,rmsXs,rmsYs=ar.loadCubes(loadPref,source,frameID)
    fymi,fyma=ar.bounds(cubeFYs.mean(axis=1),2)
    sc1=ax[0,0].scatter(xAv,yAv,c=cubeFYs.mean(axis=1),vmin=fymi,vmax=fyma)
    ax[0,0].set_title(str(int(frameID))+" ["+st+"]")
    fig.colorbar(sc1,ax=ax[0,0])

    frameID=frames[1]
    off=offset[1]
    st=ar.loadPhys("Phys/nphys_",frameID)
    cubeX,cubeY,cubeFX,cubeFY,cubeP,cubeFXs,cubeFYs,cubePs,cubeXD,cubeYD,cubeXDs,cubeYDs,xAv,yAv,fxAv,fxAv,fyAv,pAv,fxAvs,fyAvs,pAvs,rmsVal,rmsX,rmsY,rmsVals,rmsXs,rmsYs=ar.loadCubes(loadPref,source,frameID)
    fymi,fyma=ar.bounds(cubeFYs.mean(axis=1),2)
    sc1=ax[0,1].scatter(xAv,yAv,c=cubeFYs.mean(axis=1),vmin=fymi,vmax=fyma)
    ax[0,1].set_title(str(int(frameID))+" ["+st+"]")
    fig.colorbar(sc1,ax=ax[0,1])

    frameID=frames[2]
    off=offset[2]
    st=ar.loadPhys("Phys/nphys_",frameID)
    cubeX,cubeY,cubeFX,cubeFY,cubeP,cubeFXs,cubeFYs,cubePs,cubeXD,cubeYD,cubeXDs,cubeYDs,xAv,yAv,fxAv,fxAv,fyAv,pAv,fxAvs,fyAvs,pAvs,rmsVal,rmsX,rmsY,rmsVals,rmsXs,rmsYs=ar.loadCubes(loadPref,source,frameID)
    fymi,fyma=ar.bounds(cubeFYs.mean(axis=1),2)
    sc1=ax[1,0].scatter(xAv,yAv,c=cubeFYs.mean(axis=1),vmin=fymi,vmax=fyma)
    ax[1,0].set_title(str(int(frameID))+" ["+st+"]")
    fig.colorbar(sc1,ax=ax[1,0])

    frameID=frames[3]
    off=offset[3]
    st=ar.loadPhys("Phys/nphys_",frameID)
    cubeX,cubeY,cubeFX,cubeFY,cubeP,cubeFXs,cubeFYs,cubePs,cubeXD,cubeYD,cubeXDs,cubeYDs,xAv,yAv,fxAv,fxAv,fyAv,pAv,fxAvs,fyAvs,pAvs,rmsVal,rmsX,rmsY,rmsVals,rmsXs,rmsYs=ar.loadCubes(loadPref,source,frameID)
    sc1=ax[1,1].scatter(xAv,yAv,c=cubeFYs.mean(axis=1),vmin=fymi,vmax=fyma)
    fymi,fyma=ar.bounds(cubeFYs.mean(axis=1),2)
    ax[1,1].set_title(str(int(frameID))+" ["+st+"]")
    fig.colorbar(sc1,ax=ax[1,1])

    plt.suptitle("Average Spot Size (y) with Rotation (smooth component subtracted)")
    plt.savefig("fyRotS.png")

    #----------------------------------------------------------------------

def byFrameFX():


    """
    Plot the frame to frame variation in spot size (x) for four consecutive frames  
    for four different exposure times.
    Fiugres 6-8
    """

    loadPref="dbSet"
    source="w"

    frameBase=1762200
    offsets=[5,6,7,8]



    cubeX,cubeY,cubeFX,cubeFY,cubeP,cubeFXs,cubeFYs,cubePs,cubeXD,cubeYD,cubeXDs,cubeYDs,xAv,yAv,fxAv,fxAv,fyAv,pAv,fxAvs,fyAvs,pAvs,rmsVal,rmsX,rmsY,rmsVals,rmsXs,rmsYs=ar.loadCubes(loadPref,source,frameBase)
    st=ar.loadPhys("Phys/nphys_",frameBase)

    
    fig,ax=plt.subplots(2,2,figsize=(14,10))

    frameID=frameBase+offsets[0]
    offset=offsets[0]
    sc1=ax[0,0].scatter(xAv,yAv,c=cubeFX[:,offset],vmin=0.6,vmax=1.4)
    ax[0,0].set_title(str(int(frameID)))
    fig.colorbar(sc1,ax=ax[0,0])

    frameID=frameBase+offsets[1]
    offset=offsets[1]
    sc1=ax[0,1].scatter(xAv,yAv,c=cubeFX[:,offset],vmin=0.6,vmax=1.4)
    ax[0,1].set_title(str(int(frameID)))
    fig.colorbar(sc1,ax=ax[0,1])

    frameID=frameBase+offsets[2]
    offset=offsets[2]
    sc1=ax[1,0].scatter(xAv,yAv,c=cubeFX[:,offset],vmin=0.6,vmax=1.4)
    ax[1,0].set_title(str(int(frameID)))
    fig.colorbar(sc1,ax=ax[1,0])

    frameID=frameBase+offsets[3]
    offset=offsets[3]
    sc1=ax[1,1].scatter(xAv,yAv,c=cubeFX[:,offset],vmin=0.6,vmax=1.4)
    ax[1,1].set_title(str(int(frameID)))
    fig.colorbar(sc1,ax=ax[1,1])

    plt.suptitle("Frame to Frame variation of Spot Size (x) ["+st+"]")
    plt.savefig("fxVar05.png")

    #----------------------------------------

    frameBase=1762300
    offsets=[0,1,2,3]



    cubeX,cubeY,cubeFX,cubeFY,cubeP,cubeFXs,cubeFYs,cubePs,cubeXD,cubeYD,cubeXDs,cubeYDs,xAv,yAv,fxAv,fxAv,fyAv,pAv,fxAvs,fyAvs,pAvs,rmsVal,rmsX,rmsY,rmsVals,rmsXs,rmsYs=ar.loadCubes(loadPref,source,frameBase)
    st=ar.loadPhys("Phys/nphys_",frameBase)

    
    fig,ax=plt.subplots(2,2,figsize=(14,10))

    frameID=frameBase+offsets[0]
    offset=offsets[0]
    sc1=ax[0,0].scatter(xAv,yAv,c=cubeFX[:,offset],vmin=0.6,vmax=1.4)
    ax[0,0].set_title(str(int(frameID)))
    fig.colorbar(sc1,ax=ax[0,0])

    frameID=frameBase+offsets[1]
    offset=offsets[1]
    sc1=ax[0,1].scatter(xAv,yAv,c=cubeFX[:,offset],vmin=0.6,vmax=1.4)
    ax[0,1].set_title(str(int(frameID)))
    fig.colorbar(sc1,ax=ax[0,1])

    frameID=frameBase+offsets[2]
    offset=offsets[2]
    sc1=ax[1,0].scatter(xAv,yAv,c=cubeFX[:,offset],vmin=0.6,vmax=1.4)
    ax[1,0].set_title(str(int(frameID)))
    fig.colorbar(sc1,ax=ax[1,0])

    frameID=frameBase+offsets[3]
    offset=offsets[3]
    sc1=ax[1,1].scatter(xAv,yAv,c=cubeFX[:,offset],vmin=0.6,vmax=1.4)
    ax[1,1].set_title(str(int(frameID)))
    fig.colorbar(sc1,ax=ax[1,1])

    plt.suptitle("Frame to Frame variation of Spot Size (x) ["+st+"]")
    plt.savefig("fxVar10.png")

    #----------------------------------------

    frameBase=1762500
    offsets=[0,1,2,3]



    cubeX,cubeY,cubeFX,cubeFY,cubeP,cubeFXs,cubeFYs,cubePs,cubeXD,cubeYD,cubeXDs,cubeYDs,xAv,yAv,fxAv,fxAv,fyAv,pAv,fxAvs,fyAvs,pAvs,rmsVal,rmsX,rmsY,rmsVals,rmsXs,rmsYs=ar.loadCubes(loadPref,source,frameBase)
    st=ar.loadPhys("Phys/nphys_",frameBase)

    
    fig,ax=plt.subplots(2,2,figsize=(14,10))

    frameID=frameBase+offsets[0]
    offset=offsets[0]
    sc1=ax[0,0].scatter(xAv,yAv,c=cubeFX[:,offset],vmin=0.6,vmax=1.4)
    ax[0,0].set_title(str(int(frameID)))
    fig.colorbar(sc1,ax=ax[0,0])

    frameID=frameBase+offsets[1]
    offset=offsets[1]
    sc1=ax[0,1].scatter(xAv,yAv,c=cubeFX[:,offset],vmin=0.6,vmax=1.4)
    ax[0,1].set_title(str(int(frameID)))
    fig.colorbar(sc1,ax=ax[0,1])

    frameID=frameBase+offsets[2]
    offset=offsets[2]
    sc1=ax[1,0].scatter(xAv,yAv,c=cubeFX[:,offset],vmin=0.6,vmax=1.4)
    ax[1,0].set_title(str(int(frameID)))
    fig.colorbar(sc1,ax=ax[1,0])

    frameID=frameBase+offsets[3]
    offset=offsets[3]
    sc1=ax[1,1].scatter(xAv,yAv,c=cubeFX[:,offset],vmin=0.6,vmax=1.4)
    ax[1,1].set_title(str(int(frameID)))
    fig.colorbar(sc1,ax=ax[1,1])

    plt.suptitle("Frame to Frame variation of Spot Size (x) ["+st+"]")
    plt.savefig("fxVar50.png")

    #----------------------------------------

def byFrameFY():

    
    """
    Plot the frame to frame variation in spot size (y) for four consecutive frames 
    for four different exposure times.
    Figures 9-11
    """

    loadPref="dbSet"
    source="w"

    frameBase=1762200
    offsets=[5,6,7,8]



    cubeX,cubeY,cubeFX,cubeFY,cubeP,cubeFXs,cubeFYs,cubePs,cubeXD,cubeYD,cubeXDs,cubeYDs,xAv,yAv,fxAv,fxAv,fyAv,pAv,fxAvs,fyAvs,pAvs,rmsVal,rmsX,rmsY,rmsVals,rmsXs,rmsYs=ar.loadCubes(loadPref,source,frameBase)
    st=ar.loadPhys("Phys/nphys_",frameBase)

    
    fig,ax=plt.subplots(2,2,figsize=(14,10))

    frameID=frameBase+offsets[0]
    offset=offsets[0]
    sc1=ax[0,0].scatter(xAv,yAv,c=cubeFY[:,offset],vmin=0.6,vmax=1.4)
    ax[0,0].set_title(str(int(frameID)))
    fig.colorbar(sc1,ax=ax[0,0])

    frameID=frameBase+offsets[1]
    offset=offsets[1]
    sc1=ax[0,1].scatter(xAv,yAv,c=cubeFY[:,offset],vmin=0.6,vmax=1.4)
    ax[0,1].set_title(str(int(frameID)))
    fig.colorbar(sc1,ax=ax[0,1])

    frameID=frameBase+offsets[2]
    offset=offsets[2]
    sc1=ax[1,0].scatter(xAv,yAv,c=cubeFY[:,offset],vmin=0.6,vmax=1.4)
    ax[1,0].set_title(str(int(frameID)))
    fig.colorbar(sc1,ax=ax[1,0])

    frameID=frameBase+offsets[3]
    offset=offsets[3]
    sc1=ax[1,1].scatter(xAv,yAv,c=cubeFY[:,offset],vmin=0.6,vmax=1.4)
    ax[1,1].set_title(str(int(frameID)))
    fig.colorbar(sc1,ax=ax[1,1])

    plt.suptitle("Frame to Frame variation of Spot Size (y) ["+st+"]")
    plt.savefig("fyVar05.png")

    #----------------------------------------

    frameBase=1762300
    offsets=[0,1,2,3]



    cubeX,cubeY,cubeFX,cubeFY,cubeP,cubeFXs,cubeFYs,cubePs,cubeXD,cubeYD,cubeXDs,cubeYDs,xAv,yAv,fxAv,fxAv,fyAv,pAv,fxAvs,fyAvs,pAvs,rmsVal,rmsX,rmsY,rmsVals,rmsXs,rmsYs=ar.loadCubes(loadPref,source,frameBase)
    st=ar.loadPhys("Phys/nphys_",frameBase)

    
    fig,ax=plt.subplots(2,2,figsize=(14,10))

    frameID=frameBase+offsets[0]
    offset=offsets[0]
    sc1=ax[0,0].scatter(xAv,yAv,c=cubeFY[:,offset],vmin=0.6,vmax=1.4)
    ax[0,0].set_title(str(int(frameID)))
    fig.colorbar(sc1,ax=ax[0,0])

    frameID=frameBase+offsets[1]
    offset=offsets[1]
    sc1=ax[0,1].scatter(xAv,yAv,c=cubeFY[:,offset],vmin=0.6,vmax=1.4)
    ax[0,1].set_title(str(int(frameID)))
    fig.colorbar(sc1,ax=ax[0,1])

    frameID=frameBase+offsets[2]
    offset=offsets[2]
    sc1=ax[1,0].scatter(xAv,yAv,c=cubeFY[:,offset],vmin=0.6,vmax=1.4)
    ax[1,0].set_title(str(int(frameID)))
    fig.colorbar(sc1,ax=ax[1,0])

    frameID=frameBase+offsets[3]
    offset=offsets[3]
    sc1=ax[1,1].scatter(xAv,yAv,c=cubeFY[:,offset],vmin=0.6,vmax=1.4)
    ax[1,1].set_title(str(int(frameID)))
    fig.colorbar(sc1,ax=ax[1,1])

    plt.suptitle("Frame to Frame variation of Spot Size (y) ["+st+"]")
    plt.savefig("fyVar10.png")

    #----------------------------------------

    frameBase=1762500
    offsets=[0,1,2,3]



    cubeX,cubeY,cubeFX,cubeFY,cubeP,cubeFXs,cubeFYs,cubePs,cubeXD,cubeYD,cubeXDs,cubeYDs,xAv,yAv,fxAv,fxAv,fyAv,pAv,fxAvs,fyAvs,pAvs,rmsVal,rmsX,rmsY,rmsVals,rmsXs,rmsYs=ar.loadCubes(loadPref,source,frameBase)
    st=ar.loadPhys("Phys/nphys_",frameBase)

    
    fig,ax=plt.subplots(2,2,figsize=(14,10))

    frameID=frameBase+offsets[0]
    offset=offsets[0]
    sc1=ax[0,0].scatter(xAv,yAv,c=cubeFY[:,offset],vmin=0.6,vmax=1.4)
    ax[0,0].set_title(str(int(frameID)))
    fig.colorbar(sc1,ax=ax[0,0])

    frameID=frameBase+offsets[1]
    offset=offsets[1]
    sc1=ax[0,1].scatter(xAv,yAv,c=cubeFY[:,offset],vmin=0.6,vmax=1.4)
    ax[0,1].set_title(str(int(frameID)))
    fig.colorbar(sc1,ax=ax[0,1])

    frameID=frameBase+offsets[2]
    offset=offsets[2]
    sc1=ax[1,0].scatter(xAv,yAv,c=cubeFY[:,offset],vmin=0.6,vmax=1.4)
    ax[1,0].set_title(str(int(frameID)))
    fig.colorbar(sc1,ax=ax[1,0])

    frameID=frameBase+offsets[3]
    offset=offsets[3]
    sc1=ax[1,1].scatter(xAv,yAv,c=cubeFY[:,offset],vmin=0.6,vmax=1.4)
    ax[1,1].set_title(str(int(frameID)))
    fig.colorbar(sc1,ax=ax[1,1])

    plt.suptitle("Frame to Frame variation of Spot Size (y) ["+st+"]")
    plt.savefig("fyVar50.png")

    #----------------------------------------

def byFrameP():

    """
    Plot the frame to frame variation in peak brightness for four consecutive frames  
    for four different exposure times.
    Figures 12-14
    """

    
    loadPref="dbSet"
    source="l"

    frameBase=1762200
    offsets=[5,6,7,8]

    
    cubeX,cubeY,cubeFX,cubeFY,cubeP,cubeFXs,cubeFYs,cubePs,cubeXD,cubeYD,cubeXDs,cubeYDs,xAv,yAv,fxAv,fxAv,fyAv,pAv,fxAvs,fyAvs,pAvs,rmsVal,rmsX,rmsY,rmsVals,rmsXs,rmsYs=ar.loadCubes(loadPref,source,frameBase)
    st=ar.loadPhys("Phys/nphys_",frameBase)

    
    fig,ax=plt.subplots(2,2,figsize=(14,10))

    frameID=frameBase+offsets[0]
    offset=offsets[0]
    sc1=ax[0,0].scatter(xAv,yAv,c=cubeP[:,offset],vmin=1500,vmax=5000)
    ax[0,0].set_title(str(int(frameID)))
    fig.colorbar(sc1,ax=ax[0,0])

    frameID=frameBase+offsets[1]
    offset=offsets[1]
    sc1=ax[0,1].scatter(xAv,yAv,c=cubeP[:,offset],vmin=1500,vmax=5000)
    ax[0,1].set_title(str(int(frameID)))
    fig.colorbar(sc1,ax=ax[0,1])

    frameID=frameBase+offsets[2]
    offset=offsets[2]
    sc1=ax[1,0].scatter(xAv,yAv,c=cubeP[:,offset],vmin=1500,vmax=5000)
    ax[1,0].set_title(str(int(frameID)))
    fig.colorbar(sc1,ax=ax[1,0])

    frameID=frameBase+offsets[3]
    offset=offsets[3]
    sc1=ax[1,1].scatter(xAv,yAv,c=cubeP[:,offset],vmin=1500,vmax=5000)
    ax[1,1].set_title(str(int(frameID)))
    fig.colorbar(sc1,ax=ax[1,1])

    plt.suptitle("Frame to Frame variation of Brightness ["+st+"]")
    plt.savefig("pVar05.png")

    #----------------------------------------

    frameBase=1762300
    offsets=[0,1,2,3]



    cubeX,cubeY,cubeFX,cubeFY,cubeP,cubeFXs,cubeFYs,cubePs,cubeXD,cubeYD,cubeXDs,cubeYDs,xAv,yAv,fxAv,fxAv,fyAv,pAv,fxAvs,fyAvs,pAvs,rmsVal,rmsX,rmsY,rmsVals,rmsXs,rmsYs=ar.loadCubes(loadPref,source,frameBase)
    st=ar.loadPhys("Phys/nphys_",frameBase)

    
    fig,ax=plt.subplots(2,2,figsize=(14,10))

    frameID=frameBase+offsets[0]
    offset=offsets[0]
    sc1=ax[0,0].scatter(xAv,yAv,c=cubeP[:,offset],vmin=3500,vmax=7000)
    ax[0,0].set_title(str(int(frameID)))
    fig.colorbar(sc1,ax=ax[0,0])

    frameID=frameBase+offsets[1]
    offset=offsets[1]
    sc1=ax[0,1].scatter(xAv,yAv,c=cubeP[:,offset],vmin=3500,vmax=7000)
    ax[0,1].set_title(str(int(frameID)))
    fig.colorbar(sc1,ax=ax[0,1])

    frameID=frameBase+offsets[2]
    offset=offsets[2]
    sc1=ax[1,0].scatter(xAv,yAv,c=cubeP[:,offset],vmin=3500,vmax=7000)
    ax[1,0].set_title(str(int(frameID)))
    fig.colorbar(sc1,ax=ax[1,0])

    frameID=frameBase+offsets[3]
    offset=offsets[3]
    sc1=ax[1,1].scatter(xAv,yAv,c=cubeP[:,offset],vmin=3500,vmax=7000)
    ax[1,1].set_title(str(int(frameID)))
    fig.colorbar(sc1,ax=ax[1,1])

    plt.suptitle("Frame to Frame variation of Brightness["+st+"]")
    plt.savefig("pVar10.png")

    #----------------------------------------

    frameBase=1762500
    offsets=[0,1,2,3]



    cubeX,cubeY,cubeFX,cubeFY,cubeP,cubeFXs,cubeFYs,cubePs,cubeXD,cubeYD,cubeXDs,cubeYDs,xAv,yAv,fxAv,fxAv,fyAv,pAv,fxAvs,fyAvs,pAvs,rmsVal,rmsX,rmsY,rmsVals,rmsXs,rmsYs=ar.loadCubes(loadPref,source,frameBase)
    st=ar.loadPhys("Phys/nphys_",frameBase)

    
    fig,ax=plt.subplots(2,2,figsize=(14,10))

    frameID=frameBase+offsets[0]
    offset=offsets[0]
    sc1=ax[0,0].scatter(xAv,yAv,c=cubeP[:,offset],vmin=6000,vmax=20000)
    ax[0,0].set_title(str(int(frameID)))
    fig.colorbar(sc1,ax=ax[0,0])

    frameID=frameBase+offsets[1]
    offset=offsets[1]
    sc1=ax[0,1].scatter(xAv,yAv,c=cubeP[:,offset],vmin=6000,vmax=20000)
    ax[0,1].set_title(str(int(frameID)))
    fig.colorbar(sc1,ax=ax[0,1])

    frameID=frameBase+offsets[2]
    offset=offsets[2]
    sc1=ax[1,0].scatter(xAv,yAv,c=cubeP[:,offset],vmin=6000,vmax=20000)
    ax[1,0].set_title(str(int(frameID)))
    fig.colorbar(sc1,ax=ax[1,0])

    frameID=frameBase+offsets[3]
    offset=offsets[3]
    sc1=ax[1,1].scatter(xAv,yAv,c=cubeP[:,offset],vmin=6000,vmax=20000)
    ax[1,1].set_title(str(int(frameID)))
    fig.colorbar(sc1,ax=ax[1,1])

    plt.suptitle("Frame to Frame variation of Brightness ["+st+"]")
    plt.savefig("pVar50.png")

    #----------------------------------------


def fits():

    """

    demonstrate the effects of subtracting the smooth component. 
    Figure 16

    """
    
    loadPref="dbSet/"
    source="l"
    
    frames=np.array([15854,15864,15874,15884]).astype('int')
    offset=[4,1,5,0]


    frameID=frames[0]
    
    
    cubeX,cubeY,cubeFX,cubeFY,cubeP,cubeFXs,cubeFYs,cubePs,cubeXD,cubeYD,cubeXDs,cubeYDs,xAv,yAv,fxAv,fxAv,fyAv,pAv,fxAvs,fyAvs,pAvs,rmsVal,rmsX,rmsY,rmsVals,rmsXs,rmsYs=ar.loadCubes(loadPref,source,frameID)

    st=ar.loadPhys("Phys/nphys_",frames[0])

    fig,ax=plt.subplots(2,2,figsize=(14,10))
    
    xmi,xma=ar.bounds(cubeXD[:,4],2)
    ymi,yma=ar.bounds(cubeFX[:,4],2)
    fxmi,fxma=ar.bounds(cubeXD[:,4]-cubeXDs[:,4],2)
    fymi,fyma=ar.bounds(cubeFX[:,4]-cubeFXs[:,4],2)
    
    sc1=ax[0,0].scatter(xAv,yAv,c=cubeXD[:,4],vmin=xmi,vmax=xma)
    ax[0,0].set_title("dx: raw values")
    fig.colorbar(sc1,ax=ax[0,0])

    sc1=ax[0,1].scatter(xAv,yAv,c=cubeFX[:,4],vmin=ymi,vmax=yma)
    ax[0,1].set_title("Spot Size (x): raw values")
    fig.colorbar(sc1,ax=ax[0,1])

    sc1=ax[1,0].scatter(xAv,yAv,c=cubeXD[:,4]-cubeXDs[:,4],vmin=fxmi,vmax=fxma)
    ax[1,0].set_title("dx: smooth component only")
    fig.colorbar(sc1,ax=ax[1,0])

    sc1=ax[1,1].scatter(xAv,yAv,c=cubeFX[:,4]-cubeFXs[:,4],vmin=fymi,vmax=fyma)
    ax[1,1].set_title("Spot Size (x): smooth component only")
    fig.colorbar(sc1,ax=ax[1,1])

    plt.suptitle("Effects of Subtraction - "+str(int(frameID+4))+" ["+st+"]")
    plt.savefig("fitX.png")


def elevation():

    """
    demonstrate the variation of hte pattern with elevation
    Figure 16
    """
    
    loadPref="dbSet"
    source="w"
    
    frames=np.array([1762500,1762900]).astype('int')
    offset=[4,1,5,0]
     
    
    fig,ax=plt.subplots(2,2,figsize=(14,10))

    frameID=frames[0]
    st=ar.loadPhys("Phys/nphys_",frames[0])
    
    cubeX,cubeY,cubeFX,cubeFY,cubeP,cubeFXs,cubeFYs,cubePs,cubeXD,cubeYD,cubeXDs,cubeYDs,xAv,yAv,fxAv,fxAv,fyAv,pAv,fxAvs,fyAvs,pAvs,rmsVal,rmsX,rmsY,rmsVals,rmsXs,rmsYs=ar.loadCubes(loadPref,source,frameID)

    xmi,xma=ar.bounds(cubeXDs[:,4],2)
    fxmi,fxma=ar.bounds(cubeFXs[:,4],1.5)
    
    sc1=ax[0,0].scatter(xAv,yAv,c=cubeXDs[:,1])
    ax[0,0].set_title("dx (smooth component subtracted) [El=30]")
    fig.colorbar(sc1,ax=ax[0,0])

    sc1=ax[0,1].scatter(xAv,yAv,c=cubeFXs[:,1],vmin=fxmi,vmax=fxma)
    ax[0,1].set_title("spot size (smooth component subtracted) [El=30]")
    fig.colorbar(sc1,ax=ax[0,1])

    
    frameID=frames[1]
    cubeX,cubeY,cubeFX,cubeFY,cubeP,cubeFXs,cubeFYs,cubePs,cubeXD,cubeYD,cubeXDs,cubeYDs,xAv,yAv,fxAv,fxAv,fyAv,pAv,fxAvs,fyAvs,pAvs,rmsVal,rmsX,rmsY,rmsVals,rmsXs,rmsYs=ar.loadCubes(loadPref,source,frameID)

    #xmi,xma=ar.bounds(cubeXDs[:,4],2)
    #fxmi,fxma=ar.bounds(cubeFXs[:,4],1.5)
    
    sc1=ax[1,0].scatter(xAv,yAv,c=cubeXDs[:,1])
    ax[1,0].set_title("dx (smooth component subtracted) [El=45]")
    fig.colorbar(sc1,ax=ax[1,0])

    sc1=ax[1,1].scatter(xAv,yAv,c=cubeFXs[:,1],vmin=fxmi,vmax=fxma)
    ax[1,1].set_title("spot size (smooth component subtracted) [El=45]")
    fig.colorbar(sc1,ax=ax[1,1])

    plt.suptitle("Elevation Effects: 30 and 45 degrees")

    plt.savefig("elevation.png")

def byFrameDX():

    """
    plot the delta x (difference from mean position) for four consecutive frames 
    for four exposure times.
    Figures 17-19
    """

    
    source="l"
    loadPref="dbSet"

    frameBase=1762200
    offsets=[5,6,7,8]
            

    cubeX,cubeY,cubeFX,cubeFY,cubeP,cubeFXs,cubeFYs,cubePs,cubeXD,cubeYD,cubeXDs,cubeYDs,xAv,yAv,fxAv,fxAv,fyAv,pAv,fxAvs,fyAvs,pAvs,rmsVal,rmsX,rmsY,rmsVals,rmsXs,rmsYs=ar.loadCubes(loadPref,source,frameBase)
    st=ar.loadPhys("Phys/nphys_",frameBase)

    
    fig,ax=plt.subplots(2,2,figsize=(14,10))

    frameID=frameBase+offsets[0]
    offset=offsets[0]
    sc1=ax[0,0].scatter(xAv,yAv,c=cubeXD[:,offset])
    ax[0,0].set_title(str(int(frameID)))
    fig.colorbar(sc1,ax=ax[0,0])

    frameID=frameBase+offsets[1]
    offset=offsets[1]
    sc1=ax[0,1].scatter(xAv,yAv,c=cubeXD[:,offset])
    ax[0,1].set_title(str(int(frameID)))
    fig.colorbar(sc1,ax=ax[0,1])

    frameID=frameBase+offsets[2]
    offset=offsets[2]
    sc1=ax[1,0].scatter(xAv,yAv,c=cubeXD[:,offset])
    ax[1,0].set_title(str(int(frameID)))
    fig.colorbar(sc1,ax=ax[1,0])

    frameID=frameBase+offsets[3]
    offset=offsets[3]
    sc1=ax[1,1].scatter(xAv,yAv,c=cubeXD[:,offset])
    ax[1,1].set_title(str(int(frameID)))
    fig.colorbar(sc1,ax=ax[1,1])

    plt.suptitle("Frame to Frame variation of dx ["+st+"]")
    plt.savefig("dxVar05.png")


    #----------------------------------------
    frameBase=1762300
    offsets=[0,1,2,3]
            

    cubeX,cubeY,cubeFX,cubeFY,cubeP,cubeFXs,cubeFYs,cubePs,cubeXD,cubeYD,cubeXDs,cubeYDs,xAv,yAv,fxAv,fxAv,fyAv,pAv,fxAvs,fyAvs,pAvs,rmsVal,rmsX,rmsY,rmsVals,rmsXs,rmsYs=ar.loadCubes(loadPref,source,frameBase)
    st=ar.loadPhys("Phys/nphys_",frameBase)

    
    fig,ax=plt.subplots(2,2,figsize=(14,10))

    frameID=frameBase+offsets[0]
    offset=offsets[0]
    sc1=ax[0,0].scatter(xAv,yAv,c=cubeXD[:,offset])
    ax[0,0].set_title(str(int(frameID)))
    fig.colorbar(sc1,ax=ax[0,0])

    frameID=frameBase+offsets[1]
    offset=offsets[1]
    sc1=ax[0,1].scatter(xAv,yAv,c=cubeXD[:,offset])
    ax[0,1].set_title(str(int(frameID)))
    fig.colorbar(sc1,ax=ax[0,1])

    frameID=frameBase+offsets[2]
    offset=offsets[2]
    sc1=ax[1,0].scatter(xAv,yAv,c=cubeXD[:,offset])
    ax[1,0].set_title(str(int(frameID)))
    fig.colorbar(sc1,ax=ax[1,0])

    frameID=frameBase+offsets[3]
    offset=offsets[3]
    sc1=ax[1,1].scatter(xAv,yAv,c=cubeXD[:,offset])
    ax[1,1].set_title(str(int(frameID)))
    fig.colorbar(sc1,ax=ax[1,1])

    plt.suptitle("Frame to Frame variation of dx ["+st+"]")
    plt.savefig("dxVar10.png")


    #----------------------------------------
    frameBase=1762500
    offsets=[0,1,2,3]
            

    cubeX,cubeY,cubeFX,cubeFY,cubeP,cubeFXs,cubeFYs,cubePs,cubeXD,cubeYD,cubeXDs,cubeYDs,xAv,yAv,fxAv,fxAv,fyAv,pAv,fxAvs,fyAvs,pAvs,rmsVal,rmsX,rmsY,rmsVals,rmsXs,rmsYs=ar.loadCubes(loadPref,source,frameBase)
    st=ar.loadPhys("Phys/nphys_",frameBase)

    
    fig,ax=plt.subplots(2,2,figsize=(14,10))

    frameID=frameBase+offsets[0]
    offset=offsets[0]
    sc1=ax[0,0].scatter(xAv,yAv,c=cubeXD[:,offset])
    ax[0,0].set_title(str(int(frameID)))
    fig.colorbar(sc1,ax=ax[0,0])

    frameID=frameBase+offsets[1]
    offset=offsets[1]
    sc1=ax[0,1].scatter(xAv,yAv,c=cubeXD[:,offset])
    ax[0,1].set_title(str(int(frameID)))
    fig.colorbar(sc1,ax=ax[0,1])

    frameID=frameBase+offsets[2]
    offset=offsets[2]
    sc1=ax[1,0].scatter(xAv,yAv,c=cubeXD[:,offset])
    ax[1,0].set_title(str(int(frameID)))
    fig.colorbar(sc1,ax=ax[1,0])

    frameID=frameBase+offsets[3]
    offset=offsets[3]
    sc1=ax[1,1].scatter(xAv,yAv,c=cubeXD[:,offset])
    ax[1,1].set_title(str(int(frameID)))
    fig.colorbar(sc1,ax=ax[1,1])

    plt.suptitle("Frame to Frame variation of dx ["+st+"]")
    plt.savefig("dxVar50.png")


    #----------------------------------------

def byFrameDY():

    """
    plot the delta y (difference from mean position) for four consecutive frames 
    for four exposure times.
    Fiugres 20-22
    """

    source="l"
    loadPref="dbSet"

    frameBase=1762200
    offsets=[5,6,7,8]
            

    cubeX,cubeY,cubeFX,cubeFY,cubeP,cubeFXs,cubeFYs,cubePs,cubeXD,cubeYD,cubeXDs,cubeYDs,xAv,yAv,fxAv,fxAv,fyAv,pAv,fxAvs,fyAvs,pAvs,rmsVal,rmsX,rmsY,rmsVals,rmsXs,rmsYs=ar.loadCubes(loadPref,source,frameBase)
    st=ar.loadPhys("Phys/nphys_",frameBase)

    
    fig,ax=plt.subplots(2,2,figsize=(14,10))

    frameID=frameBase+offsets[0]
    offset=offsets[0]
    sc1=ax[0,0].scatter(xAv,yAv,c=cubeYD[:,offset])
    ax[0,0].set_title(str(int(frameID)))
    fig.colorbar(sc1,ax=ax[0,0])

    frameID=frameBase+offsets[1]
    offset=offsets[1]
    sc1=ax[0,1].scatter(xAv,yAv,c=cubeYD[:,offset])
    ax[0,1].set_title(str(int(frameID)))
    fig.colorbar(sc1,ax=ax[0,1])

    frameID=frameBase+offsets[2]
    offset=offsets[2]
    sc1=ax[1,0].scatter(xAv,yAv,c=cubeYD[:,offset])
    ax[1,0].set_title(str(int(frameID)))
    fig.colorbar(sc1,ax=ax[1,0])

    frameID=frameBase+offsets[3]
    offset=offsets[3]
    sc1=ax[1,1].scatter(xAv,yAv,c=cubeYD[:,offset])
    ax[1,1].set_title(str(int(frameID)))
    fig.colorbar(sc1,ax=ax[1,1])

    plt.suptitle("Frame to Frame variation of dy ["+st+"]")
    plt.savefig("dyVar05.png")


    #----------------------------------------
    frameBase=1762300
    offsets=[0,1,2,3]
            

    cubeX,cubeY,cubeFX,cubeFY,cubeP,cubeFXs,cubeFYs,cubePs,cubeXD,cubeYD,cubeXDs,cubeYDs,xAv,yAv,fxAv,fxAv,fyAv,pAv,fxAvs,fyAvs,pAvs,rmsVal,rmsX,rmsY,rmsVals,rmsXs,rmsYs=ar.loadCubes(loadPref,source,frameBase)
    st=ar.loadPhys("Phys/nphys_",frameBase)

    
    fig,ax=plt.subplots(2,2,figsize=(14,10))

    frameID=frameBase+offsets[0]
    offset=offsets[0]
    sc1=ax[0,0].scatter(xAv,yAv,c=cubeYD[:,offset])
    ax[0,0].set_title(str(int(frameID)))
    fig.colorbar(sc1,ax=ax[0,0])

    frameID=frameBase+offsets[1]
    offset=offsets[1]
    sc1=ax[0,1].scatter(xAv,yAv,c=cubeYD[:,offset])
    ax[0,1].set_title(str(int(frameID)))
    fig.colorbar(sc1,ax=ax[0,1])

    frameID=frameBase+offsets[2]
    offset=offsets[2]
    sc1=ax[1,0].scatter(xAv,yAv,c=cubeYD[:,offset])
    ax[1,0].set_title(str(int(frameID)))
    fig.colorbar(sc1,ax=ax[1,0])

    frameID=frameBase+offsets[3]
    offset=offsets[3]
    sc1=ax[1,1].scatter(xAv,yAv,c=cubeYD[:,offset])
    ax[1,1].set_title(str(int(frameID)))
    fig.colorbar(sc1,ax=ax[1,1])

    plt.suptitle("Frame to Frame variation of dy ["+st+"]")
    plt.savefig("dyVar10.png")


    #----------------------------------------
    frameBase=1762500
    offsets=[0,1,2,3]
            

    cubeX,cubeY,cubeFX,cubeFY,cubeP,cubeFXs,cubeFYs,cubePs,cubeXD,cubeYD,cubeXDs,cubeYDs,xAv,yAv,fxAv,fxAv,fyAv,pAv,fxAvs,fyAvs,pAvs,rmsVal,rmsX,rmsY,rmsVals,rmsXs,rmsYs=ar.loadCubes(loadPref,source,frameBase)
    st=ar.loadPhys("Phys/nphys_",frameBase)

    
    fig,ax=plt.subplots(2,2,figsize=(14,10))

    frameID=frameBase+offsets[0]
    offset=offsets[0]
    sc1=ax[0,0].scatter(xAv,yAv,c=cubeYD[:,offset])
    ax[0,0].set_title(str(int(frameID)))
    fig.colorbar(sc1,ax=ax[0,0])

    frameID=frameBase+offsets[1]
    offset=offsets[1]
    sc1=ax[0,1].scatter(xAv,yAv,c=cubeYD[:,offset])
    ax[0,1].set_title(str(int(frameID)))
    fig.colorbar(sc1,ax=ax[0,1])

    frameID=frameBase+offsets[2]
    offset=offsets[2]
    sc1=ax[1,0].scatter(xAv,yAv,c=cubeYD[:,offset])
    ax[1,0].set_title(str(int(frameID)))
    fig.colorbar(sc1,ax=ax[1,0])

    frameID=frameBase+offsets[3]
    offset=offsets[3]
    sc1=ax[1,1].scatter(xAv,yAv,c=cubeYD[:,offset])
    ax[1,1].set_title(str(int(frameID)))
    fig.colorbar(sc1,ax=ax[1,1])

    plt.suptitle("Frame to Frame variation of dy ["+st+"]")
    plt.savefig("dyVar50.png")


    #----------------------------------------


def diagPlot():

    """
    plot an exmaple of a set of diagnostic plots (spot sizes, dx, dy) for a frame. 
    Figure 28
    """
    
    loadPref="dbSet"
    source="w"
    
    frameID=1762500
    offset=1
    #frames=np.array([1762500,1762900]).astype('int')
    
    fig,ax=plt.subplots(2,2,figsize=(14,10))
    

    cubeX,cubeY,cubeFX,cubeFY,cubeP,cubeFXs,cubeFYs,cubePs,cubeXD,cubeYD,cubeXDs,cubeYDs,xAv,yAv,fxAv,fxAv,fyAv,pAv,fxAvs,fyAvs,pAvs,rmsVal,rmsX,rmsY,rmsVals,rmsXs,rmsYs=ar.loadCubes(loadPref,source,frameID)
    st=ar.loadPhys("Phys/nphys_",frameID)
                   
    xmi,xma=ar.bounds(cubeFXs[:,0],1.5)
    ymi,yma=ar.bounds(cubeFYs[:,0],1.5)
    dxmi,dxma=ar.bounds(cubeXDs[:,0],1.5)
    dymi,dyma=ar.bounds(cubeYDs[:,0],1.5)

    sc1=ax[0,0].scatter(xAv,yAv,c=cubeFXs[:,0],vmin=xmi,vmax=xma)
    ax[0,0].set_title("Spot Size (X) (smooth component subtracted)")
    fig.colorbar(sc1,ax=ax[0,0])

    sc1=ax[0,1].scatter(xAv,yAv,c=cubeFYs[:,0],vmin=ymi,vmax=yma)
    ax[0,1].set_title("Spot Size (Y) (smooth component subtracted)")
    fig.colorbar(sc1,ax=ax[0,1])
    source="l"
    cubeX,cubeY,cubeFX,cubeFY,cubeP,cubeFXs,cubeFYs,cubePs,cubeXD,cubeYD,cubeXDs,cubeYDs,xAv,yAv,fxAv,fxAv,fyAv,pAv,fxAvs,fyAvs,pAvs,rmsVal,rmsX,rmsY,rmsVals,rmsXs,rmsYs=ar.loadCubes(loadPref,source,frameID)
    xdmi,xdma=ar.bounds(cubeXDs[:,offset],1.5)
    ydmi,ydma=ar.bounds(cubeYDs[:,offset],1.5)

    
    sc1=ax[1,0].scatter(xAv,yAv,c=cubeXDs[:,0],vmin=xdmi,vmax=xdma)
    ax[1,0].set_title("dx (smooth component subtracted)")
    fig.colorbar(sc1,ax=ax[1,0])

    sc1=ax[1,1].scatter(xAv,yAv,c=cubeYDs[:,0],vmin=ydmi,vmax=ydma)
    ax[1,1].set_title("dy (smooth component subtracted)")
    fig.colorbar(sc1,ax=ax[1,1])
    
    plt.suptitle("Diagnostic Plots for Frame "+str(frameID)+" ["+st+"]")
    plt.savefig("ydir.png")


def yDir():
    """
    plot an exmaple of a set of diagnostic plots (spot sizes, dx, dy) for a frame. 
    Figure 23
    """

    loadPref="dbSet"

    source="w"
    frameID=15854
    offset=4
    #frames=np.array([1762500,1762900]).astype('int')
    
    fig,ax=plt.subplots(2,2,figsize=(14,10))

    cubeX,cubeY,cubeFX,cubeFY,cubeP,cubeFXs,cubeFYs,cubePs,cubeXD,cubeYD,cubeXDs,cubeYDs,xAv,yAv,fxAv,fxAv,fyAv,pAv,fxAvs,fyAvs,pAvs,rmsVal,rmsX,rmsY,rmsVals,rmsXs,rmsYs=ar.loadCubes(loadPref,source,frameID)
    st=ar.loadPhys("Phys/nphys_",frameID)
    
    xmi,xma=ar.bounds(cubeFXs[:,offset],1.5)
    ymi,yma=ar.bounds(cubeFYs[:,offset],1.5)

    sc1=ax[0,0].scatter(xAv,yAv,c=cubeFXs[:,offset],vmin=xmi,vmax=xma)
    ax[0,0].set_title("Spot Size (X) (smooth component subtracted)")
    fig.colorbar(sc1,ax=ax[0,0])

    sc1=ax[0,1].scatter(xAv,yAv,c=cubeFYs[:,offset],vmin=ymi,vmax=yma)
    ax[0,1].set_title("Spot Size (Y) (smooth component subtracted)")
    fig.colorbar(sc1,ax=ax[0,1])

    source="l"
    cubeX,cubeY,cubeFX,cubeFY,cubeP,cubeFXs,cubeFYs,cubePs,cubeXD,cubeYD,cubeXDs,cubeYDs,xAv,yAv,fxAv,fxAv,fyAv,pAv,fxAvs,fyAvs,pAvs,rmsVal,rmsX,rmsY,rmsVals,rmsXs,rmsYs=ar.loadCubes(loadPref,source,frameID)
    xdmi,xdma=ar.bounds(cubeXDs[:,offset],1.5)
    ydmi,ydma=ar.bounds(cubeYDs[:,offset],1.5)

    
    sc1=ax[1,0].scatter(xAv,yAv,c=cubeXDs[:,offset],vmin=xdmi,vmax=xdma)
    ax[1,0].set_title("dx (smooth component subtracted)")
    fig.colorbar(sc1,ax=ax[1,0])

    sc1=ax[1,1].scatter(xAv,yAv,c=cubeYDs[:,offset],vmin=ydmi,vmax=ydma)
    ax[1,1].set_title("dy (smooth component subtracted)")
    fig.colorbar(sc1,ax=ax[1,1])
    
    plt.suptitle("Diagnostic Plots for Frame "+str(frameID+offset)+" ["+st+"]")
    plt.savefig("diagPlot.png")

  

def rmsRot():

    """
    show the variation in spot motion rms as a function of rotation
    Figure 29
    """
    
    frames=np.array([15854,15864,15874,15884]).astype('int')

    loadPref="dbSet"
    source="l"

    fig,ax=plt.subplots(2,2,figsize=(14,10))
    frameID=frames[0]
    st=ar.loadPhys("Phys/nphys_",frameID)
    cubeX,cubeY,cubeFX,cubeFY,cubeP,cubeFXs,cubeFYs,cubePs,cubeXD,cubeYD,cubeXDs,cubeYDs,xAv,yAv,fxAv,fxAv,fyAv,pAv,fxAvs,fyAvs,pAvs,rmsVal,rmsX,rmsY,rmsVals,rmsXs,rmsYs=ar.loadCubes(loadPref,source,frameID)
    vmi,vma=ar.bounds(rmsXs,2)
    
    sc1=ax[0,0].scatter(xAv,yAv,c=rmsXs)
    ax[0,0].set_title(str(int(frameID))+" ["+st+"]")
    fig.colorbar(sc1,ax=ax[0,0])

    frameID=frames[1]
    st=ar.loadPhys("Phys/nphys_",frameID)
    cubeX,cubeY,cubeFX,cubeFY,cubeP,cubeFXs,cubeFYs,cubePs,cubeXD,cubeYD,cubeXDs,cubeYDs,xAv,yAv,fxAv,fxAv,fyAv,pAv,fxAvs,fyAvs,pAvs,rmsVal,rmsX,rmsY,rmsVals,rmsXs,rmsYs=ar.loadCubes(loadPref,source,frameID)
    vmi,vma=ar.bounds(rmsXs,2)
    sc1=ax[0,1].scatter(xAv,yAv,c=rmsXs,vmin=vmi,vmax=vma)
    ax[0,1].set_title(str(int(frameID))+" ["+st+"]")
    fig.colorbar(sc1,ax=ax[0,1])

    frameID=frames[2]
    st=ar.loadPhys("Phys/nphys_",frameID)
    cubeX,cubeY,cubeFX,cubeFY,cubeP,cubeFXs,cubeFYs,cubePs,cubeXD,cubeYD,cubeXDs,cubeYDs,xAv,yAv,fxAv,fxAv,fyAv,pAv,fxAvs,fyAvs,pAvs,rmsVal,rmsX,rmsY,rmsVals,rmsXs,rmsYs=ar.loadCubes(loadPref,source,frameID)
    vmi,vma=ar.bounds(rmsXs,2)
    sc1=ax[1,0].scatter(xAv,yAv,c=rmsXs,vmin=vmi,vmax=vma)
    ax[1,0].set_title(str(int(frameID))+" ["+st+"]")
    fig.colorbar(sc1,ax=ax[1,0])

    frameID=frames[3]
    st=ar.loadPhys("Phys/nphys_",frameID)
    cubeX,cubeY,cubeFX,cubeFY,cubeP,cubeFXs,cubeFYs,cubePs,cubeXD,cubeYD,cubeXDs,cubeYDs,xAv,yAv,fxAv,fxAv,fyAv,pAv,fxAvs,fyAvs,pAvs,rmsVal,rmsX,rmsY,rmsVals,rmsXs,rmsYs=ar.loadCubes(loadPref,source,frameID)
    vmi,vma=ar.bounds(rmsXs,2)
    sc1=ax[1,1].scatter(xAv,yAv,c=rmsXs,vmin=vmi,vmax=vma)
    ax[1,1].set_title(str(int(frameID))+" ["+st+"]")
    fig.colorbar(sc1,ax=ax[1,1])

    plt.suptitle("Spot Motion RMS (x) with Rotation")
    plt.savefig("rmsRotX.png")

def rmsTime():

    """
    show the spot motion rms variation as a function of exposure time
    Figure 30
    """

    loadPref="dbSet"
    source="l"

    frames=np.array([1762200,1762300,1762400,1762500]).astype('int')

    fig,ax=plt.subplots(2,2,figsize=(14,10))

    frameID=frames[0]
    st=ar.loadPhys("Phys/nphys_",frameID)
    cubeX,cubeY,cubeFX,cubeFY,cubeP,cubeFXs,cubeFYs,cubePs,cubeXD,cubeYD,cubeXDs,cubeYDs,xAv,yAv,fxAv,fxAv,fyAv,pAv,fxAvs,fyAvs,pAvs,rmsVal,rmsX,rmsY,rmsVals,rmsXs,rmsYs=ar.loadCubes(loadPref,source,frameID)
    vmi,vma=ar.bounds(rmsXs,2)

    sc1=ax[0,0].scatter(xAv,yAv,c=rmsXs,vmin=vmi,vmax=vma)
    ax[0,0].set_title(str(int(frameID))+" ["+st+"]")
    fig.colorbar(sc1,ax=ax[0,0])

    frameID=frames[1]
    st=ar.loadPhys("Phys/nphys_",frameID)
    cubeX,cubeY,cubeFX,cubeFY,cubeP,cubeFXs,cubeFYs,cubePs,cubeXD,cubeYD,cubeXDs,cubeYDs,xAv,yAv,fxAv,fxAv,fyAv,pAv,fxAvs,fyAvs,pAvs,rmsVal,rmsX,rmsY,rmsVals,rmsXs,rmsYs=ar.loadCubes(loadPref,source,frameID)
    vmi,vma=ar.bounds(rmsXs,2)

    sc1=ax[0,1].scatter(xAv,yAv,c=rmsXs,vmin=vmi,vmax=vma)
    ax[0,1].set_title(str(int(frameID))+" ["+st+"]")
    fig.colorbar(sc1,ax=ax[0,1])

    frameID=frames[2]
    st=ar.loadPhys("Phys/nphys_",frameID)
    cubeX,cubeY,cubeFX,cubeFY,cubeP,cubeFXs,cubeFYs,cubePs,cubeXD,cubeYD,cubeXDs,cubeYDs,xAv,yAv,fxAv,fxAv,fyAv,pAv,fxAvs,fyAvs,pAvs,rmsVal,rmsX,rmsY,rmsVals,rmsXs,rmsYs=ar.loadCubes(loadPref,source,frameID)
    vmi,vma=ar.bounds(rmsXs,2)

    sc1=ax[1,0].scatter(xAv,yAv,c=rmsXs,vmin=vmi,vmax=vma)
    ax[1,0].set_title(str(int(frameID))+" ["+st+"]")
    fig.colorbar(sc1,ax=ax[1,0])

    frameID=frames[3]
    st=ar.loadPhys("Phys/nphys_",frameID)
    cubeX,cubeY,cubeFX,cubeFY,cubeP,cubeFXs,cubeFYs,cubePs,cubeXD,cubeYD,cubeXDs,cubeYDs,xAv,yAv,fxAv,fxAv,fyAv,pAv,fxAvs,fyAvs,pAvs,rmsVal,rmsX,rmsY,rmsVals,rmsXs,rmsYs=ar.loadCubes(loadPref,source,frameID)
    vmi,vma=ar.bounds(rmsXs,2)

    sc1=ax[1,1].scatter(xAv,yAv,c=rmsXs,vmin=vmi,vmax=vma)
    ax[1,1].set_title(str(int(frameID))+" ["+st+"]")
    fig.colorbar(sc1,ax=ax[1,1])

    plt.suptitle("Spot Motion RMS (x) with Exposure Time")
    plt.savefig("rmsTime.png")

def positions():

    """
    A simple plot showing spot locations
    Used in POwerpoint
    """

    frameID=15874
    loadPref="dbSet"
    source="l"

    fig,ax=plt.subplots()
    st=ar.loadPhys("Phys/nphys_",frameID)
    cubeX,cubeY,cubeFX,cubeFY,cubeP,cubeFXs,cubeFYs,cubePs,cubeXD,cubeYD,cubeXDs,cubeYDs,xAv,yAv,fxAv,fxAv,fyAv,pAv,fxAvs,fyAvs,pAvs,rmsVal,rmsX,rmsY,rmsVals,rmsXs,rmsYs=ar.loadCubes(loadPref,source,frameID)

    sc1=ax.scatter(xAv,yAv)
    ax.set_xlabel("Position (mm)")
    ax.set_ylabel("Position (mm)")
    plt.savefig("positions.png")


def smoothCpt1():

    """
    show the smooth variation of parameters across the field, 
    for different rotations, showing a set of parameters on each plot.


    Figure 6

    """
    
    frames=np.array([15854,15864,15874,15884]).astype('int')

    loadPref="dbSet"
    source="w"

    fig,ax=plt.subplots(2,2,figsize=(14,10))


    frameID=frames[2]
    st=ar.loadPhys("Phys/nphys_",frameID)
    cubeX,cubeY,cubeFX,cubeFY,cubeP,cubeFXs,cubeFYs,cubePs,cubeXD,cubeYD,cubeXDs,cubeYDs,xAv,yAv,fxAv,fxAv,fyAv,pAv,fxAvs,fyAvs,pAvs,rmsVal,rmsX,rmsY,rmsVals,rmsXs,rmsYs=ar.loadCubes(loadPref,source,frameID)

    scpt=fxAv-fxAvs+fxAv.mean()
    vmi,vma=ar.bounds(scpt,2)
    sc1=ax[0,0].scatter(xAv,yAv,c=scpt,vmin=vmi,vmax=vma)
    ax[0,0].set_title("Spot Size (x)")
    fig.colorbar(sc1,ax=ax[0,0])

    scpt=fyAv-fyAvs+fyAv.mean()
    vmi,vma=ar.bounds(scpt,2)
    sc1=ax[0,1].scatter(xAv,yAv,c=scpt,vmin=vmi,vmax=vma)
    ax[0,1].set_title("Spot Size (y)")
    fig.colorbar(sc1,ax=ax[0,1])

    scpt=pAv-pAvs+pAv.mean()
    vmi,vma=ar.bounds(scpt,2)
    sc1=ax[1,0].scatter(xAv,yAv,c=scpt,vmin=vmi,vmax=vma)
    ax[1,0].set_title("Spot Brightness")
    fig.colorbar(sc1,ax=ax[1,0])
    fig.delaxes(ax[1,1])


    plt.suptitle("Variations of Paramters across Field ("+str(int(frameID))+") ["+st+"]")
    plt.savefig("smoothSummary.png")


def smoothCpt():

    """
    show the smooth variation of parameters across the field, for different rotations, 
    showing a set of different rotations 
    for a single parameter in each plot. 
    Not used in report
    """

    
    frames=np.array([15854,15864,15874,15884]).astype('int')

    #frames=np.array([1762200,1762300,1762400,1762500]).astype('int')
    loadPref="dbSet"
    source="w"

    fig,ax=plt.subplots(2,2,figsize=(14,10))

    frameID=frames[0]
    st=ar.loadPhys("Phys/nphys_",frameID)
    cubeX,cubeY,cubeFX,cubeFY,cubeP,cubeFXs,cubeFYs,cubePs,cubeXD,cubeYD,cubeXDs,cubeYDs,xAv,yAv,fxAv,fxAv,fyAv,pAv,fxAvs,fyAvs,pAvs,rmsVal,rmsX,rmsY,rmsVals,rmsXs,rmsYs=ar.loadCubes(loadPref,source,frameID)

    scpt=fxAv-fxAvs+fxAv.mean()
    vmi,vma=ar.bounds(scpt,2)
    sc1=ax[0,0].scatter(xAv,yAv,c=scpt,vmin=vmi,vmax=vma)
    ax[0,0].set_title(str(int(frameID))+" ["+st+"]")
    fig.colorbar(sc1,ax=ax[0,0])

    frameID=frames[1]
    st=ar.loadPhys("Phys/nphys_",frameID)
    cubeX,cubeY,cubeFX,cubeFY,cubeP,cubeFXs,cubeFYs,cubePs,cubeXD,cubeYD,cubeXDs,cubeYDs,xAv,yAv,fxAv,fxAv,fyAv,pAv,fxAvs,fyAvs,pAvs,rmsVal,rmsX,rmsY,rmsVals,rmsXs,rmsYs=ar.loadCubes(loadPref,source,frameID)

    scpt=fxAv-fxAvs+fxAv.mean()
    vmi,vma=ar.bounds(scpt,2)
    sc1=ax[0,1].scatter(xAv,yAv,c=scpt,vmin=vmi,vmax=vma)
    ax[0,1].set_title(str(int(frameID))+" ["+st+"]")
    fig.colorbar(sc1,ax=ax[0,1])

    frameID=frames[2]
    st=ar.loadPhys("Phys/nphys_",frameID)
    cubeX,cubeY,cubeFX,cubeFY,cubeP,cubeFXs,cubeFYs,cubePs,cubeXD,cubeYD,cubeXDs,cubeYDs,xAv,yAv,fxAv,fxAv,fyAv,pAv,fxAvs,fyAvs,pAvs,rmsVal,rmsX,rmsY,rmsVals,rmsXs,rmsYs=ar.loadCubes(loadPref,source,frameID)

    scpt=fxAv-fxAvs+fxAv.mean()
    vmi,vma=ar.bounds(scpt,2)
    sc1=ax[1,0].scatter(xAv,yAv,c=scpt,vmin=vmi,vmax=vma)
    ax[1,0].set_title(str(int(frameID))+" ["+st+"]")
    fig.colorbar(sc1,ax=ax[1,0])

    frameID=frames[3]
    st=ar.loadPhys("Phys/nphys_",frameID)
    cubeX,cubeY,cubeFX,cubeFY,cubeP,cubeFXs,cubeFYs,cubePs,cubeXD,cubeYD,cubeXDs,cubeYDs,xAv,yAv,fxAv,fxAv,fyAv,pAv,fxAvs,fyAvs,pAvs,rmsVal,rmsX,rmsY,rmsVals,rmsXs,rmsYs=ar.loadCubes(loadPref,source,frameID)

    scpt=fxAv-fxAvs+fxAv.mean()
    vmi,vma=ar.bounds(scpt,2)
    sc1=ax[1,1].scatter(xAv,yAv,c=scpt,vmin=vmi,vmax=vma)
    ax[1,1].set_title(str(int(frameID))+" ["+st+"]")
    fig.colorbar(sc1,ax=ax[1,1])


    fig.suptitle("Smooth component of Average Spot Size (x)")
    fig.savefig("smoothFX.png")
    fig,ax=plt.subplots(2,2,figsize=(14,10))

    frameID=frames[0]
    st=ar.loadPhys("Phys/nphys_",frameID)
    cubeX,cubeY,cubeFX,cubeFY,cubeP,cubeFXs,cubeFYs,cubePs,cubeXD,cubeYD,cubeXDs,cubeYDs,xAv,yAv,fxAv,fxAv,fyAv,pAv,fxAvs,fyAvs,pAvs,rmsVal,rmsX,rmsY,rmsVals,rmsXs,rmsYs=ar.loadCubes(loadPref,source,frameID)

    scpt=fyAv-fyAvs+fyAv.mean()
    vmi,vma=ar.bounds(scpt,2)
    sc1=ax[0,0].scatter(xAv,yAv,c=scpt,vmin=vmi,vmax=vma)
    ax[0,0].set_title(str(int(frameID))+" ["+st+"]")
    fig.colorbar(sc1,ax=ax[0,0])

    frameID=frames[1]
    st=ar.loadPhys("Phys/nphys_",frameID)
    cubeX,cubeY,cubeFX,cubeFY,cubeP,cubeFXs,cubeFYs,cubePs,cubeXD,cubeYD,cubeXDs,cubeYDs,xAv,yAv,fxAv,fxAv,fyAv,pAv,fxAvs,fyAvs,pAvs,rmsVal,rmsX,rmsY,rmsVals,rmsXs,rmsYs=ar.loadCubes(loadPref,source,frameID)

    scpt=fyAv-fyAvs+fyAv.mean()
    vmi,vma=ar.bounds(scpt,2)
    sc1=ax[0,1].scatter(xAv,yAv,c=scpt,vmin=vmi,vmax=vma)
    ax[0,1].set_title(str(int(frameID))+" ["+st+"]")
    fig.colorbar(sc1,ax=ax[0,1])

    frameID=frames[2]
    st=ar.loadPhys("Phys/nphys_",frameID)
    cubeX,cubeY,cubeFX,cubeFY,cubeP,cubeFXs,cubeFYs,cubePs,cubeXD,cubeYD,cubeXDs,cubeYDs,xAv,yAv,fxAv,fxAv,fyAv,pAv,fxAvs,fyAvs,pAvs,rmsVal,rmsX,rmsY,rmsVals,rmsXs,rmsYs=ar.loadCubes(loadPref,source,frameID)

    scpt=fyAv-fyAvs+fyAv.mean()
    vmi,vma=ar.bounds(scpt,2)
    sc1=ax[1,0].scatter(xAv,yAv,c=scpt,vmin=vmi,vmax=vma)
    ax[1,0].set_title(str(int(frameID))+" ["+st+"]")
    fig.colorbar(sc1,ax=ax[1,0])

    frameID=frames[3]
    st=ar.loadPhys("Phys/nphys_",frameID)
    cubeX,cubeY,cubeFX,cubeFY,cubeP,cubeFXs,cubeFYs,cubePs,cubeXD,cubeYD,cubeXDs,cubeYDs,xAv,yAv,fxAv,fxAv,fyAv,pAv,fxAvs,fyAvs,pAvs,rmsVal,rmsX,rmsY,rmsVals,rmsXs,rmsYs=ar.loadCubes(loadPref,source,frameID)

    scpt=fyAv-fyAvs+fyAv.mean()
    vmi,vma=ar.bounds(scpt,2)
    sc1=ax[1,1].scatter(xAv,yAv,c=scpt,vmin=vmi,vmax=vma)
    ax[1,1].set_title(str(int(frameID))+" ["+st+"]")
    fig.colorbar(sc1,ax=ax[1,1])


    fig.suptitle("Smooth component of Average Spot Size (y)")
    fig.savefig("smoothFY.png")

    
    fig,ax=plt.subplots(2,2,figsize=(14,10))

    frameID=frames[0]
    st=ar.loadPhys("Phys/nphys_",frameID)
    cubeX,cubeY,cubeFX,cubeFY,cubeP,cubeFXs,cubeFYs,cubePs,cubeXD,cubeYD,cubeXDs,cubeYDs,xAv,yAv,fxAv,fxAv,fyAv,pAv,fxAvs,fyAvs,pAvs,rmsVal,rmsX,rmsY,rmsVals,rmsXs,rmsYs=ar.loadCubes(loadPref,source,frameID)

    scpt=pAv-pAvs+pAv.mean()
    vmi,vma=ar.bounds(scpt,2)
    sc1=ax[0,0].scatter(xAv,yAv,c=scpt,vmin=vmi,vmax=vma)
    ax[0,0].set_title(str(int(frameID))+" ["+st+"]")
    fig.colorbar(sc1,ax=ax[0,0])

    frameID=frames[1]
    st=ar.loadPhys("Phys/nphys_",frameID)
    cubeX,cubeY,cubeFX,cubeFY,cubeP,cubeFXs,cubeFYs,cubePs,cubeXD,cubeYD,cubeXDs,cubeYDs,xAv,yAv,fxAv,fxAv,fyAv,pAv,fxAvs,fyAvs,pAvs,rmsVal,rmsX,rmsY,rmsVals,rmsXs,rmsYs=ar.loadCubes(loadPref,source,frameID)

    scpt=pAv-pAvs+pAv.mean()
    vmi,vma=ar.bounds(scpt,2)
    sc1=ax[0,1].scatter(xAv,yAv,c=scpt,vmin=vmi,vmax=vma)
    ax[0,1].set_title(str(int(frameID))+" ["+st+"]")
    fig.colorbar(sc1,ax=ax[0,1])

    frameID=frames[2]
    st=ar.loadPhys("Phys/nphys_",frameID)
    cubeX,cubeY,cubeFX,cubeFY,cubeP,cubeFXs,cubeFYs,cubePs,cubeXD,cubeYD,cubeXDs,cubeYDs,xAv,yAv,fxAv,fxAv,fyAv,pAv,fxAvs,fyAvs,pAvs,rmsVal,rmsX,rmsY,rmsVals,rmsXs,rmsYs=ar.loadCubes(loadPref,source,frameID)

    scpt=pAv-pAvs+pAv.mean()
    vmi,vma=ar.bounds(scpt,2)
    sc1=ax[1,0].scatter(xAv,yAv,c=scpt,vmin=vmi,vmax=vma)
    ax[1,0].set_title(str(int(frameID))+" ["+st+"]")
    fig.colorbar(sc1,ax=ax[1,0])

    frameID=frames[3]
    st=ar.loadPhys("Phys/nphys_",frameID)
    cubeX,cubeY,cubeFX,cubeFY,cubeP,cubeFXs,cubeFYs,cubePs,cubeXD,cubeYD,cubeXDs,cubeYDs,xAv,yAv,fxAv,fxAv,fyAv,pAv,fxAvs,fyAvs,pAvs,rmsVal,rmsX,rmsY,rmsVals,rmsXs,rmsYs=ar.loadCubes(loadPref,source,frameID)

    scpt=pAv-pAvs+pAv.mean()
    vmi,vma=ar.bounds(scpt,2)
    sc1=ax[1,1].scatter(xAv,yAv,c=scpt,vmin=vmi,vmax=vma)
    ax[1,1].set_title(str(int(frameID))+" ["+st+"]")
    fig.colorbar(sc1,ax=ax[1,1])


    fig.suptitle("Smooth component of Average Peak Value")
    fig.savefig("smoothP.png")
    
def methodsRMS(frameID):

    """
    show maps of the RMS for differnet methods 
    figure 45,46
    """

    loadPref="dbSet"
    source1="w"
    source2="l"
    source3="p"
    st=ar.loadPhys("Phys/nphys_",frameID)
    cubeX1,cubeY1,cubeFX1,cubeFY1,cubeP1,cubeFXs1,cubeFYs1,cubePs1,cubeXD1,cubeYD1,cubeXDs1,cubeYDs1,xAv1,yAv1,fxAv1,fxAv1,fyAv1,pAv1,fxAvs1,fyAvs1,pAvs1,rmsVal1,rmsX1,rmsY1,rmsVals1,rmsXs1,rmsYs1=ar.loadCubes(loadPref,source1,frameID)


    cubeX2,cubeY2,cubeFX2,cubeFY2,cubeP2,cubeFXs2,cubeFYs2,cubePs2,cubeXD2,cubeYD2,cubeXDs2,cubeYDs2,xAv2,yAv2,fxAv2,fxAv2,fyAv2,pAv2,fxAvs2,fyAvs2,pAvs2,rmsVal2,rmsX2,rmsY2,rmsVals2,rmsXs2,rmsYs2=ar.loadCubes(loadPref,source2,frameID)

    cubeX3,cubeY3,cubeFX3,cubeFY3,cubeP3,cubeFXs3,cubeFYs3,cubePs3,cubeXD3,cubeYD3,cubeXDs3,cubeYDs3,xAv3,yAv3,fxAv3,fxAv3,fyAv3,pAv3,fxAvs3,fyAvs3,pAvs3,rmsVal3,rmsX3,rmsY3,rmsVals3,rmsXs3,rmsYs3=ar.loadCubes(loadPref,source3,frameID)

    
    xdiffwl=cubeX1-cubeX2
    ydiffwl=cubeY1-cubeY2

    xdiffpl=cubeX3-cubeX2
    ydiffpl=cubeY3-cubeY2
    

    xdrwl=xdiffwl.std(axis=1)
    ydrwl=ydiffwl.std(axis=1)
    xdrpl=xdiffpl.std(axis=1)
    ydrpl=ydiffpl.std(axis=1)

    fig,ax=plt.subplots(2,2,figsize=(14,10))

    vmi,vma=ar.bounds(xdrwl,1.5)
    sc1=ax[0,0].scatter(xAv1,yAv1,c=xdrwl,vmin=vmi,vmax=vma)
    fig.colorbar(sc1,ax=ax[0,0])
    ax[0,0].set_title("RMS Windowed-MCS (X)")

    vmi,vma=ar.bounds(ydrwl,1.5)
    sc1=ax[0,1].scatter(xAv1,yAv1,c=ydrwl,vmin=vmi,vmax=vma)
    fig.colorbar(sc1,ax=ax[0,1])
    ax[0,1].set_title("RMS Windowed-MCS (Y)")

    vmi,vma=ar.bounds(xdrpl,1.5)
    sc1=ax[1,0].scatter(xAv1,yAv1,c=xdrpl,vmin=vmi,vmax=vma)
    fig.colorbar(sc1,ax=ax[1,0])
    ax[1,0].set_title("RMS Sextractor-MCS (X)")

    vmi,vma=ar.bounds(ydrpl,1.5)
    sc1=ax[1,1].scatter(xAv1,yAv1,c=ydrpl,vmin=vmi,vmax=vma)
    fig.colorbar(sc1,ax=ax[1,1])
    ax[1,1].set_title("RMS Sextractor-MCS (Y)")

    plt.suptitle("Direct Comparison Between Methods for Set "+str(int(frameID))+" ["+st+"]")
    plt.savefig("methodRMS_"+str(int(frameID))+".png")
    
def methods():

    """
    show maps of the spot sizes for different methods
    figure 32-34
    """
    
    frameID=15854
    loadPref="dbSet"
    st=ar.loadPhys("Phys/nphys_",frameID)

    offset=3
    source="p"

    cubeX,cubeY,cubeFX,cubeFY,cubeP,cubeFXs,cubeFYs,cubePs,cubeXD,cubeYD,cubeXDs,cubeYDs,xAv,yAv,fxAv,fxAv,fyAv,pAv,fxAvs,fyAvs,pAvs,rmsVal,rmsX,rmsY,rmsVals,rmsXs,rmsYs=ar.loadCubes(loadPref,source,frameID)

    pFX=cubeFX[:,offset]
    pFY=cubeFY[:,offset]
    pFX=cubeFXs[:,offset]
    pFY=cubeFYs[:,offset]
    pXD=cubeXD[:,offset]
    pYD=cubeYD[:,offset]
    pXD=cubeXDs[:,offset]
    pYD=cubeYDs[:,offset]
    prmsX=rmsX
    
    source="w"

    cubeX,cubeY,cubeFX,cubeFY,cubeP,cubeFXs,cubeFYs,cubePs,cubeXD,cubeYD,cubeXDs,cubeYDs,xAv,yAv,fxAv,fxAv,fyAv,pAv,fxAvs,fyAvs,pAvs,rmsVal,rmsX,rmsY,rmsVals,rmsXs,rmsYs=ar.loadCubes(loadPref,source,frameID)

    wFX=cubeFX[:,offset]
    wFY=cubeFY[:,offset]
    wFX=cubeFXs[:,offset]
    wFY=cubeFYs[:,offset]
    wXD=cubeXD[:,offset]
    wYD=cubeYD[:,offset]
    wXD=cubeXDs[:,offset]
    wYD=cubeYDs[:,offset]
    wrmsX=rmsX
    
    source="l"

    cubeX,cubeY,cubeFX,cubeFY,cubeP,cubeFXs,cubeFYs,cubePs,cubeXD,cubeYD,cubeXDs,cubeYDs,xAv,yAv,fxAv,fxAv,fyAv,pAv,fxAvs,fyAvs,pAvs,rmsVal,rmsX,rmsY,rmsVals,rmsXs,rmsYs=ar.loadCubes(loadPref,source,frameID)

    lFX=cubeFX[:,offset]
    lFY=cubeFY[:,offset]
    lFX=cubeFXs[:,offset]
    lFY=cubeFYs[:,offset]
    lXD=cubeXD[:,offset]
    lYD=cubeYD[:,offset]
    lXD=cubeXDs[:,offset]
    lYD=cubeYDs[:,offset]
    lrmsX=rmsX
    
    source="t"

    cubeX,cubeY,cubeFX,cubeFY,cubeP,cubeFXs,cubeFYs,cubePs,cubeXD,cubeYD,cubeXDs,cubeYDs,xAv,yAv,fxAv,fxAv,fyAv,pAv,fxAvs,fyAvs,pAvs,rmsVal,rmsX,rmsY,rmsVals,rmsXs,rmsYs=ar.loadCubes(loadPref,source,frameID)

    mFX=cubeFX[:,offset]
    mFY=cubeFY[:,offset]
    mFX=cubeFXs[:,offset]
    mFY=cubeFYs[:,offset]
    mXD=cubeXD[:,offset]
    mYD=cubeYD[:,offset]
    mXD=cubeXDs[:,offset]
    mYD=cubeYDs[:,offset]
    mrmsX=rmsX
    
    fig,ax=plt.subplots(2,2,figsize=(14,10))
    vmi,vma=ar.bounds(pFX,1.5)
    #sc1=ax[0,0].scatter(xAv,yAv,c=pFX,vmin=vmi,vmax=vma)
    #fig.colorbar(sc1,ax=ax[0,0])
    #ax[0,0].set_title("PSF Fitting")

    vmi,vma=ar.bounds(wFX,1.5)
    sc1=ax[0,1].scatter(xAv,yAv,c=wFX,vmin=vmi,vmax=vma)
    fig.colorbar(sc1,ax=ax[0,1])
    ax[0,1].set_title("SEXtractor Windowed")

    vmi,vma=ar.bounds(lFX,1.5)
    sc1=ax[1,0].scatter(xAv,yAv,c=lFX,vmin=vmi,vmax=vma)
    fig.colorbar(sc1,ax=ax[1,0])
    ax[1,0].set_title("MCS Implemented")

    #vmi,vma=ar.bounds(mFX,1.5)
    #sc1=ax[1,1].scatter(xAv,yAv,c=mFX,vmin=vmi,vmax=vma)
    ##fig.colorbar(sc1,ax=ax[1,1])
    #ax[1,1].set_title("Moment")

    fig.delaxes(ax[1,1])
    fig.delaxes(ax[0,0])
    
    plt.suptitle("Spot Size (x) with different methods ("+str(frameID)+") ["+st+"]")
    plt.savefig("methodFX.png")

    #-------------------------------------------
    
    fig,ax=plt.subplots(2,2,figsize=(14,10))
    vmi,vma=ar.bounds(pXD,1.5)
    sc1=ax[0,0].scatter(xAv,yAv,c=pXD,vmin=vmi,vmax=vma)
    fig.colorbar(sc1,ax=ax[0,0])
    ax[0,0].set_title("PSF Fitting")

    vmi,vma=ar.bounds(wXD,1.5)
    sc1=ax[0,1].scatter(xAv,yAv,c=wXD,vmin=vmi,vmax=vma)
    fig.colorbar(sc1,ax=ax[0,1])
    ax[0,1].set_title("SEXtractor Windowed")

    vmi,vma=ar.bounds(lXD,1.5)
    sc1=ax[1,0].scatter(xAv,yAv,c=lXD,vmin=vmi,vmax=vma)
    fig.colorbar(sc1,ax=ax[1,0])
    ax[1,0].set_title("MCS Implemented")

    vmi,vma=ar.bounds(mXD,1.5)
    sc1=ax[1,1].scatter(xAv,yAv,c=mXD,vmin=vmi,vmax=vma)
    fig.colorbar(sc1,ax=ax[1,1])
    ax[1,1].set_title("Moment")

    plt.suptitle("dx with different methods ("+str(frameID)+") ["+st+"]")
    plt.savefig("methodDX.png")
    
    #-------------------------------------------
    
    fig,ax=plt.subplots(2,2,figsize=(14,10))
    vmi,vma=ar.bounds(prmsX,1.5)
    sc1=ax[0,0].scatter(xAv,yAv,c=prmsX,vmin=vmi,vmax=vma)
    fig.colorbar(sc1,ax=ax[0,0])
    ax[0,0].set_title("PSF Fitting")

    vmi,vma=ar.bounds(wrmsX,1.5)
    sc1=ax[0,1].scatter(xAv,yAv,c=wrmsX,vmin=vmi,vmax=vma)
    fig.colorbar(sc1,ax=ax[0,1])
    ax[0,1].set_title("SEXtractor Windowed")

    vmi,vma=ar.bounds(lrmsX,1.5)
    sc1=ax[1,0].scatter(xAv,yAv,c=lrmsX,vmin=vmi,vmax=vma)
    fig.colorbar(sc1,ax=ax[1,0])
    ax[1,0].set_title("MCS Implemented")

    vmi,vma=ar.bounds(mrmsX,1.5)
    sc1=ax[1,1].scatter(xAv,yAv,c=mrmsX,vmin=vmi,vmax=vma)
    fig.colorbar(sc1,ax=ax[1,1])
    ax[1,1].set_title("Moment")

    plt.suptitle("RMS with different methods ("+str(frameID)+") ["+st+"]")
    plt.savefig("methodRMS.png")
    
        
def exampleRMS2():

    wdata = np.genfromtxt("w_rmsC.dat", delimiter=None, names=True,dtype=None)
    pdata = np.genfromtxt("p_rmsC.dat", delimiter=None, names=True,dtype=None)
    tdata = np.genfromtxt("t_rmsC.dat", delimiter=None, names=True,dtype=None)
    ldata = np.genfromtxt("l_rmsC.dat", delimiter=None, names=True,dtype=None)

    xx=np.arange(0,0.005,0.0005)

    ind1=np.where(wdata['t']==0.5)
    ind2=np.where(wdata['t']==1)
    ind3=np.where(wdata['t']==2)
    ind4=np.where(wdata['t']==5)

    fig,ax=plt.subplots(figsize=(8,4))

    ax.scatter(pdata['rmsmm'][ind1],wdata['rmsmm'][ind1],color="red")
    ax.scatter(pdata['rmsmm'][ind3],wdata['rmsmm'][ind3],color="yellow")
    ax.scatter(pdata['rmsmm'][ind4],wdata['rmsmm'][ind4],color="green")
    ax.scatter(pdata['rmsmm'][ind2],wdata['rmsmm'][ind2],color="orange")
    ax.plot(xx,xx,color="black")
    ax.set_xlim(0,0.005)
    ax.set_ylim(0,0.0055)
    ax.set_ylabel("RMS Spot Motion")
    ax.set_title("SeXTractor Windowed")

    plt.savefig("eg1.png")
    fig,ax=plt.subplots(figsize=(8,4))

    ax.scatter(pdata['rmsmm'][ind1],ldata['rmsmm'][ind1],color="red")
    ax.scatter(pdata['rmsmm'][ind3],ldata['rmsmm'][ind3],color="yellow")
    ax.scatter(pdata['rmsmm'][ind4],ldata['rmsmm'][ind4],color="green")
    ax.scatter(pdata['rmsmm'][ind2],ldata['rmsmm'][ind2],color="orange")
    ax.plot(xx,xx,color="black")
    ax.set_xlim(0.0005,0.005)
    ax.set_ylim(0.0005,0.0055)
    ax.set_ylabel("RMS Spot Motion")
    ax.set_title("MCS Code")

    plt.savefig("eg2.png")

def exampleRMS():

    """
     Comparing RMS between different methods over the whole set. 
    Figures 39,40
    """

    wdata = np.genfromtxt("w_rmsC.dat", delimiter=None, names=True,dtype=None)
    pdata = np.genfromtxt("p_rmsC.dat", delimiter=None, names=True,dtype=None)
    tdata = np.genfromtxt("t_rmsC.dat", delimiter=None, names=True,dtype=None)
    ldata = np.genfromtxt("l_rmsC.dat", delimiter=None, names=True,dtype=None)

    xx=np.arange(0,0.005,0.0005)

    ind1=np.where(wdata['t']==0.5)
    ind2=np.where(wdata['t']==1)
    ind3=np.where(wdata['t']==2)
    ind4=np.where(wdata['t']==5)

    fig,ax=plt.subplots(3,figsize=(13,13))

    ax[0].scatter(pdata['rmsmm'][ind1],wdata['rmsmm'][ind1],color="red")
    ax[0].scatter(pdata['rmsmm'][ind3],wdata['rmsmm'][ind3],color="yellow")
    ax[0].scatter(pdata['rmsmm'][ind4],wdata['rmsmm'][ind4],color="green")
    ax[0].scatter(pdata['rmsmm'][ind2],wdata['rmsmm'][ind2],color="orange")
    ax[0].plot(xx,xx,color="black")
    ax[0].set_xlim(0,0.005)
    ax[0].set_ylim(0,0.0055)
    ax[0].set_ylabel("RMS Spot Motion")
    ax[0].set_title("SeXTractor Windowed")


    ax[1].scatter(pdata['rmsmm'][ind1],ldata['rmsmm'][ind1],color="red")
    ax[1].scatter(pdata['rmsmm'][ind3],ldata['rmsmm'][ind3],color="yellow")
    ax[1].scatter(pdata['rmsmm'][ind4],ldata['rmsmm'][ind4],color="green")
    ax[1].scatter(pdata['rmsmm'][ind2],ldata['rmsmm'][ind2],color="orange")
    ax[1].plot(xx,xx,color="black")
    ax[1].set_xlim(0.0005,0.005)
    ax[1].set_ylim(0.0005,0.0055)
    ax[1].set_ylabel("RMS Spot Motion")
    ax[1].set_title("MCS Code")

    ax[2].scatter(pdata['rmsmm'][ind1],tdata['rmsmm'][ind1],color="red")
    ax[2].scatter(pdata['rmsmm'][ind3],tdata['rmsmm'][ind3],color="yellow")
    ax[2].scatter(pdata['rmsmm'][ind4],tdata['rmsmm'][ind4],color="green")
    ax[2].scatter(pdata['rmsmm'][ind2],tdata['rmsmm'][ind2],color="orange")
    ax[2].plot(xx,xx,color="black")
    ax[2].set_xlim(0.0005,0.005)
    ax[2].set_ylim(0.0005,0.0055)
    ax[2].set_ylabel("RMS Spot Motion")
    ax[2].set_title("Moment")
    ax[2].set_xlabel("RMS Spot Motion (PSF Fitting)")
    
    plt.suptitle("Comparison of Average RMS Values with PSF Fitting")
    plt.savefig("rmsComp1.png")

    fig,ax=plt.subplots(3,figsize=(13,13))

    ax[0].scatter(pdata['rmsmm'][ind1],wdata['rmsmm'][ind1]-pdata['rmsmm'][ind1],color="red")
    ax[0].scatter(pdata['rmsmm'][ind3],wdata['rmsmm'][ind3]-pdata['rmsmm'][ind3],color="yellow")
    ax[0].scatter(pdata['rmsmm'][ind4],wdata['rmsmm'][ind4]-pdata['rmsmm'][ind4],color="green")
    ax[0].scatter(pdata['rmsmm'][ind2],wdata['rmsmm'][ind2]-pdata['rmsmm'][ind2],color="orange")
    ax[0].plot(xx,np.zeros(len(xx)),color="black")
    ax[0].set_xlim(0,0.005)
    ax[0].set_ylim(-0.001,0.001)
    ax[0].set_ylabel("RMS Spot Motion Differences")
    ax[0].set_title("SeXTractor Windowed")

    ax[1].scatter(pdata['rmsmm'][ind1],ldata['rmsmm'][ind1]-pdata['rmsmm'][ind1],color="red")
    ax[1].scatter(pdata['rmsmm'][ind3],ldata['rmsmm'][ind3]-pdata['rmsmm'][ind3],color="yellow")
    ax[1].scatter(pdata['rmsmm'][ind4],ldata['rmsmm'][ind4]-pdata['rmsmm'][ind4],color="green")
    ax[1].scatter(pdata['rmsmm'][ind2],ldata['rmsmm'][ind2]-pdata['rmsmm'][ind2],color="orange")
    ax[1].plot(xx,np.zeros(len(xx)),color="black")
    ax[1].set_xlim(0.0005,0.005)
    ax[1].set_ylim(-0.001,0.001)
    ax[1].set_ylabel("F")
    ax[1].set_ylabel("RMS Spot Motion Differences")
    ax[1].set_title("SeXTractor Windowed")
    ax[1].set_title("MCS Code")

    ax[2].scatter(pdata['rmsmm'][ind1],tdata['rmsmm'][ind1]-pdata['rmsmm'][ind1],color="red")
    ax[2].scatter(pdata['rmsmm'][ind3],tdata['rmsmm'][ind3]-pdata['rmsmm'][ind3],color="yellow")
    ax[2].scatter(pdata['rmsmm'][ind4],tdata['rmsmm'][ind4]-pdata['rmsmm'][ind4],color="green")
    ax[2].scatter(pdata['rmsmm'][ind2],tdata['rmsmm'][ind2]-pdata['rmsmm'][ind2],color="orange")
    ax[2].plot(xx,np.zeros(len(xx)),color="black")
    ax[2].set_xlim(0.0005,0.005)
    ax[2].set_ylim(-0.001,0.001)
    ax[2].set_title("Moment")
    ax[2].set_ylabel("RMS Spot Motion Differences")
    ax[2].set_xlabel("RMS Spot Motion (PSF Fitting)")

    plt.suptitle("Comparison of Difference in Average RMS Between PSF and other Methods")
    plt.savefig("rmsComp2.png")

def exampleRMS1():

    """
     Comparing RMS between different methods over the whole set. 
    Figures 41
    """

    
    wdata = np.genfromtxt("w_rmsC.dat", delimiter=None, names=True,dtype=None)
    pdata = np.genfromtxt("p_rmsC.dat", delimiter=None, names=True,dtype=None)
    ldata = np.genfromtxt("l_rmsC.dat", delimiter=None, names=True,dtype=None)
    tdata = np.genfromtxt("t_rmsC.dat", delimiter=None, names=True,dtype=None)

    xx=np.arange(0,0.005,0.0005)

    ind1=np.where(np.all([wdata['t']==0.5,wdata['day'] > 1],axis=0))
    ind2=np.where(np.all([wdata['t']==1,wdata['day'] > 1],axis=0))
    ind3=np.where(np.all([wdata['t']==2,wdata['day'] > 1],axis=0))
    ind4=np.where(np.all([wdata['t']==5,wdata['day'] > 1],axis=0))

    fig,ax=plt.subplots(3,figsize=(13,13))

    ax[0].scatter(ldata['fx'][ind1],wdata['rmsmm'][ind1],color="red")
    ax[0].scatter(ldata['fx'][ind3],wdata['rmsmm'][ind3],color="yellow")
    ax[0].scatter(ldata['fx'][ind4],wdata['rmsmm'][ind4],color="green")
    ax[0].scatter(ldata['fx'][ind2],wdata['rmsmm'][ind2],color="orange")
    ax[0].set_ylim(0,0.0055)
    ax[0].set_ylabel("RMS Spot Motion")
    ax[0].set_title("SeXTractor Windowed")

    ax[1].scatter(ldata['fx'][ind1],ldata['rmsmm'][ind1],color="red")
    ax[1].scatter(ldata['fx'][ind3],ldata['rmsmm'][ind3],color="yellow")
    ax[1].scatter(ldata['fx'][ind4],ldata['rmsmm'][ind4],color="green")
    ax[1].scatter(ldata['fx'][ind2],ldata['rmsmm'][ind2],color="orange")
    ax[1].set_ylim(0.0005,0.0055)
    ax[1].set_ylabel("RMS Spot Motion")
    ax[1].set_title("MCS")

    ax[2].scatter(ldata['fx'][ind1],tdata['rmsmm'][ind1],color="red")
    ax[2].scatter(ldata['fx'][ind3],tdata['rmsmm'][ind3],color="yellow")
    ax[2].scatter(ldata['fx'][ind4],tdata['rmsmm'][ind4],color="green")
    ax[2].scatter(ldata['fx'][ind2],tdata['rmsmm'][ind2],color="orange")
    ax[2].set_ylim(0.0005,0.0055)
    ax[2].set_ylabel("RMS Spot Motion")
    ax[2].set_title("Moment")
    ax[2].set_xlabel("RMS Spot Motion (PSF Fitting)")

    plt.suptitle("RMS as a function of Spot Size (x) for Different Methods")
    plt.savefig("rmsComp3.png")

    fig,ax=plt.subplots(3,figsize=(13,13))

    ax[0].scatter(ldata['fx'][ind1],wdata['rmsmm'][ind1]-pdata['rmsmm'][ind1],color="red")
    ax[0].scatter(ldata['fx'][ind3],wdata['rmsmm'][ind3]-pdata['rmsmm'][ind3],color="yellow")
    ax[0].scatter(ldata['fx'][ind4],wdata['rmsmm'][ind4]-pdata['rmsmm'][ind4],color="green")
    ax[0].scatter(ldata['fx'][ind2],wdata['rmsmm'][ind2]-pdata['rmsmm'][ind2],color="orange")
    #ax[0].plot(xx,np.zeros(len(xx)),color="black")
    #ax[0].set_xlim(0,0.005)
    ax[0].set_ylim(-0.001,0.001)
    ax[0].set_ylabel("RMS Spot Motion")
    ax[0].set_title("SeXTractor Windowed")

    ax[1].scatter(ldata['fx'][ind1],ldata['rmsmm'][ind1]-pdata['rmsmm'][ind1],color="red")
    ax[1].scatter(ldata['fx'][ind3],ldata['rmsmm'][ind3]-pdata['rmsmm'][ind3],color="yellow")
    ax[1].scatter(ldata['fx'][ind4],ldata['rmsmm'][ind4]-pdata['rmsmm'][ind4],color="green")
    ax[1].scatter(ldata['fx'][ind2],ldata['rmsmm'][ind2]-pdata['rmsmm'][ind2],color="orange")
    #ax[1].plot(xx,np.zeros(len(xx)),color="black")
    #ax[1].set_xlim(0.0005,0.005)
    ax[1].set_ylim(-0.001,0.001)
    ax[1].set_ylabel("RMS Spot Motion")
    ax[1].set_title("MCS")

    ax[2].scatter(ldata['fx'][ind1],tdata['rmsmm'][ind1]-pdata['rmsmm'][ind1],color="red")
    ax[2].scatter(ldata['fx'][ind3],tdata['rmsmm'][ind3]-pdata['rmsmm'][ind3],color="yellow")
    ax[2].scatter(ldata['fx'][ind4],tdata['rmsmm'][ind4]-pdata['rmsmm'][ind4],color="green")
    ax[2].scatter(ldata['fx'][ind2],tdata['rmsmm'][ind2]-pdata['rmsmm'][ind2],color="orange")
    #ax[2].plot(xx,np.zeros(len(xx)),color="black")
    #ax[2].set_xlim(0.0005,0.005)
    ax[2].set_ylim(-0.001,0.001)
    ax[2].set_ylabel("RMS Spot Motion")
    ax[2].set_title("Moment")
    ax[2].set_xlabel("RMS Spot Motion (PSF Fitting)")

    plt.suptitle("RMS - PSF Fitting Value as a Function of Spot Size (x)")
    plt.savefig("rmsComp4.png")

def makeCompCubes(frameIDs,loadPref,source1,source2,outPref):

    """
    comparison between methods. 
    Fiugre 44
    """
    
    frameId1=frameIDs[0]
    offset=21

    source1="l"
    source2="p"
    source3="w"
    
    #look at changes in parameters with fframe #
    
    ii=0

    #for each frame
    for frameID in frameIDs:
        valsl=np.load(loadPref+"/"+source1+"dump_"+str(frameID)+".npy")
        valsl=ma.masked_invalid(valsl)

        valsp=np.load(loadPref+"/"+source2+"dump_"+str(frameID)+".npy")
        valsp=ma.masked_invalid(valsp)

        valsw=np.load(loadPref+"/"+source3+"dump_"+str(frameID)+".npy")
        valsw=ma.masked_invalid(valsw)
        
        ind=np.where(valsl[:,12] < 130)

        #get values for this frame
        xl=ma.masked_invalid(valsl[ind,1].ravel())
        yl=ma.masked_invalid(valsl[ind,2].ravel())
        fxl=ma.masked_invalid(valsl[ind,6].ravel())
        fyl=ma.masked_invalid(valsl[ind,7].ravel())
        
        xp=ma.masked_invalid(valsp[ind,1].ravel())
        yp=ma.masked_invalid(valsp[ind,2].ravel())
        fxp=ma.masked_invalid(valsp[ind,6].ravel())
        fyp=ma.masked_invalid(valsp[ind,7].ravel())

        xw=ma.masked_invalid(valsw[ind,1].ravel())
        yw=ma.masked_invalid(valsw[ind,2].ravel())
        fxw=ma.masked_invalid(valsw[ind,6].ravel())
        fyw=ma.masked_invalid(valsw[ind,7].ravel())

        
        if(ii==0):
            sz=valsw[ind,:].shape
            cubeX_lw=np.zeros((sz[1],len(frameIDs)))
            cubeY_lw=np.zeros((sz[1],len(frameIDs)))

            cubeXD_lw=np.zeros((sz[1],len(frameIDs)))
            cubeYD_lw=np.zeros((sz[1],len(frameIDs)))

            cubeX_lp=np.zeros((sz[1],len(frameIDs)))
            cubeY_lp=np.zeros((sz[1],len(frameIDs)))

            cubeXD_lp=np.zeros((sz[1],len(frameIDs)))
            cubeYD_lp=np.zeros((sz[1],len(frameIDs)))

            
        #cubeX[:,ii]=xl
        #cubeY[:,ii]=yl
        
        cubeXD_lw[:,ii]=xw-xl
        cubeYD_lw[:,ii]=yw-yl

        cubeXD_lp[:,ii]=xp-xl
        cubeYD_lp[:,ii]=yp-yl

        ii=ii+1
        
    #relative to mean frame
    dXmi,dXma=bounds(cubeXD_lw[:,offset],1.5)
    dYmi,dYma=bounds(cubeYD_lw[:,offset],1.5)

    fig,ax=plt.subplots(2,2,figsize=(14,10))
        
    scFirst=ax[0,0].scatter(xl,yl,c=cubeXD_lw[:,offset])
    ax[0,0].set_title("dx (Windowed-MCS)")
    fig.colorbar(scFirst,ax=ax[0,0])
    scFirst=ax[0,1].scatter(xl,yl,c=cubeYD_lw[:,offset])
    ax[0,1].set_title("dy (Windowed-MCS)")
    fig.colorbar(scFirst,ax=ax[0,1])

    dXmi,dXma=bounds(cubeXD_lp[:,offset],1.5)
    dYmi,dYma=bounds(cubeYD_lp[:,offset],1.5)

    scFirst=ax[1,0].scatter(xl,yl,c=cubeXD_lp[:,offset])
    ax[1,0].set_title("dx (PSF Fitting-MCS)")
    fig.colorbar(scFirst,ax=ax[1,0])
    scFirst=ax[1,1].scatter(xl,yl,c=cubeYD_lp[:,offset])
    ax[1,1].set_title("dy (PSF Fitting-MCS)")
    fig.colorbar(scFirst,ax=ax[1,1])

    plt.suptitle("Comparison of Methods ["+str(frameIDs[offset])+"]")          
    plt.savefig("subComp.png")
        
import matplotlib.pylab as plt
import numpy as np

"""
plot spot to spot variation in rms.
Fiugre 31
"""

ff=np.loadtxt("rmsmean.txt")

frameID=ff[:,0]
t=ff[:,1]
za=ff[:,2]
inr=ff[:,3]
rmsXs=ff[:,4]
rmsYs=ff[:,5]
rmss=ff[:,6]
rmsX=ff[:,7]
rmsY=ff[:,8]
rms=ff[:,0]

fig,ax=plt.subplots()

m1=np.zeros(16)
m2=[1,1.1,1.2,1.3,2,2.1,2.2,2.3,3,3.1,3.2,3.3,4,4.1,4.2,4.3]
m3=np.zeros(16)

v1=[]
v2=[]
v3=[]
v4=[]
t1=[]
t2=[]
t3=[]
t4=[]


for i in range(len(t)):
    if(t[i]==0.5):
        tt=1
    elif(t[i]==1):
        tt=2
    elif(t[i]==2):
        tt=3
    elif(t[i]==5):
        tt=4
    else:
        tt=-20

    if(tt > 0):

        #// 
        if((180 % inr[i]) != 0):
            v1.append(rmsX[i])
            t1.append(tt)
            v2.append(rmsY[i])
            t2.append(tt+0.1)
            v3.append(rmsXs[i])
            t3.append(tt+0.2)
            v4.append(rmsYs[i])
            t4.append(tt+0.3)
            
            m1[4*(tt-1)]+=rmsX[i]
            m1[4*(tt-1)+1]+=rmsY[i]
            m1[4*(tt-1)+2]+=rmsXs[i]
            m1[4*(tt-1)+3]+=rmsYs[i]
            m3[4*(tt-1)]+=1
            m3[4*(tt-1)+1]+=1
            m3[4*(tt-1)+2]+=1
            m3[4*(tt-1)+3]+=1


ax.scatter(t1,v1,color="blue",label="RMS X")
ax.scatter(t2,v2,color="green",label="RMS Y")
ax.scatter(t3,v3,color="orange",label="RMS X - smth")
ax.scatter(t4,v4,color="red",label="RMS Y - smth")

plt.ylim([0,0.003])
plt.xlim([0.5,5])

ax.scatter(m2,m1/m3,color="black")

ax.legend()
plt.xticks([1.15,2.15,3.15,4.15],['0.5','1','2','5'])
ax.set_title("Spot to Spot Variation in RMS (mm)")
plt.savefig("rms_t1.png")

fig,ax=plt.subplots()
ax.scatter(za,rmsXs)
ax.scatter(za+2,rmsY)
ax.scatter(za+4,rmsXs)
ax.scatter(za+6,rmsYs)

plt.savefig("rms_za.png")

fig,ax=plt.subplots()
ax.scatter(inr,rmsXs)
ax.scatter(inr+10,rmsY)
ax.scatter(inr+20,rmsXs)
ax.scatter(inr+30,rmsYs)

plt.savefig("rms_inr.png")

#rotAvParm()
#rotDeltaS()
#byFrameFX()
#byFrameFY()
#byFrameP()
#fits()
#elevation()
#byFrameDX()
#byFrameDY()
#rmsRot()
#rmsTime()
#methods()
#exampleRMS()
#exampleRMS1()
#yDir()
#modelMoire()
#smoothCpt()
#smoothCpt1()
#positions()
#exampleRMS2()

#simStats(1762200)
#simStats(1762300)
#simStats(1762400)
#simStats(1762500)
#simVals()
#offsetsX("w","p","Mean dx win-PSF")
#offsetsX("w","l","Mean dx win-MCS")
#offsetsX("l","p","Mean dx MCS-PFS")
#offsetsY("w","p","Mean dy win-PSF")
#offsetsY("w","l","Mean dy win-MCS")
#offsetsY("l","p","Mean dy MCS-PFS")

#plotSimStats()
#offsets()
#simVals()
#smoothCpt1()
#diagPlot()
#fits()
#elevation()
frames=np.array([1762200,1762300,1762400,1762500]).astype('int')


#for frameID in frames:
#    methodsRMS(frameID)

#offsets()
