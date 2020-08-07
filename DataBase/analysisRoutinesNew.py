

import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import sigmaclip
import sys
from astropy.io import fits
import numpy.ma as ma
import sexRoutines as sex

import os

import scipy
import analysisTidy as at

import meetingPlotsTidy as mp



#the try/except is for using in the MHS environment, or standalone. 

try:
    import mcsActor.Visualization.mcsRoutines as mcs
except:
    import mcsRoutines as mcs

try:
    import mcsActor.Visualization.fpsRoutines as fps
except:
    import fpsRoutines as fps

try:
    import mcsActor.Visualization.visRoutines as vis
except:
    import visRoutines as vis

try:
    import mcsActor.Visualization.plotRoutines as visplot
except:
    import plotRoutines as visplot


def differences(frameIDs,loadPref,st,source1,source2,dataType,suffix,ff):

    frameID=frameIDs[0]
    cubeX1,cubeY1,cubeFX1,cubeFY1,cubeP1,cubeFXs1,cubeFYs1,cubePs1,cubeXD1,cubeYD1,cubeXDs1,cubeYDs1,xAv1,yAv1,fxAv1,fxAv1,fyAv1,pAv1,fxAvs1,fyAvs1,pAvs1,rmsVal1,rmsX1,rmsY1,rmsVals1,rmsXs1,rmsYs1=loadCubes(loadPref,source1,frameID)


    cubeX2,cubeY2,cubeFX2,cubeFY2,cubeP2,cubeFXs2,cubeFYs2,cubePs2,cubeXD2,cubeYD2,cubeXDs2,cubeYDs2,xAv2,yAv2,fxAv2,fxAv2,fyAv2,pAv2,fxAvs2,fyAvs2,pAvs2,rmsVal2,rmsX2,rmsY2,rmsVals2,rmsXs2,rmsYs2=loadCubes(loadPref,source2,frameID)

    print("%8d %s %6.3f %6.3f %6.3f %6.3f %6.3f %6.3f" % (frameIDs[0],st,1000*(xAv1-xAv2).mean(),1000*(yAv1-yAv2).mean(),1000*rmsVal1.mean(),fxAv1.mean(),fyAv1.mean(),(fxAv1/fyAv1).mean()),file=ff)

    
def rewrite(frameIDs,filesUC,files,bits='int16'):

    """
    simple wrapper to read and re-write MCS images for accessibility by SeXtractor, which
    doesn't handle compression. 

    set bits to 'int32' for exposures > 5s.

    """
  
    for f1,f2 in zip(files,filesUC):
        print(f1,f2)
        data,hdr=fits.getdata(f1,1,header=True)
        fits.writeto(f2,data.astype(bits),hdr,clobber=True)

def loadCubes(loadPref,source,frameID):

    """

    Utility routine to read in pre-computed analysis cubes; reads files created by 
    makeCubesAll.

    input: loadPref: prefix for loading (ie, directory path)
           source: key for reduction method (p, w, etc)
           frameID: frame# (including moveID if appropriate)



    Variable names are of the form

    cube: cubes of values
          X, Y - position
          FX, FY - size
          P - peak value
          XD, XD - different in position from mean value

    Averages
          ?Av -
           x,y - position
           fx,fy - size
           p - peak

    RMS
          Val - total rms
          X, y - in a asingle driection

    adding an s to the above means smooth component subtracted


    """

    #cubeD contains 3D cubes, cubeA contains 2D averages
    
    dd=np.load(loadPref+"/"+source+"cubeD_"+str(frameID)+".npy")
    aa=np.load(loadPref+"/"+source+"cubeA_"+str(frameID)+".npy")

    #extract the values
    
    cubeX=dd[0].astype('float')
    cubeY=dd[1].astype('float')
    cubeFX=dd[2].astype('float')
    cubeFY=dd[3].astype('float')
    cubeP=dd[4].astype('float')
    cubeFXs=dd[5].astype('float')
    cubeFYs=dd[6].astype('float')
    cubePs=dd[7].astype('float')
    cubeXD=dd[8].astype('float')
    cubeYD=dd[9].astype('float')
    cubeXDs=dd[10].astype('float')
    cubeYDs=dd[11].astype('float')

    xAv=aa[0].astype('float')
    yAv=aa[1].astype('float')
    fxAv=aa[2].astype('float')
    fyAv=aa[4].astype('float')
    pAv=aa[5].astype('float')
    fxAvs=aa[6].astype('float')
    fyAvs=aa[7].astype('float')
    pAvs=aa[8].astype('float')
    rmsVal=aa[9].astype('float')
    rmsX=aa[10].astype('float')
    rmsY=aa[11].astype('float')
    rmsVals=aa[12].astype('float')
    rmsXs=aa[13].astype('float')
    rmsYs=aa[14].astype('float')

    #mask invalid values

    cubeX=ma.masked_invalid(cubeX)
    cubeY=ma.masked_invalid(cubeY)
    cubeP=ma.masked_invalid(cubeP)
    cubeFX=ma.masked_invalid(cubeFX)
    cubeFY=ma.masked_invalid(cubeFY)
    cubePs=ma.masked_invalid(cubePs)
    cubeFXs=ma.masked_invalid(cubeFXs)
    cubeFYs=ma.masked_invalid(cubeFYs)
    cubeXD=ma.masked_invalid(cubeXD)
    cubeYD=ma.masked_invalid(cubeYD)
    cubeXDs=ma.masked_invalid(cubeXDs)
    cubeYDs=ma.masked_invalid(cubeYDs)

    xAv=ma.masked_invalid(xAv)
    yAv=ma.masked_invalid(yAv)
    pAv=ma.masked_invalid(pAv)
    fxAv=ma.masked_invalid(fxAv)
    fyAv=ma.masked_invalid(fyAv)
    pAvs=ma.masked_invalid(pAvs)
    fxAvs=ma.masked_invalid(fxAvs)
    fyAvs=ma.masked_invalid(fyAvs)
    rmsVal=ma.masked_invalid(rmsVal)
    rmsX=ma.masked_invalid(rmsX)
    rmsY=ma.masked_invalid(rmsY)
    rmsVals=ma.masked_invalid(rmsVals)
    rmsXs=ma.masked_invalid(rmsXs)
    rmsYs=ma.masked_invalid(rmsYs)
    
    #and return
    
    return cubeX,cubeY,cubeFX,cubeFY,cubeP,cubeFXs,cubeFYs,cubePs,cubeXD,cubeYD,cubeXDs,cubeYDs,xAv,yAv,fxAv,fxAv,fyAv,pAv,fxAvs,fyAvs,pAvs,rmsVal,rmsX,rmsY,rmsVals,rmsXs,rmsYs
 
def calcSmooth(A,x,y,F):

    Fmod = A[0]+ x*A[1] + y*A[2] + x**2*A[3] + x**2*y*A[4] + x**2*y**2*A[5] + y**2*A[6] + x*y**2*A[7] + x*y*A[8]

    diff = Fmod - F

    return diff.flatten()

def subSmoothMP(xAv,yAv,parm):


    X=xAv[~xAv.mask].flatten()
    Y=yAv[~xAv.mask].flatten()
    F=parm[~xAv.mask].flatten()

    X1=xAv.flatten()
    Y1=yAv.flatten()

    A1=np.array([X1*0+1, X1, Y1, X1**2, X1**2*Y1, X1**2*Y1**2, Y1**2, X1*Y1**2, X1*Y1]).T

    result = scipy.optimize.leastsq(calcSmooth,np.ones((9)),args=(X,Y,F))
    smooth=np.zeros((len(xAv)))

    result=np.array(result[0])
    for i in range(len(result)):
        smooth=smooth+A1.T[i]*result[i]

    parm=parm-smooth
    return parm,smooth
    
def subSmooth(xAv,yAv,parm,order):

    #fit a smooth surface to a set of x,y points and values

    #first order fits linear in x and y,
    #second order fits quadratic in X and Y

    #values

    X=xAv[~xAv.mask].flatten()
    Y=yAv[~xAv.mask].flatten()
    F=parm[~xAv.mask].flatten()

    X1=xAv.flatten()
    Y1=yAv.flatten()

    parmMean=F.mean()

    
    #array of X Axes for fit

    if(order==1):
        A=np.array([ X*0+1,  X,  Y]).T
        A1=np.array([X1*0+1, X1, Y1]).T
    if(order==2):
        A=np.array([ X*0+1,  X,  Y,  X**2,  X**2*Y,   X**2* Y**2,  Y**2,  X*Y**2,   X* Y]).T
        A1=np.array([X1*0+1, X1, Y1, X1**2, X1**2*Y1, X1**2*Y1**2, Y1**2, X1*Y1**2, X1*Y1]).T

    #do fit
    coeff,r,rank,s=np.linalg.lstsq(A,F,rcond=-1)
    smooth=np.zeros((len(xAv)))
    #create smooth map
    for i in range(len(coeff)):
        smooth=smooth+coeff[i]*A1.T[i]

    parmSub=parm-smooth


    return parmSub,smooth



def bounds(val,rms):


    vmi=val.mean()-rms*val.std()
    vma=val.mean()+rms*val.std()

    return vmi,vma


def simStats(frameIDs,loadPref,source):

    cubeX,cubeY,cubeFX,cubeFY,cubeP,cubeFXs,cubeFYs,cubePs,cubeXD,cubeYD,cubeXDs,cubeYDs,xAv,yAv,fxAv,fxAv,fyAv,pAv,fxAvs,fyAvs,pAvs,rmsVal,rmsX,rmsY,rmsVals,rmsXs,rmsYs=loadCubes(loadPref,source,frameIDs[0])

    print("%2s %8d %5.2f %5.3f %5.2f" %(source,frameIDs[0],1000*abs(cubeXD).mean(),1000*abs(cubeYD).mean(),1000*rmsVal.mean()))

    
def rmsVals(frameIDs,st,loadPref,outPref,source,vmi,vma,ofile,day,doPlot):

    frameId1=frameIDs[0]
    #look at changes in parameters with fframe #
    130
    ii=0

    #load teh data nad put it in an array
    ii=0

    phys=at.readPhysTidy(frameIDs[0],"NewSet")
    for frameID in frameIDs:

        vals=np.load(loadPref+"/"+source+"dump_"+str(int(frameID))+".npy")
        ind=np.where(vals[:,12] < 130)

        if(ii==0):
            xfirst=np.array(vals[ind,1].ravel())
            yfirst=np.array(vals[ind,2].ravel())
            xArray=np.empty((len(xfirst),len(frameIDs)))
            yArray=np.empty((len(xfirst),len(frameIDs)))
            fxArray=np.empty((len(xfirst),len(frameIDs)))
            fyArray=np.empty((len(xfirst),len(frameIDs)))
            peakArray=np.empty((len(xfirst),len(frameIDs)))

        xArray[:,ii]=vals[ind,1].ravel()
        yArray[:,ii]=vals[ind,2].ravel()
        fxArray[:,ii]=vals[ind,6].ravel()
        fyArray[:,ii]=vals[ind,7].ravel()
        peakArray[:,ii]=vals[ind,8].ravel()

        ii+=1

    #mask nans
    xfirst=ma.masked_invalid(xfirst)
    yfirst=ma.masked_invalid(yfirst)
    xArray=ma.masked_invalid(xArray)
    yArray=ma.masked_invalid(yArray)
    fxArray=ma.masked_invalid(fxArray)
    fyArray=ma.masked_invalid(fyArray)
    peakArray=ma.masked_invalid(peakArray)
    backArray=ma.masked_invalid(peakArray)
               
    #get the rms
    
    dd=np.zeros(xArray.shape)
    xd=np.zeros(xArray.shape)
    yd=np.zeros(xArray.shape)
    for i in range(xArray.shape[0]):
        dd[i,:]=np.sqrt((xArray[i,:]-xfirst[i])**2+(yArray[i,:]-yfirst[i])**2)
        xd[i,:]=np.sqrt((xArray[i,:]-xfirst[i])**2)
        yd[i,:]=np.sqrt((yArray[i,:]-yfirst[i])**2)
 
    dd=ma.masked_where(((dd <=0) | (xArray.mask == True)), dd)
    xd=ma.masked_where(((dd <=0) | (xArray.mask == True)), xd)
    yd=ma.masked_where(((dd <=0) | (xArray.mask == True)), yd)
    xAv=xArray.mean(axis=1)
    yAv=yArray.mean(axis=1)

    rmsVal=dd.std(axis=1)
    rmsX=xd.std(axis=1)
    rmsY=yd.std(axis=1)
    fxAv=fxArray.mean(axis=1)
    fyAv=fyArray.mean(axis=1)

    #remove smoothed component

    szR,smR=subSmoothMP(xAv,yAv,rmsVal)
    szRx,smRx=subSmoothMP(xAv,yAv,rmsX)
    szRy,smRy=subSmoothMP(xAv,yAv,rmsY)
    szFX,smFX=subSmoothMP(xAv,yAv,fxAv)
    
    #do the plots
    plotRange=None
    nbins=30

    prefix=str(frameIDs[0])
    inter=0
    stitle=""
    #print(rmsVal.min(),rmsVal.mean(),rmsVal.max())

    print(frameIDs[0],frameIDs[-1],phys['t'],phys['za'],phys['inr'],day,phys['date'],phys['hst'],rmsVal.mean(),rmsX.mean(),rmsY.mean(),fxArray.mean(),fyArray.mean(),mp.calcReltime(phys['hst']),file=ofile)

    if(doPlot==1):
        stitle="  ["+st+"] (-smooth)"
        plotRange=[0.002,0.006]
        plotRange=None
        visplot.pairPlot(xfirst,yfirst,rmsVal,rmsVal.ravel(),plotRange,"RMS",prefix,"_rms"+source,"RMS","mm",nbins,inter,stitle=stitle)

        rXmi,rXma=bounds(szRx,1.5)
        rYmi,rYma=bounds(szRy,1.5)
        rmi,rma=bounds(szR,1.5)
        FXmi,FXma=bounds(szFX,1.5)
    
        fig,ax=plt.subplots(2,2,figsize=(14,10))
        scFirst=ax[0,0].scatter(xAv,yAv,c=szRx,vmin=rXmi,vmax=rXma)
        fig.colorbar(scFirst,ax=ax[0,0])
        ax[0,0].set_title("RMS X (sm sub)")
        scFirst=ax[0,1].scatter(xAv,yAv,c=szRy,vmin=rYmi,vmax=rYma)
        fig.colorbar(scFirst,ax=ax[0,1])
        ax[0,1].set_title("RMS Y (sm sub)")
        scFirst=ax[1,0].scatter(xAv,yAv,c=szR,vmin=rmi,vmax=rma)
        fig.colorbar(scFirst,ax=ax[1,0])
        ax[1,0].set_title("RMS (sm sub)")
        scFirst=ax[1,1].scatter(xAv,yAv,c=szFX,vmin=FXmi,vmax=FXma)
        fig.colorbar(scFirst,ax=ax[1,1])
        ax[1,1].set_title("FWHM X (sm sub)")

        plt.suptitle("RMS - "+str(frameIDs[0])+" "+st+" (mm)")
        plt.savefig(outPref+"/"+source+"rms_"+str(frameIDs[0])+"_sm.png")
        
        rXmis,rXmas=bounds(rmsX,1.5)
        rYmis,rYmas=bounds(rmsY,1.5)
        rmis,rmas=bounds(rmsVal,1.5)
        #FXmi,FXma=bounds(szFX,1.5)
    
        fig,ax=plt.subplots(2,2,figsize=(14,10))
        scFirst=ax[0,0].scatter(xAv,yAv,c=rmsX,vmin=rXmis,vmax=rXmas)
        fig.colorbar(scFirst,ax=ax[0,0])
        ax[0,0].set_title("RMS X")
        scFirst=ax[0,1].scatter(xAv,yAv,c=rmsY,vmin=rYmis,vmax=rYmas)
        fig.colorbar(scFirst,ax=ax[0,1])
        ax[0,1].set_title("RMS Y")
        scFirst=ax[1,0].scatter(xAv,yAv,c=rmsVal,vmin=rmis,vmax=rmas)
        fig.colorbar(scFirst,ax=ax[1,0])
        ax[1,0].set_title("RMS")
        scFirst=ax[1,1].scatter(xAv,yAv,c=szFX,vmin=FXmi,vmax=FXma)
        fig.colorbar(scFirst,ax=ax[1,1])
        ax[1,1].set_title("FWHM X")

        plt.suptitle("RMS - "+str(frameIDs[0])+" "+st+" (mm)")
        plt.savefig(outPref+"/"+source+"rms_"+str(frameIDs[0])+".png")

def movies(frameIDs,st,loadPref,movPref,source,vma,delFiles,smth):

    plt.close('all')
    frameId1=frameIDs[0]
    #look at changes in parameters with fframe #

    ii=0

    prefix=source+"_"
    if(vma[0]==0):
        prefix=prefix+"scaled_"
    else:
        prefix=prefix+"fixed_"
    if(smth==1):
        prefix=prefix+"smooth"
    else:
        prefix=prefix+"raw"

    cubeX,cubeY,cubeFX,cubeFY,cubeP,cubeFXs,cubeFYs,cubePs,cubeXD,cubeYD,cubeXDs,cubeYDs,xAv,yAv,fxAv,fxAv,fyAv,pAv,fxAvs,fyAvs,pAvs,rmsVal,rmsX,rmsY,rmsVals,rmsXs,rmsYs=loadCubes(loadPref,source,frameIDs[0])

    #for each frame
    for frameID in frameIDs:

        print(frameID,cubeYDs.mean())
        fig,ax=plt.subplots(2,2,figsize=(14,10))


        FXmi,FXma=bounds(cubeFXs[:,ii],1.5)
        FYmi,FYma=bounds(cubeFYs[:,ii],1.5)
        
        
        if(smth==0):
            

            if(vma[0]==0):
                #relative to mean frame
                dXmi,dXma=bounds(cubeXD[:,ii],1.5)
                dYmi,dYma=bounds(cubeYD[:,ii],1.5)
                
                scFirst=ax[0,0].scatter(xAv,yAv,c=cubeXD[:,ii])
                ax[0,0].set_title("dx")
                fig.colorbar(scFirst,ax=ax[0,0])
                scFirst=ax[0,1].scatter(xAv,yAv,c=cubeYD[:,ii])
                ax[0,1].set_title("dy")
                fig.colorbar(scFirst,ax=ax[0,1])
            else:
                #relative to mean frame
                scFirst=ax[0,0].scatter(xAv,yAv,c=cubeXD[:,ii],vmin=-vma[0],vmax=vma[0])
                ax[0,0].set_title("dx")
                fig.colorbar(scFirst,ax=ax[0,0])
                scFirst=ax[0,1].scatter(xAv,yAv,c=cubeYD[:,ii],vmin=-vma[0],vmax=vma[0])
                ax[0,1].set_title("dy")
                fig.colorbar(scFirst,ax=ax[0,1])
        else:

            dXmi,dXma=bounds(cubeXDs[:,ii],1.5)
            dYmi,dYma=bounds(cubeYDs[:,ii],1.5)

            if(vma[0]==0):
                #relative to mean frame
                scFirst=ax[0,0].scatter(xAv,yAv,c=cubeXDs[:,ii],vmin=dXmi,vmax=dXma)
                ax[0,0].set_title("dx")
                fig.colorbar(scFirst,ax=ax[0,0])
                scFirst=ax[0,1].scatter(xAv,yAv,c=cubeYDs[:,ii],vmin=dYmi,vmax=dYma)
                ax[0,1].set_title("dy")
                fig.colorbar(scFirst,ax=ax[0,1])
            else:
                #relative to mean frame
                scFirst=ax[0,0].scatter(xAv,yAv,c=cubeXDs[:,ii],vmin=-vma[0],vmax=vma[0])
                ax[0,0].set_title("dx")
                fig.colorbar(scFirst,ax=ax[0,0])
                scFirst=ax[0,1].scatter(xAv,yAv,c=cubeYDs[:,ii],vmin=-vma[0],vmax=vma[0])
                ax[0,1].set_title("dy")
                fig.colorbar(scFirst,ax=ax[0,1])

        if(len(vma)==1):
            scFWHM=ax[1,0].scatter(xAv,yAv,c=cubeFXs[:,ii])
            fig.colorbar(scFWHM,ax=ax[1,0])
            ax[1,0].set_title("FWHM (x)")
            scFWHM=ax[1,1].scatter(xAv,yAv,c=cubeFYs[:,ii])
            fig.colorbar(scFWHM,ax=ax[1,1])
            ax[1,1].set_title("FWHM (y)")
        else:
            scFWHM=ax[1,0].scatter(xAv,yAv,c=cubeFXs[:,ii],vmax=vma[1])
            fig.colorbar(scFWHM,ax=ax[1,0])
            ax[1,0].set_title("FWHM (x)")
            scFWHM=ax[1,1].scatter(xAv,yAv,c=cubeFYs[:,ii],vmax=vma[2])
            fig.colorbar(scFWHM,ax=ax[1,1])
            ax[1,1].set_title("FWHM (y)")
 

        if(source=="l"):
            plt.suptitle("MCS Algorithm - "+str(frameIDs[0])+" "+st+" (mm)")
        if(source=="m"):
            plt.suptitle("Minus Mean: Gaussian smoothed moment - "+str(frameIDs[0])+" "+st+" (mm)")
        if(source=="t"):
            plt.suptitle("Minus Mean: top hat smoothed moment - "+str(frameIDs[0])+" "+st+" (mm)")
        if(source=="u"):
            plt.suptitle("Minus Mean: Unsmoothed smoothed moment - "+str(frameIDs[0])+" "+st+" (mm)")
        if(source=="p"):
            plt.suptitle("Minus Mean: PSFEx fitting - "+str(frameIDs[0])+" "+st+" (mm)")
        if(source=="w"):
            plt.suptitle("Minus Mean: windowed moments - "+str(frameIDs[0])+" "+st+" (mm)")
        if(source=="s"):
            plt.suptitle("Minus Mean: windowed moments (sims) - "+str(frameIDs[0])+" "+st+" (mm)")
        if(source=="n"):
            plt.suptitle("Minus Mean: standard - "+str(frameIDs[0])+" "+st+" (mm)")
        
        plt.savefig(movPref+"/"+source+"_"+prefix+"_"+str(frameId1)+"_"+str(int(ii)).zfill(2)+".png")
        
        ii=ii+1
            
    os.system("ffmpeg -framerate 3 -i "+movPref+"/"+source+"_"+prefix+"_"+str(frameId1)+"_%02d.png "+source+"_"+prefix+"_"+str(frameId1)+".mp4")

    #if(delFiles==1):
    #    os.system("rm "+outPref+"_*.png")
def moviesRedo(frameIDs,st,loadPref,movPref,source,vma,delFiles,smth):


    frameId1=frameIDs[0]

    prefix=source+"_"
    if(smth==1):
        prefix=prefix+"smooth_"
    else:
        prefix=prefix+"raw_"
    if(vma[0]==0):
        prefix=prefix+"scaled"
    else:
        prefix=prefix+"fixed"
    os.system("ffmpeg -framerate 3 -i "+movPref+"/"+prefix+"_"+str(frameId1)+"_%02d.png "+prefix+"_"+str(frameId1)+".mp4")

def makeCompCubes(frameIDs,loadPref,source1,source2,outPref):

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
        
        
def makeCubesAll(frameIDs,loadPref,source):

    #cubeX - x differences
    #cubeYs - x differences, smoothed component subtracted
    #cubeFX - fwhms
    #cube FXx - fwhms, smoothed subtracted 

    frameId1=frameIDs[0]
    #look at changes in parameters with fframe #
    
    ii=0

    #for each frame
    for frameID in frameIDs:
        vals=np.load(loadPref+"/"+source+"dump_"+str(frameID)+".npy")
        vals=ma.masked_invalid(vals)

        #trim by position
        ind=np.where(vals[:,12] < 130)

        #get values for this frame
        x=ma.masked_invalid(vals[ind,1].ravel())
        y=ma.masked_invalid(vals[ind,2].ravel())
        fx=ma.masked_invalid(vals[ind,6].ravel())
        fy=ma.masked_invalid(vals[ind,7].ravel())
        peak=ma.masked_invalid(vals[ind,8].ravel())

        #set up arrays on first pass
        if(ii==0):
            xfirst=np.array(vals[ind,1].ravel())
            yfirst=np.array(vals[ind,2].ravel())
            
            sz=vals[ind,:].shape
            cubeX=np.zeros((sz[1],len(frameIDs)))
            cubeY=np.zeros((sz[1],len(frameIDs)))

            cubeFX=np.zeros((sz[1],len(frameIDs)))
            cubeFY=np.zeros((sz[1],len(frameIDs)))
            cubeFXs=np.zeros((sz[1],len(frameIDs)))
            cubeFYs=np.zeros((sz[1],len(frameIDs)))

            cubeXD=np.zeros((sz[1],len(frameIDs)))
            cubeYD=np.zeros((sz[1],len(frameIDs)))
            cubeXDs=np.zeros((sz[1],len(frameIDs)))
            cubeYDs=np.zeros((sz[1],len(frameIDs)))

            cubeP=np.zeros((sz[1],len(frameIDs)))
            cubePs=np.zeros((sz[1],len(frameIDs)))

        #smoothed versions of fwhm
        avfx=ma.masked_invalid(fx).mean()
        avfy=ma.masked_invalid(fy).mean()
        avp=ma.masked_invalid(peak).mean()
        szX,smX=subSmoothMP(x,y,fx)
        szY,smY=subSmoothMP(x,y,fy)
        szP,smP=subSmoothMP(x,y,peak)
        szX=szX+avfx
        szY=szY+avfy
        szP=szP+avp
        
        #base values
        cubeX[:,ii]=x
        cubeY[:,ii]=y
        cubeP[:,ii]=peak
        cubePs[:,ii]=szP
        
        cubeFX[:,ii]=ma.masked_invalid(fx)
        cubeFY[:,ii]=ma.masked_invalid(fy)
        cubeFXs[:,ii]=ma.masked_invalid(szX)
        cubeFYs[:,ii]=ma.masked_invalid(szY)
        ii=ii+1

    #mean values of positions
    xM=ma.masked_invalid(cubeX).mean(axis=1)
    yM=ma.masked_invalid(cubeY).mean(axis=1)
    fxM=ma.masked_invalid(cubeFX).mean(axis=1)
    fyM=ma.masked_invalid(cubeFY).mean(axis=1)
    fxMs=ma.masked_invalid(cubeFXs).mean(axis=1)
    fyMs=ma.masked_invalid(cubeFYs).mean(axis=1)
    pM=ma.masked_invalid(cubeP).mean(axis=1)
    pMs=ma.masked_invalid(cubePs).mean(axis=1)

    cubeX=ma.masked_invalid(cubeX)

    #get the rms
    dd=np.zeros(cubeX.shape)
    xd=np.zeros(cubeX.shape)
    yd=np.zeros(cubeX.shape)
    for i in range(cubeX.shape[0]):
        dd[i,:]=np.sqrt((cubeX[i,:]-xfirst[i])**2+(cubeY[i,:]-yfirst[i])**2)
        xd[i,:]=np.sqrt((cubeX[i,:]-xfirst[i])**2)
        yd[i,:]=np.sqrt((cubeY[i,:]-yfirst[i])**2)

    dd=ma.masked_where(((dd <=0) | (cubeX.mask == True)), dd)
    xd=ma.masked_where(((dd <=0) | (cubeX.mask == True)), xd)
    yd=ma.masked_where(((dd <=0) | (cubeX.mask == True)), yd)
    
    rmsVal=dd.std(axis=1)
    rmsX=xd.std(axis=1)
    rmsY=yd.std(axis=1)

    rVa=rmsVal.mean()
    rXa=rmsX.mean()
    rYa=rmsY.mean()

    szR,smR=subSmoothMP(xM,yM,rmsVal)
    szRx,smRx=subSmoothMP(xM,yM,rmsX)
    szRy,smRy=subSmoothMP(xM,yM,rmsY)

    szR=szR+rVa
    szRx=szRx+rXa
    szRy=szRx+rYa
    
    #now a second pass to calculate differences in position to mean
    ii=0
    for frameID in frameIDs:
        print(loadPref+"/"+source+"dump_"+str(frameID)+".npy")
        vals=0
        vals=np.load(loadPref+"/"+source+"dump_"+str(frameID)+".npy")
        vals=ma.masked_invalid(vals)

        ind=np.where(vals[:,12] < 130)
        
        x2=ma.masked_invalid(vals[ind,1].ravel())
        y2=ma.masked_invalid(vals[ind,2].ravel())

        xdMean=x2-xM
        ydMean=y2-yM

        ddMean=np.sqrt(xdMean**2+ydMean**2)

        cubeXD[:,ii]=xdMean
        cubeYD[:,ii]=ydMean

        avXD=xdMean.mean()
        avYD=ydMean.mean()
            
        #an estimate of the FWHM size
        #fwhm=np.sqrt(vals[ind,6]**2+vals[ind,7]**2)
            
        szX,smX=subSmoothMP(xM,yM,xdMean)
        szY,smY=subSmoothMP(xM,yM,ydMean)

        cubeXDs[:,ii]=szX+avXD
        cubeYDs[:,ii]=szY+avYD

        ii=ii+1

    
    #now make the dump arrays
    
    dumpArrayD=[cubeX,cubeY,cubeFX,cubeFY,cubeP,cubeFXs,cubeFYs,cubePs,cubeXD,cubeYD,cubeXDs,cubeYDs]
    dumpArrayA=[xM,yM,fxM,fxM,fyM,pM,fxMs,fyMs,pMs,rmsVal,rmsX,rmsY,szR,szRx,szRy]
    
    np.save(loadPref+"/"+source+"cubeD_"+str(frameIDs[0])+".npy",dumpArrayD,allow_pickle=False)
    np.save(loadPref+"/"+source+"cubeA_"+str(frameIDs[0])+".npy",dumpArrayA,allow_pickle=False)
    
def patterns(frameIDs,st,loadPref,source,vma,delFiles):

    dumpArray=np.load(loadPref+"/"+source+"cube1_"+str(frameIDs[0])+".npy")

    cubeX=dumpArray[0].astype('float')
    cubeY=dumpArray[1].astype('float')
    cubeFXs=dumpArray[2].astype('float')
    cubeFYs=dumpArray[3].astype('float')

    cubeX=ma.masked_invalid(cubeX)
    cubeF=ma.masked_invalid(cubeY)
    cubeFXs=ma.masked_invalid(cubeFXs)
    cubeFYs=ma.masked_invalid(cubeFYs)

    dumpArray=np.load(loadPref+"/"+source+"cube2_"+str(frameIDs[0])+".npy")

    cubeXD=dumpArray[0].astype('float')
    cubeYD=dumpArray[1].astype('float')
    cubeXDs=dumpArray[2].astype('float')
    cubeYDs=dumpArray[3].astype('float')

    cubeXD=ma.masked_invalid(cubeXD)
    cubeYD=ma.masked_invalid(cubeYD)
    cubeXDs=ma.masked_invalid(cubeXDs)
    cubeYDs=ma.masked_invalid(cubeYDs)

    dumpArray=np.load(loadPref+"/"+source+"cube3_"+str(frameIDs[0])+".npy")

    x1=dumpArray[0].astype('float')
    y1=dumpArray[1].astype('float')

    fx=cubeFXs.mean(axis=1)
    fy=cubeFYs.mean(axis=1)
        
    xdMean=cubeXD.mean(axis=1)
    ydMean=cubeYD.mean(axis=1)
    xdSMean=np.abs(cubeXDs).mean(axis=1)
    ydSMean=np.abs(cubeYDs).mean(axis=1)
    xdMed=np.median(cubeXD,axis=1)
    ydMed=np.median(cubeYD,axis=1)
    xdSMed=np.median(cubeXDs,axis=1)
    ydSMed=np.median(cubeYDs,axis=1)

    xmi=fx.mean()-1.5*cubeFXs.std()
    xma=fx.mean()+1.5*cubeFXs.std()
    ymi=fy.mean()-1.5*cubeFYs.std()
    yma=fy.mean()+1.5*cubeFYs.std()

    
    #relative to first frame
    fig,ax=plt.subplots(2,2,figsize=(14,10))
    scFirst=ax[0,0].scatter(x1,y1,c=xdSMean)
    ax[0,0].set_title("xdSMed")
    fig.colorbar(scFirst,ax=ax[0,0])
    scFirst=ax[0,1].scatter(x1,y1,c=ydSMean)
    ax[0,1].set_title("ydSMed")
    fig.colorbar(scFirst,ax=ax[0,1])
    scFWHM=ax[1,0].scatter(x1,y1,c=fx,vmin=xmi,vmax=xma)
    fig.colorbar(scFWHM,ax=ax[1,0])
    ax[1,0].set_title("FWHM (x)")
    scFWHM=ax[1,1].scatter(x1,y1,c=fy,vmin=ymi,vmax=yma)
    fig.colorbar(scFWHM,ax=ax[1,1])
    ax[1,1].set_title("FWHM (y)")

    if(source=="m"):
        plt.suptitle("Averages: Gaussian smoothed moment - "+str(frameIDs[0])+" "+st+" (mm)")
    if(source=="t"):
        plt.suptitle("Averages: top hat smoothed moment - "+str(frameIDs[0])+" "+st+" (mm)")
    if(source=="u"):
        plt.suptitle("Averages: Unsmoothed smoothed moment - "+str(frameIDs[0])+" "+st+" (mm)")
    if(source=="p"):
        plt.suptitle("Averages: PSFEx fitting - "+str(frameIDs[0])+" "+st+" (mm)")
    if(source=="w"):
        plt.suptitle("Averages: windowed moments - "+str(frameIDs[0])+" "+st+" (mm)")
    if(source=="w"):
        plt.suptitle("Averages: windowed moments (sims) - "+str(frameIDs[0])+" "+st+" (mm)")
    if(source=="w"):
        plt.suptitle("Averages: windowed moments (sims1) - "+str(frameIDs[0])+" "+st+" (mm)")
        
    plt.savefig(outPref+source+"pattern_"+str(frameIDs[0])+".png")
    

def loadPhys(loadPref,frameID):

    aa=np.load(loadPref+str(frameID)+".npy")
    
    
    t=float(aa[0])
    el=float(aa[1])
    inr=float(aa[2])
    za=90-int(15*round(el/15))
    inr=int(90 * round(inr/90))
    st="t={:.1f} za={:d} inr={:d}".format(t,za,inr)
    
    return st


def checkFiles(frameIDs,loadPref,source,check):

    print("fframe=",frameIDs[0])
    for frameID in frameIDs:
        fname=loadPref+"/"+source+"dump_"+str(frameID)+".npy"

        if(check=="y"):
            if(os.path.isfile(fname)):
                print(fname)
        else:
            if(not os.path.isfile(fname)):
                print(fname)
