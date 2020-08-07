
from astropy.io import fits
import numpy as np
import os

import mcsRoutines as mcs
import fpsRoutines as fps
import visRoutines as vis
import astromatic_wrapper as aw

import analysisRoutines as ar
import fileNames as fn
import centroid as centroid

def getParameters(parm1,parm2,sourceDir,tpe,dataType,suffix):

    """

    All the various parameters needed for the set - telescope information, predicted postion of science and fiducial fibres,
    filenames, frameIDs. 

    """
    
    #for file names
    fPref="PFSC"
    inter=0
    cenFlag='local'
    config="aug19"

    
    frameId1=parm1

        
    if(frameId1==9298):
        frameSkip=[9298]
    if(frameId1==8777):
        frameSkip=[8805,8806]
    elif(frameId1==11602):
        frameSkip=[11620]
    elif(frameId1==13219):
        frameSkip=[13247]
    elif(frameId1==16831):
        frameSkip=[16853,16852,16851,16850,16849]   
    else:
        frameSkip=[]


    #get filenames
    files,frameIDs=fn.getFileNamesAll(parm1,parm2,frameSkip,sourceDir,fPref,dataType,tpe,suffix)

    #get telescope information
    ff=ar.readPhysTidy(frameIDs[0],"Phys")
    za=ff['za']
    inr=ff['inr']
    t=ff['t']

    #get instrument parameters
    rotCent,offset=vis.loadInstParams('aug19')
    
    #get spot positions in PFI mm
    fiducials,scienceFibres=fps.getFieldDefinition('')

    #fibres=fps.getAllFibrePos(fiducials,scienceFibres,za,inr,rotCent,offset)
    fibres,allPos=fps.getAllFibrePos1(fiducials,scienceFibres,za,inr,rotCent,offset)
    
    #get the geometry for a given image
    fidPos,sciPos=fps.getFibrePos(fiducials,scienceFibres,za,inr,rotCent,offset)


    return fidPos,sciPos,fiducials,scienceFibres,za,inr,rotCent,offset,frameIDs,fibres,allPos,files


def getAllCentroids(files,thresh1,thresh2,fwhmx,fwhmy,boxFind,boxCent,nmin,nmax,maxIt,frameIDs,loadPref,source):        

    """

    wrapper to do the centroiding for a set of files. Writes to a text
    output file.  Returns the last set of xy positions (used for the
    next routine.

    input: 
    
    files: list of FITS files
    outfile: name of output text file
    FWHM, boxsize, thresh, rl,rh,sl,sh: parameters for centroiding routine

    output: writes to a file with each line of the format
       framenum, x, y, fx, fy, peak, back, qual

    """

    nfiles=len(files)
    
    #cycle through files
    filenum=0

    print(str(nfiles)+" Frames. Centroiding ",end="")
    for (fname,fid) in zip(files,frameIDs):

        #read in image
        image=vis.getImage(fname)
        print("got image")
        
        #centroid and reshape the results
    
        a=centroid.centroid_only(image.astype('<i4'),fwhmx,fwhmy,thresh1,thresh2,boxFind,boxCent,nmin,nmax,maxIt,0)

        centroids=np.frombuffer(a,dtype='<f8')
        centroids=np.reshape(centroids,(len(centroids)//7,7))

        npoint=centroids.shape[0]

        tCentroids = np.zeros((npoint,8))
        tCentroids[:,1:]=centroids
        tCentroids[:,0]=np.zeros((npoint))+fid

        np.save(loadPref+"/"+source+"cent_"+str(fid)+".npy",tCentroids,allow_pickle=False)
        filenum+=1
            
def runCentroids(parm1,parm2,sourceDir,tpe,redo,loadPref,source,suffix,boxCent=None):

    """

    Runs the centroid routine and does the matching

    """

    dataType="pinhole"

    sigmaFind=50   #threshold for finding spots
    sigmaCent=15   #threshold for calculating moments
    boxFind=10   #box size for finding spots
    nmin=10     #minimum acceptable # of pixels in spot
    nmax=90     #maximum acceptable # of pixels in spot
    maxIt=20    #maximum # interations for centroiding
    sigmaThresh=4
    threshFact=2
    fwhmx=0
    fwhmy=0
    fact=1

    if(boxCent==None):
        boxCent=4   #box size for centroiding

    frameId1=parm1
        
    if(frameId1==9298):
        frameSkip=[9298]
    if(frameId1==8777):
        frameSkip=[8805,8806]
    elif(frameId1==11602):
        frameSkip=[11620]
    elif(frameId1==13219):
        frameSkip=[13247]
    elif(frameId1==16831):
        frameSkip=[16853,16852,16851,16850,16849]   
    else:
        frameSkip=[]

    fidPos,sciPos,fiducials,scienceFibres,za,inr,rotCent,offset,frameIDs,fibres,allPos,files=getParameters(parm1,parm2,sourceDir,tpe,dataType,suffix)

    if(redo==1):

        
        image=vis.getImage(files[0])
        findThresh,centThresh,xrange,yrange=mcs.getThresh(image,'fieldID',sigmaThresh,threshFact,sigmaFind,sigmaCent,fibres)
        getAllCentroids(files,findThresh*fact,centThresh,fwhmx,fwhmy,boxFind,boxCent,nmin,nmax,maxIt,frameIDs,loadPref,source)

    ii=0
    for fid in frameIDs:
        getPhys(files[i],fid)

        #load the centroids
        centroids=np.load(loadPref+"/"+source+"cent_"+str(int(fid))+".npy")

        #radius for matching
        tol=50
            
        #assemble array
        cc=np.array([range(len(centroids[:,1].ravel())),centroids[:,1].ravel(),centroids[:,2].ravel(),centroids[:,3].ravel(),centroids[:,4].ravel(),centroids[:,5].ravel(),centroids[:,6].ravel(),centroids[:,7].ravel()]).T
            
        #transform centroids to mm
        ccTrans=fps.transformCentroids(cc,za,inr,rotCent,offset)

        #match fiducials
        fCentroid=mcs.findHomes(ccTrans,fiducials,tol)
            
        #get affine transformation
        trans,xd,yd,sx,sy,rot=vis.getTransform(fCentroid[:,1],fCentroid[:,2],fiducials[:,1],fiducials[:,2],1)
        #apply affine
        xx1,yy1=vis.transformPointsNew(ccTrans[:,1],ccTrans[:,2],xd,yd,rot,sx,sy)

        #match science fibres
        tol=6
        
        aCent=mcs.findHomes(np.array([cc[:,0],xx1,yy1,cc[:,0],cc[:,1],cc[:,2],cc[:,3],cc[:,4],cc[:,5],cc[:,6],cc[:,7]]).T,scienceFibres,tol)
        dx,dy = fps.getDiff(aCent,scienceFibres)

        #dump to files
        dumpArray=np.array([aCent[:,0], #number
                            aCent[:,1], #x mm
                            aCent[:,2], #y mm
                            aCent[:,3], #number (original)
                            aCent[:,4], #x (pix)
                            aCent[:,5], #y (pix)
                            aCent[:,6], #fx
                            aCent[:,7], #fy
                            aCent[:,8], #peak
                            aCent[:,9], #ellip
                            aCent[:,10],
                            scienceFibres[:,1],
                            scienceFibres[:,2],
                            dx,dy,]).T

        np.save(loadPref+"/"+source+"dump_"+str(fid)+".npy",dumpArray,allow_pickle=False)
        np.save(loadPref+"/"+source+"trans_"+str(fid)+".npy",trans,allow_pickle=False)
        ii=ii+1

def getPhys(fname,frameID):

    hdu=fits.open(fname)
    hdr = hdu[0].header
    
    t=hdr['exptime']
    el=hdr['altitude']
    inr=hdr['inr-str']
    adc=hdr['adc-str']
    date=hdr['date-obs']
    ut=hdr['ut']
    hst=hdr['hst']
    hum=hdr['dom-hum']
    tmp=hdr['dom-tmp']
    prs=hdr['dom-prs']
    wnd=hdr['dom-wnd']
    mcm1t=hdr['w_mcm1t']
    mctopt=hdr['w_mctopt']
    mccftt=hdr['w_mccftt']
    mccovt=hdr['w_mccovt']
    mccint=hdr['w_mccint']
    mccott=hdr['w_mccott']
    mcelet=hdr['w_mcelet']
    mcflow=hdr['w_mcflow']
    ohum=hdr['out-hum']
    otmp=hdr['out-tmp']
    oprs=hdr['out-prs']
    ownd=hdr['out-wnd']

    dumpArray=[t,el,inr,adc,date,ut,hst,hum,tmp,prs,ohum,otmp,oprs,ownd,mcm1t,mctopt,mccftt,mccovt,mccint,mccott,mcelet,mcflow]
    np.save("Phys/nphys_"+str(frameID)+".npy",dumpArray,allow_pickle=False)

def runGetPhys(plist1,plist2,sourceDir,tpe,redo):

    for parm1,parm2 in zip(plist1,plist2):
        #print(parm1,parm2,parm3,tpe,sourceDir)
        frameIDs,centroidFile,files=ar.getFrameIDs(parm1,parm2,parm3,tpe,sourceDir)
        print(frameIDs)
        for fname,frameID in zip(files,frameIDs):
            getPhys(fname,frameID)

