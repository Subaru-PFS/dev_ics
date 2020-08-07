import math
import galsim
from galsim import des
import numpy as np
from scipy.stats import sigmaclip


def getFileNamesUC(frameId1,frameId2,frameSkip,sourceDir,fPref,dataType):

    """

    Generate a list of file names with complete path.

    the dataType flag toggles between different input sources
    (Taichung lab vs pinhole mask), and can be added to as needed.
    This will need to be updated when there are exposure + move ids.
    

    Input
       frameId1, frameId2 - first and last frames (inclusive)
       frameSkip - list of frames to be excluded
       sourceDir - full path of source directory
       fPref - file prefix (PFSC for MCS)
       dataType - 'pinhole' - commissioning run format (MCS, move id = 0)
                - 'taichung' - various lab data from testing
       
    Returns
       files = list of fileanmes (full path)
       prefix - prefix for output plot file names
       centroidFile - filename for saving centroids in local operation
       frameIDs - list of frameIds
        

    """
    
    #generate a list of frame IDs, delete any skipped ones
    
    frameIDs=list(np.arange(frameId1,frameId2+1))
    for ids in frameSkip:
        frameIDs.remove(ids)

    files=[]

    #different input data, assemle the file names

    if(dataType=='taichung'):
        for i in frameIDs:
            files.append(sourceDir+fPref+str(i).zfill(4)+"_uc.fits")

        prefix=fPref
        centroidFile=prefix+"_ncentroids.dat"

    elif(dataType=='pinhole'):
        for i in frameIDs:
            files.append(sourceDir+fPref+str(i).zfill(6)+"00_uc.fits")

        prefix="see_"+str(frameId1).zfill(5)+"_"+str(frameId2).zfill(5)
        
    centroidFile=prefix+"_centroidsSEX.dat"

    return files,prefix,centroidFile,frameIDs

def getFileNamesAllUC(parm1,parm2,parm3,frameSkip,sourceDir,fPref,dataType,tpe):
    """
    variation of filename routines
    """
    
    if(tpe==0):
        frameId1=parm1
        frameId2=parm2
        files,prefix,centroidFile,frameIDs=getFileNamesUC(frameId1,frameId2,frameSkip,sourceDir,fPref,dataType)

    elif(tpe==1):
        frameId1=parm1
        moveId1=parm2
        moveId2=parm3
        files,prefix,centroidFile,frameIDs=getFileNamesVisUC(frameId1,moveId1,moveId2,frameSkip,sourceDir,fPref,dataType)

    return files,prefix,centroidFile,frameIDs

def getFileNamesVisUC(frameId1,moveId1,moveId2,frameSkip,sourceDir,fPref,dataType):

    """
    variation of filename routines
    """

    frameIDs=frameId1*100+np.arange(moveId1,moveId2+1)
    files=[]

    #different input data, assemle the file names

    for i in frameIDs:
        files.append(sourceDir+fPref+str(i).zfill(8)+"_uc.fits")

        prefix="see_"+str(frameId1).zfill(6)+"_"+str(moveId1).zfill(2)+"_"+str(moveId1).zfill(2)
        
    centroidFile=prefix+"_centroidsSEX.dat"

    return files,prefix,centroidFile,frameIDs
       

def getFileNamesVis(frameId1,moveId1,moveId2,frameSkip,sourceDir,fPref,dataType):
    """
    variation of filename routines
    """

    frameIDs=frameId1*100+np.arange(moveId1,moveId2+1)
    files=[]

    #different input data, assemle the file names

    for i in frameIDs:
        files.append(sourceDir+fPref+str(i).zfill(8)+".fits")

        prefix="see_"+str(frameId1).zfill(6)+"_"+str(moveId1).zfill(2)+"_"+str(moveId1).zfill(2)
        
    centroidFile=prefix+"_centroids.dat"

    return files,prefix,centroidFile,frameIDs

def getFileNames(frameId1,frameId2,frameSkip,sourceDir,fPref,dataType):

    """

    Generate a list of file names with complete path.

    the dataType flag toggles between different input sources
    (Taichung lab vs pinhole mask), and can be added to as needed.
    This will need to be updated when there are exposure + move ids.
    

    Input
       frameId1, frameId2 - first and last frames (inclusive)
       frameSkip - list of frames to be excluded
       sourceDir - full path of source directory
       fPref - file prefix (PFSC for MCS)
       dataType - 'pinhole' - commissioning run format (MCS, move id = 0)
                - 'taichung' - various lab data from testing
       
    Returns
       files = list of fileanmes (full path)
       prefix - prefix for output plot file names
       centroidFile - filename for saving centroids in local operation
       frameIDs - list of frameIds
        

    """
    
    #generate a list of frame IDs, delete any skipped ones

    frameIDs=list(np.arange(frameId1,frameId2+1))
    for ids in frameSkip:
        frameIDs.remove(ids)

    files=[]

    #different input data, assemle the file names

    if(dataType=='taichung'):
        for i in frameIDs:
            files.append(sourceDir+fPref+str(i).zfill(4)+".fits")

        prefix=fPref
        centroidFile=prefix+"_ncentroids.dat"

    elif(dataType=='pinhole'):
        for i in frameIDs:
            files.append(sourceDir+fPref+str(i).zfill(6)+"00.fits")

        prefix="see_"+str(frameId1).zfill(5)+"_"+str(frameId2).zfill(5)
        
    centroidFile=prefix+"_centroids.dat"

    return files,prefix,centroidFile,frameIDs


def getFileNamesAll(parm1,parm2,parm3,frameSkip,sourceDir,fPref,dataType,tpe):
    """
    variation of filename routines
    """

    print(parm1,parm2,frameSkip)
    if(tpe==0):
        frameId1=parm1
        frameId2=parm2
        files,prefix,centroidFile,frameIDs=getFileNames(frameId1,frameId2,frameSkip,sourceDir,fPref,dataType)

    elif(tpe==1):
        frameId1=parm1
        moveId1=parm2
        moveId2=parm3
        files,prefix,centroidFile,frameIDs=getFileNamesVis(frameId1,moveId1,moveId2,frameSkip,sourceDir,fPref,dataType)

    return files,prefix,centroidFile,frameIDs


def simImage(sourceDir,imFile,catFile,psfFile,outFile):

    """ 

    Create a simulated image using PSFEx modelled PSFs and noise properties of the source image 

    Input: sourceDir: input directory data
           imFile: imput image file name
           catFile: catalogue file (output from SeXtractor)
           psfFile: psf model (output from PSFEx)
           outFile: name of output file for image

    Output: writes to fits file. 

    The catFile must contain the fields X_IMAGE, Y_IMAGE, FLUX_APER (or the code to 
    be changed to equivalent for positions of sources and integrated flux). 

    """
    #load necessary stuff from files.

    #NOte that the MCS image files have two HDUs, one with
    #the WCS information, one with the image information. 
    
    galHdr1 = galsim.FitsHeader(imFile, dir=sourceDir, hdu=0)
    galHdr2 = galsim.FitsHeader(imFile, dir=sourceDir, hdu=1)
    cat = galsim.Catalog(catFile, hdu=2, dir=sourceDir, file_type="FITS")
    psfex=des.DES_PSFEx(psfFile,imFile,dir=sourceDir)
    image=galsim.fits.read(imFile,sourceDir,hdu=1)
    
    #get setup the image. match the (currently trivial) WCS with the image, and
    #create a blank image
    
    wcs = galsim.FitsWCS(header=galHdr1)
    xSize=galHdr2['NAXIS1']
    ySize=galHdr2['NAXIS2']
    simImage = galsim.Image(xSize, ySize, wcs=wcs)

    #some definitions for extracting catalogue columsn
    xCol="X_IMAGE"
    yCol="Y_IMAGE"
    fluxCol="FLUX_APER"

    #get noise statistics. Read in the catalogue positions to estimate the centre
    #of the image in whatever rotation it has. This is so we get the noise statistics
    #from teh mask region, excludin the rest of the image. 
    
    xVals=cat.data[xCol]
    yVals=cat.data[yCol]
    xMean=int(xVals.mean())
    yMean=int(yVals.mean())
    radius=1800

    subIm=image.array[int(xMean-radius):int(xMean+radius),int(yMean-radius):int(yMean+radius)]
    im,a,b=sigmaclip(subIm.ravel(),5,5)

    skyLevel=im.mean()
    skySigma=im.std()
    gain = skyLevel / skySigma**2  #this definition from teh galsim tutorials
    nobj = cat.nobjects

    print('Catalog has ',nobj,' objects. Sky level is ',int(skyLevel),' Sky sigma is ',int(skySigma))

    #now cycle over the catalogue.
    
    for k in range(nobj):

        #get position and flux
        x = cat.getFloat(k,xCol)
        y = cat.getFloat(k,yCol)
        flux = cat.getFloat(k,fluxCol)*5

        #some position calculation for the galsim routines
        # + 0.5 to account for even-size postage stamps
        x=x+0.5
        y=y+0.5
        ix = int(math.floor(x+0.5))
        iy = int(math.floor(y+0.5))
        dx = x-ix
        dy = y-iy

        imagePos = galsim.PositionD(x,y)
        offset = galsim.PositionD(dx,dy)

        #calculate PSF for given position and flux
        psf=psfex.getPSF(imagePos).withFlux(flux)

        #make image
        stamp = psf.drawImage(wcs=wcs.local(imagePos), offset=offset, method='no_pixel')
        stamp.setCenter(ix,iy)

        #and place on image, taking into consideration edges
        bounds = stamp.bounds & simImage.bounds
        simImage[bounds] += stamp[bounds]

    #now that we've done all the spots, add noise

    #background sky level
    simImage += skyLevel

    #CCD noise
    random_seed = 1339201
    rng=galsim.BaseDeviate(random_seed)
    noise = galsim.CCDNoise(rng, gain=gain)

    #poisonnian noise
    simImage.addNoise(noise)
    noise = galsim.PoissonNoise(rng, sky_level=skyLevel)
    simImage.addNoise(noise)

    #and dump to a file. Will overwrite existing file.
    simImage.write(outFile,clobber=True)


def runSimImage(plist1,plist2,plist3,sourceDir,tpe,redo):
    

    delFiles=0
    redo=0
    dataType='pinhole'
    fPref="PFSC"
    delFiles=0
    
    for parm1,parm2,parm3 in zip(plist1,plist2,plist3):

       
        frameId1=parm1

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


        files,prefix,centroidFile,frameIDs=getFileNamesAll(parm1,parm2,parm3,frameSkip,"",fPref,dataType,tpe)
        filesUC,prefix,centroidFile,frameIDs=getFileNamesAllUC(parm1,parm2,parm3,frameSkip,"",fPref,dataType,tpe)


        for imFile,pref in zip(files,filesUC):
            psfFile=pref+"_cat1.psf"
            catFile=pref+"_cat3.dat"
            outFile="/Volumes/M5/"+imFile+"_sim.fits"
            simImage(sourceDir,imFile,catFile,psfFile,outFile)



#-------------------------------------------------------

redo=0
sourceDir="/Volumes/M5/Aug19Data/Day1"
tpe=0

fframe=[8450,8500,8550,8573,8623,8677,8727,8777,8843,8893,8943,8993,9043,9093,9148,9198,9248,9298,9348,9398,9448,9498,9548,9598,9648,9698,9748,9798,9848,9898,9948,9998,10048,10098,10148,10198,10248,10298,10348,10398,10448,10498,10548,10598,10648,10698,10748,10798,10848,10898,10948,10998,11048,11098,11148,11198,11248]

lframe=[8499,8549,8570,8622,8672,8726,8776,8836,8892,8942,8992,9042,9092,9142,9168,9247,9297,9347,9397,9447,9497,9547,9597,9647,9697,9747,9797,9847,9897,9947,9997,10047,10097,10147,10197,10247,10297,10347,10397,10447,10497,10547,10597,10647,10697,10747,10797,10847,10897,10947,10997,11047,11097,11147,11197,11247,11297]

fframe=[8843,8893,8943,8993,9043,9093,9148,9198,9248,9298,9348,9398,9448,9498,9548,9598,9648,9698,9748,9798,9848,9898,9948,9998,10048,10098,10148,10198,10248,10298,10348,10398,10448,10498,10548,10598,10648,10698,10748,10798,10848,10898,10948,10998,11048,11098,11148,11198,11248]

lframe=[8892,8942,8992,9042,9092,9142,9168,9247,9297,9347,9397,9447,9497,9547,9597,9647,9697,9747,9797,9847,9897,9947,9997,10047,10097,10147,10197,10247,10297,10347,10397,10447,10497,10547,10597,10647,10697,10747,10797,10847,10897,10947,10997,11047,11097,11147,11197,11247,11297]

#fframe=[10748,10798,10848,10898,10948,10998,11048,11098,11148,11198,11248]
#lframe=[10797,10847,10897,10947,10997,11047,11097,11147,11197,11247,11297]

runSimImage(fframe,lframe,lframe,sourceDir,tpe,redo)

#-------------------------------------------------------

sourceDir="/Volumes/M5/Aug19Data/Day2"


fframe=[11311,11352,11402,11452,11552,11602,11652,11702,11752,11802,11852,11902,11952,12002,12052,12102,12152,12200,12250,12300,12351,12401,12451,12501,12554,12604,12654,12711,12969,13019,13069,13119,13169,13219,13269,13319,13369,13419,13469,13519,13569,13619,13669,13719,13769,13819,13869,13919,13969,14019,14069,14119,14169,14219]

fframe=[11311,11352,11402,11452,11552,11602,11652,11702,11752,11802,11852,11902,11952,12002,12052,12102,12152,12200,12250,12300,12351,12401,12451,12501,12554,12604,12654,12711,12969,13019,13069,13119,13169,13219,13269,13319,13369,13419,13469,13519,13569,13619,13669,13819,13869,13919,13969,14019,14069,14119,14169,14219]


#fframe=[13619,13669,13819,13869,13919,13969,14019,14069,14119,14169,14219]


lframe=[11347,11401,11451,11501,11601,11651,11701,11751,11801,11851,11901,11951,12001,12051,12101,12151,12198,12249,12299,12349,12400,12450,12500,12550,12603,12653,12703,12760,13018,13068,13118,13168,13218,13268,13318,13368,13418,13468,13518,13568,13618,13668,13718,13868,13918,13968,14018,14068,14118,14168,14218,14268]

#lframe=[13668,13718,13868,13918,13968,14018,14068,14118,14168,14218,14268]


runSimImage(fframe,lframe,lframe,sourceDir,tpe,redo)

#-------------------------------------------------------

sourceDir="/Volumes/M5/Aug19Data/Day3"

fframe=[14301,14351,14401,14451,14501,14551,14601,14651,14701,14751,14801,14851,14901,14951,15001,15051,15101,15151,15201,15251,15301,15351,15401,15451,15501,15551,15601,15651,15701,15752]

lframe=[14350,14400,14450,14500,14550,14600,14650,14700,14750,14800,14850,14900,14950,15000,15050,15100,15150,15200,15250,15300,15350,15400,15450,15500,15550,15600,15650,15700,15750,15781]

runSimImage(fframe,lframe,lframe,sourceDir,tpe,redo)



#-------------------------------------------------------

sourceDir="/Volumes/M5/Aug19Data/Day3"
#16104
fframe=[15854,15864,15874,15884,15894,15904,15914,15924,15934,15944,15954,15964,15974,15984,15994,16004,16014,16024,16034,16044,16054,16064,16074,16084,16094,16104,16114,16124,16134,16144,16154,16164,16174,16184,16194,16204,16214,16224,16234,16244,16257,16304,16354,16404,16454]

lframe=[15863,15873,15883,15893,15903,15913,15923,15933,15943,15953,15963,15973,15983,15993,16003,16013,16023,16033,16043,16053,16063,16073,16083,16093,16103,16113,16123,16133,16143,16153,16163,16173,16183,16193,16203,16213,16223,16233,16243,16253,16303,16353,16403,16453,16503]

fframe=[16244,16257,16304,16354,16404,16454]

lframe=[16253,16303,16353,16403,16453,16503]


runSimImage(fframe,lframe,lframe,sourceDir,tpe,redo)

#-------------------------------------------------------

sourceDir="/Volumes/M5/Aug19Data/Day4"
#17493

fframe=[17530,17535,17540,17545,17550,17555,17560,17565,17570,17575,17580,17585]


lframe=[17534,17539,17544,17549,17554,17559,17564,17569,17574,17579,17584,17589]

runSimImage(fframe,lframe,lframe,sourceDir,tpe,redo)


#-------------------------------------------------------

tpe=1

sourceDir="/Volumes/M5/Aug19Data/Day5"
fframe=[17604,17605,17606,17607,17622,17623,17624,17625,17626,17627,17628,17629,17630,17631,17632,17633,17634,17640,17641]
mv1=[       0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0,    0]
mv2=[      49,   49,   49,   49,   49,   49,   49,   49,   49,   49,   49,   49,   49,   49,   49,   49,   49,   29,   29]


#fframe=[17630,17631,17632,17633]
#mv1=[0,0,0,0]
#mv2=[49,49,49,49]

#fframe=[17641]
#mv1=[0]
#mv2=[29]
runSimImage(fframe,mv1,mv2,sourceDir,tpe,redo)
