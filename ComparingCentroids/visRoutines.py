
import numpy as np
import numpy.ma as ma
from scipy.stats import sigmaclip
from astropy.io.fits import getdata
import cv2

try:
    import mcsActor.windowedCentroid.centroid as centroid
except:
    import centroid as centroid


     
def loadInstParams(config):

    """

    load instrument parameters. Dummy function, update!

    """

    
    if(config=='oct18'):
        rotCent=[4691.5,3095.7]
        offset=[0,-85]
    if(config=='aug19'):
        rotCent=[4470, 2873]
        offset=[0,-85]
        
    return rotCent,offset

def getImage(filename):

    """

    Simple wrapper to read an image from file

    """
    
    image=getdata(filename)

    return image
    
def toFits(filename):

    """

    Quick routine to convert raw image to fits.

    """
    
    a=np.fromfile(filename,dtype='uint16')
    image=a.reshape([5778,8960])

    pf.writeto(filename+".fits",image)

def getTransform(x,y,xx,yy,getVals):

    """

    given two sets of registered points, estimate the rigid transformation
    this is in a separate routine mostly for bookkeeping purposes. 
    Returns transformation matrix, and if getVals == 1 returns the 
    extracted parameters (rotation, translation, scale) as well. 

    input:
    x,y: input positions
    xx,yy: transformed positions
    getVales: if ==1, return parameters too

    output: 
    transformation: matrix 
    xd,yd: translations
    sx,sy: scalings
    rotation: rotation (radians)

    """

    #turn data into right form
    pts1=np.zeros((1,len(x),2))
    pts2=np.zeros((1,len(x),2))

    pts1[0,:,0]=x
    pts1[0,:,1]=y

    pts2[0,:,0]=xx
    pts2[0,:,1]=yy

    #float32 is needed
    pts1=np.float32(pts1)
    pts2=np.float32(pts2)

    #calculate the transformation
    transformation = cv2.estimateRigidTransform(pts1, pts2, False)

    afCoeff,inlier=cv2.estimateAffinePartial2D(pts1, pts2)
    
    if(getVals == 0):
        return transformation
    
    if(getVals == 1):

        #extract the parameters


        sx=np.sqrt(afCoeff[0,0]**2+afCoeff[0,1]**2)
        sy=np.sqrt(afCoeff[1,0]**2+afCoeff[1,1]**2)
        
        xd=afCoeff[0,2]
        yd=afCoeff[1,2]
         
        rotation=np.arctan2(afCoeff[1,0]/np.sqrt(afCoeff[0,0]**2+afCoeff[0,1]**2),
                                    afCoeff[1,1]/np.sqrt(afCoeff[1,0]**2+afCoeff[1,1]**2))

        
        #sx=np.sqrt(transformation[0,0]**2+transformation[0,1]**2)
        #sy=np.sqrt(transformation[1,0]**2+transformation[1,1]**2)
        #xd=transformation[0,2]
        #yd=transformation[1,2]
        #rotation = np.arctan2(transformation[1,0]/sx,transformation[1,1]/sy)

        return afCoeff,xd,yd,sx,sy,rotation

def transformPointsNew(x,y,xd,yd,theta,sx,sy):

    """
    Apply a rigid transformation to the mask (trans, scale, rot). Mostly bookkeeping
    stuff. 

    input:
    x,y: mask positions
    xd,yd: translation
    theta: rotation (radians)
    s: scale

    output: transformed x,y

    """
    
    #create transformation matrix
    matrix=np.zeros((2,3))
    matrix[0,0]=np.cos(theta)*sx
    matrix[0,1]=-np.sin(theta)*sy
    matrix[1,0]=np.sin(theta)*sx
    matrix[1,1]=np.cos(theta)*sy
    matrix[0,2]=xd
    matrix[1,2]=yd

    #bookkeeping for coordinate format
    pts=np.zeros((1,len(x),2))

    pts[0,:,0]=x
    pts[0,:,1]=y

    pts=np.float32(pts)

    #do the transform
    pts1=cv2.transform(pts,matrix)

    #more bookkeeping
    xx=pts1[0,:,0]
    yy=pts1[0,:,1]

    return xx,yy
        

