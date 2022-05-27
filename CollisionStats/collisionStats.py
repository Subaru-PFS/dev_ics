import pfs.utils.coordinates.transform as transformUtils
import matplotlib.pylab as plt
import pandas as pd
import numpy as np
from opdb import opdb
from scipy import spatial

def collQuery(db,mcsFrameId,cameraName):

    """
    pull the necessary info from the DB
    """
    
    pfsVisitId = mcsFrameId // 100
    iteration = mcsFrameId % 100
    
    if('rmod' in cameraName):
        insrot = 0
        altitude = 90
    else:
        sql=f'select altitude, insrot from mcs_exposure where mcs_frame_id={mcsFrameId}'
        df = db.fetch_query(sql)

        altitude=df['altitude'].values[0]
        insrot=df['insrot'].values[0]
    
    sql=f'select mcs_frame_id,spot_id,mcs_center_x_pix,mcs_center_y_pix from mcs_data where mcs_frame_id={mcsFrameId}'
    dfPix = db.fetch_query(sql)
    
    sql = f'select pfs_visit_id*100+iteration as mcs_frame_id, cobra_id, pfi_center_x_mm, pfi_center_y_mm from cobra_match where pfs_visit_id={pfsVisitId} and iteration={iteration}'
    dfMM = db.fetch_query(sql)
    
    return dfPix,dfMM,altitude,insrot

def transformFiducials(mcsData,fids,cameraName,altitude,insrot):

    """
    match and transform the fiducials - we use transformation rather than the true values 
    so that any variations in the transformations are uniform across the data set (as we're
    looking for close neighbours)

    """
    
    pfiTransform = transformUtils.fromCameraName(cameraName, 
                altitude=altitude, insrot=insrot,nsigma=0, alphaRot=1)

    # get outer outer ring and good fibres list
    outerRingIds = [29, 30, 31, 61, 62, 64, 93, 94, 95, 96]
    fidsOuterRing = fids[fids.fiducialId.isin(outerRingIds)]
    badFids = [1,32,34,61,68,75,88,89,2,4,33,36,37,65,66,67,68,69]
    fidList=list(fids['fiducialId'].values)
    for i in badFids:
        try:
            fidList.remove(i)
        except:
            pass

    goodFid = np.zeros(len(fids), dtype=bool)   
    for i in fidList:
        goodFid[fids.fiducialId == i] = True

    # update transform
    ffids0,dist0=pfiTransform.updateTransform(mcsData, fidsOuterRing, matchRadius=8.0, nMatchMin=0.1)
    ffids1,dist1=pfiTransform.updateTransform(mcsData, fids[goodFid], matchRadius=4.2,nMatchMin=0.1)
    ffids1,dist1=pfiTransform.updateTransform(mcsData, fids[goodFid], matchRadius=2,nMatchMin=0.1)

    # transform f
    
    xMM,yMM = pfiTransform.mcsToPfi(mcsData['mcs_center_x_pix'].values,mcsData['mcs_center_y_pix'].values)

    #match the fids
    fidXmm,fidYmm = matchFids(xMM,yMM,fids['x_mm'].values,fids['y_mm'].values)
    return fidXmm, fidYmm,xMM,yMM

def matchFids(xP,yP,xF,yF):

    """
    match fiducials and measured points
    """
    
    mX = []
    mY = []
    for i in range(len(xF)):
        dd=np.sqrt((xF[i]-xP)**2+(yF[i]-yP)**2)
        ind=np.nanargmin(dd)
        if(dd[ind] < 2):
            mX.append(xP[ind])
            mY.append(yP[ind])
        else:
            mX.append(0)
            mY.append(0)
    return mX,mY
    
def checkFidCollision(x,y,fidX,fidY,collRad):
    
    #find potential collisions between cobras and fiducials
    #use cKDTree for speed
    
    # create a tree based on the spot positions: the np.c_ is needed for contiguous memory
    # allocation, which is needed to run the code properly
    
    pp = np.c_[x.ravel(), y.ravel()]
    tree = spatial.cKDTree(pp)
    
    # query the tree for  fiducial positions within collRad of any spot positioon
    # returns a list of lists, one for each point
    potColl=tree.query_ball_point(np.c_[fidX,fidY],collRad)
    
    return potColl

def checkCobraCollision(x,y,collRad):
    
    # find potential collisions between cobras and fiducials. 
    # almost the same as checkFidCollision, except we need to exclude 
    # points matching with themselves
    #use cKDTree for speed

    # create a tree based on the spot positions: the np.c_ is needed for contiguous memory
    # allocation, which is needed to run the code properly
    
    pp = np.c_[x.ravel(), y.ravel()]
    tree = spatial.cKDTree(pp)

    # query the tree for  fiducial positions within collRad of any spot positioon
    # returns a list of lists, one for each point
    potColl=tree.query_ball_point(pp,collRad)

    # remove the self match from each query
    for i in range(len(potColl)):
        potColl[i].remove(i)
        
    return potColl


def collCheck(db,mcsFrameId,cameraName,fids):

    """
    run the sequence of commands

    """
    
    dfPix,dfMM,altitude,insrot = collQuery(db,7275200,cameraName)
    fidXmm, fidYmm, xMM, yMM =transformFiducials(dfPix,fids,cameraName,altitude,insrot)
    potColl1=checkFidCollision(dfMM['pfi_center_x_mm'],dfMM['pfi_center_y_mm'],fidXmm,fidYmm,1.5)
    potColl2=checkCobraCollision(dfMM['pfi_center_x_mm'],dfMM['pfi_center_y_mm'],1.5)

    return potColl1,potColl2
