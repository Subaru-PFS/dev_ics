
import numpy as np

def getFileNames(frameId1,frameId2,frameSkip,sourceDir,fPref,dataType,suffix):

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

    #different input data, assemle the file names. Each has two cases, for before and after moveid

    #basic file name. 
    if(dataType=='pinhole'):

        if(frameIDs[0] < 1700000): 
            for i in frameIDs:
                files.append(sourceDir+"/"+fPref+str(i).zfill(6)+"00.fits")
        else:
            for i in frameIDs:
                files.append(sourceDir+"/"+fPref+str(i).zfill(8)+".fits")
                
    #uncompressed files (for sextractor)
    elif(dataType=='uncomp'):
        if(frameIDs[0] < 1700000): 
            for i in frameIDs:
                files.append(sourceDir+"/"+fPref+str(i).zfill(6)+"00_uc.fits")
        else:
            for i in frameIDs:
                files.append(sourceDir+"/"+fPref+str(i).zfill(8)+"_uc.fits")
            
    #simulated version 1        
    elif(dataType=='sim1'):
        if(frameIDs[0] < 1700000): 
            for i in frameIDs:
                files.append(sourceDir+"/"+fPref+str(i).zfill(6)+"00.fits_"+suffix+".fits")
        else:
            for i in frameIDs:
                files.append(sourceDir+"/"+fPref+str(i).zfill(8)+".fits_"+suffix+".fits")

    elif(dataType=="sim2"):
        ii=0
        if(frameIDs[0] < 1700000): 
            for i in frameIDs:
                files.append(sourceDir+"/"+fPref+str(frameIDs[0]).zfill(6)+"00.fits_"+suffix+str(int(ii)).zfill(2)+".fits")
                ii=ii+1

        else:
            for i in frameIDs:
                files.append(sourceDir+"/"+fPref+str(frameIDs[0]).zfill(8)+".fits_"+suffix+str(int(ii)).zfill(2)+".fits")
                ii=ii+1

    return files,frameIDs


def getFileNamesAll(parm1,parm2,frameSkip,sourceDir,fPref,dataType,tpe,suffix):
    
    frameId1=parm1
    frameId2=parm2
    files,frameIDs=getFileNames(frameId1,frameId2,frameSkip,sourceDir,fPref,dataType,suffix)
    return files,frameIDs
