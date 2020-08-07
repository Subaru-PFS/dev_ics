"""
Populate the database, reading data from locally created numpy files. 

There are switches depending on whether the frameid is one per image, 
or a set of images for each frameid, one for each move, as this was 
implemented part way through the run. 

"""


import numpy as np
import psycopg2
import fileNames as fr
from datetime import datetime

def doComm(cur,comm,inter):

    """

    debugging routine for testing sql command creation.

    Args:
    cur: cursor object connected to database
    comm: string containing the sql command string
    inter: flag. inter=1 executes the command, =1 print to stdout and execute, otherwise just print
    
    """
    
    if(inter == 1):
        cur.execute(comm)
    elif(inter == 2):
        print(comm)
        cur.execute(comm)
    else:
        print(comm)

def insertVisit(frameID,cur):

    """
    The visit insertion is simple. 

    """

    comm="""insert into pfs_visit (pfs_visit_id) VALUES ("""+str(int(frameID // 100))+""")"""
    cur.execute(comm)
        
def readPhys(frameID,loadPref):

    """
    Retrieve telemetry, exposure time, etc. from numpy file (replaced by database equivalent).
    Accounts for the effects of midnight and the date change. 

    Input: frameID
           loadPref: prefix for file
    """

    if(frameID < 1760000):
        ff=np.load(loadPref+"/nphys_"+str(int(frameID/100))+".npy")
    else:
        ff=np.load(loadPref+"/nphys_"+str(int(frameID))+".npy")

    phys={}

    phys["pfs_visit_id"]=str(int(frameID // 100))
    phys["mcs_exptime"]=ff[0]
    phys["altitude"]=ff[1]
    phys["azimuth"]=ff[2]
    phys["insrot"]=ff[3]
    phys["adc_pa"]=ff[4]
    phys["dome_pressure"]=ff[5]
    phys["dome_temperature"]=ff[6]
    phys["dome_humidity"]=ff[7]
    phys["outside_temperature"]=ff[8]
    phys["outside_pressure"]=ff[9]
    phys["outside_humidity"]=ff[10]
    phys["mcs_cover_temperature"]=ff[11]
    phys["mcs_m1_temperature"]=ff[12]

    #format the date string (UT)
    phys['taken_at']="'"+ff[15]+" "+ff[13]+"'"

    obs1=ff[14].split(":")
    obs=ff[15].split("-")
    #HST has  a different date than UT if before midnight
    if(float(obs1[0]) > 10):
        obs[2]=str(int(obs[2])-1)

    phys["taken_in_hst_at"]="'"+obs[0]+"-"+obs[1]+"-"+obs[2]+" "+ff[14]+"'"


    return phys

def readFPS(frameID,loadPref):

    """
    read FPS data from numpy files; a bit of massaging to get them in the right form. 
    """

    if(frameID < 1760000):
        ff=np.load(loadPref+"/ldump_"+str(int(frameID/100))+".npy")
    else:
        ff=np.load(loadPref+"/ldump_"+str(int(frameID))+".npy")

    sz=ff.shape
    
    fps={}

    if(frameID < 1760000):
        fps["mcs_frame_id"]=np.repeat(int(frameID),sz[0])
        fps["iteration"]=np.repeat(0,sz[0])
    else:
        fps["mcs_frame_id"]=np.repeat((frameID),sz[0])
        fps["iteration"]=np.repeat(frameID % 100,sz[0])
        
    fps["spot_id"]=np.where(np.isnan(ff[:,3].ravel()),None,ff[:,3].ravel())
    fps["cobra_id"]=ff[:,0].astype('int')

    fps["pfi_center_x_mm"]=np.where(np.isnan(ff[:,1].ravel()),None,ff[:,1].ravel())
    fps["pfi_center_y_mm"]=np.where(np.isnan(ff[:,2].ravel()),None,ff[:,2].ravel())

    
    
    return fps,sz
    
    
    
def insertFPS(frameID,loadPref,cur):

    """
    insert fps data into database. 
    """
    
    table="cobra_status"

    fps,sz=readFPS(frameID,loadPref)

    for i in range(sz[0]):
        
        comm="""insert into cobra_status ("""

        for key in fps:
            comm=comm+key+""","""
        comm=comm[:-1]+""") VALUES ("""
        #comm=comm+str(int(frameID))+","


        for key in fps:
            if(fps[key][i]==None):
                comm=comm+"NULL"+""","""
            else:
                comm=comm+str(fps[key][i])+""","""
        comm=comm[:-1]+""")"""
        cur.execute(comm)
    
def readCentroids(frameID,loadPref):

    """
    read centroid data from numpy files; a bit of massage for the rihgt format. 
    """

    if(frameID < 1760000):
        ff=np.load(loadPref+"/lcent_"+str(int(frameID/100))+".npy")
    else:
        ff=np.load(loadPref+"/lcent_"+str(int(frameID))+".npy")

    
    #ff=np.load(loadPref+"/lcent_"+str(frameID)+".npy")
    sz=ff.shape

    #fps["spot_id"]=np.where(np.isnan(ff[:,3].ravel()),None,ff[:,3].ravel())

    cents={}

    if(frameID < 1760000):
        cents['mcs_frame_id']=ff[:,0].ravel().astype('int')*100
    else:
        cents['mcs_frame_id']=ff[:,0].ravel().astype('int')
    cents['spot_id']=np.arange(sz[0]).astype('int')
    cents['mcs_center_x_pix']=np.where(np.isnan(ff[:,1].ravel()),None,ff[:,1].ravel())
    cents['mcs_center_y_pix']=np.where(np.isnan(ff[:,2].ravel()),None,ff[:,2].ravel())
    cents['mcs_second_moment_x_pix']=np.where(np.isnan(ff[:,3].ravel()),None,ff[:,3].ravel())
    cents['mcs_second_moment_y_pix']=np.where(np.isnan(ff[:,4].ravel()),None,ff[:,4].ravel())
    cents['peakvalue']=np.where(np.isnan(ff[:,5].ravel()),None,ff[:,5].ravel())
    cents['mcs_second_moment_xy_pix']=np.where(np.isnan(ff[:,6].ravel()),None,ff[:,6].ravel())
    cents['bgvalue']=np.where(np.isnan(ff[:,7].ravel()),None,ff[:,7].ravel())

    return cents,sz

    
def insertCentroids(frameID,loadPref,cur):

    """
    write centroids to database. 
    """

    table="mcs_data"

    cents,sz=readCentroids(frameID,loadPref)

    for i in range(sz[0]):
        
        comm="""insert into mcs_data ("""

        for key in cents:
            comm=comm+key+""","""
        comm=comm[:-1]+""") VALUES ("""
        #comm=comm+str(int(frameID))+","
        for key in cents:
            if(cents[key][i]==None):
                comm=comm+"NULL"+""","""
            else:
                comm=comm+str(cents[key][i])+""","""
        comm=comm[:-1]+""")"""
        cur.execute(comm)

def insertPhys(frameID,loadPref,cur):

    """
    insert telescope information to database
    """

    phys=readPhys(frameID,loadPref)

    comm="""insert into mcs_exposure ("""

    comm=comm+"mcs_frame_id,"
    for key in phys:
        comm=comm+key+""","""
    comm=comm[:-1]+""") VALUES ("""

    comm=comm+str(int(frameID))+","
    for key in phys:
        comm=comm+phys[key]+""","""
    comm=comm[:-1]+""")"""
    cur.execute(comm)

def insertSet(plist1,plist2,plist3,sourceDir,tpe,loadPref,switch,cur):

    """
    wrapper to insert a set of data. 
    """
    
    delFiles=0
    redo=0
    dataType='pinhole'
    fPref="PFSC"

    for parm1,parm2,parm3 in zip(plist1,plist2,plist3):
        #bookkeeping - a few sets have missing frames
        frameId1=parm1
        print(frameId1)
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

        files,filesUC,prefix,centroidFile,frameIDs=fr.getFileNamesAll(parm1,parm2,parm3,frameSkip,sourceDir,fPref,dataType,tpe)

        
        if(frameID[0] < 1760000):
            for frameID in frameIDs:

                frameID=frameID*100
                if(switch=="phys"):
                    insertPhys(frameID,loadPref,cur)
                if(switch=="cent"):
                    insertCentroids(frameID,loadPref,cur)
                if(switch=="fps"):
                    insertFPS(frameID,loadPref,cur)
                if(switch=="visit"):
                    insertVisit(frameID,cur)
        else:
            if(switch=="visit"):
                insertVisit(frameIDs[0],cur)
            if(switch=="phys"):
                insertPhys(frameIDs[0],loadPref,cur)
            if(switch=="cent"):
                for frameID in frameIDs:
                    insertCentroids(frameID,loadPref,cur)
            if(switch=="fps"):
                for frameID in frameIDs:
                    insertFPS(frameID,loadPref,cur)
