
"""
Run through the whole set of data, make sure the data are there for each set, and have the 
expected dimensions. 

"""

import numpy as np
import psycopg2



import dbStuff as db
import pandas as pd
import pandas.io.sql as sqlio
import dbStuff as db
import fileNames as fr
#9148-9168 instead of 9148-9197

def runAll(switch,cur):
    
    loadPref="dbSet"

    sourceDir="/Volumes/M5/Aug19Data/Day1"
    day=1
    tpe=0
        
    fframe=[8450,8500,8550,8573,8623,8677,8727,8777,8843,8893,8943,8993,9043,9093,9148,9198,9248,9298,9348,9398,9448,9498,9548,9598,9648,9698,9748,9798,9848,9898,9948,9998,10048,10098,10148,10198,10248,10298,10348,10398,10448,10498,10548,10598,10648,10698,10748,10798,10848,10898,10948,10998,11048,11098,11148,11198,11248]


    lframe=[8499,8549,8570,8622,8672,8726,8776,8836,8892,8942,8992,9042,9092,9142,9168,9247,9297,9347,9397,9447,9497,9547,9597,9647,9697,9747,9797,9847,9897,9947,9997,10047,10097,10147,10197,10247,10297,10347,10397,10447,10497,10547,10597,10647,10697,10747,10797,10847,10897,10947,10997,11047,11097,11147,11197,11247,11297]

    #fframe=[8677]
    #lframe=[8726]
   
    checkSet(fframe,lframe,lframe,sourceDir,tpe,loadPref,switch,cur)
    
    #---------------------------------------------------------------------
    
    sourceDir="/Volumes/M5/Aug19Data/Day2"
    day=2
    tpe=0

    fframe=[11311,11352,11402,11452,11552,11602,11652,11702,11752,11802,11852,11902,11952,12002,12052,12102,12152,12200,12250,12300,12351,12401,12451,12501,12554,12604,12654,12711,12969,13019,13069,13119,13169,13219,13269,13319,13369,13419,13469,13519,13619,13669,13569,13819,13869,13919,13969,14019,14069,14119,14169,14219]

    lframe=[11347,11401,11451,11501,11601,11651,11701,11751,11801,11851,11901,11951,12001,12051,12093,12151,12198,12249,12299,12349,12400,12450,12500,12550,12603,12653,12703,12760,13018,13068,13118,13168,13218,13268,13318,13368,13418,13468,13518,13568,13668,13718,13618,13868,13918,13968,14018,14068,14118,14168,14218,14268]

    #fframe=[13819,13869,13919,13969,14019,14069,14119,14169,14219]

    #lframe=[13868,13918,13968,14018,14068,14118,14168,14218,14268]

    
    #fframe=[12052,14219]
    #lframe=[12093,14268]
    
    #db.insertSet(fframe,lframe,lframe,sourceDir,tpe,loadPref,switch,cur)
    checkSet(fframe,lframe,lframe,sourceDir,tpe,loadPref,switch,cur)

    #---------------------------------------------------------------------

    sourceDir="/Volumes/M5/Aug19Data/Day3"
    day=3
    tpe=0

    fframe=[14301,14351,14401,14451,14501,14551,14601,14651,14701,14751,14801,14851,14901,14951,15001,15051,15101,15151,15201,15251,15301,15351,15401,15451,15501,15551,15601,15651,15701,15752,15854,15864,15874,15884,15894,15904,15914,15924,15934,15944,15954,15964,15974,15984,15994,16004,16014,16024,16034,16044,16054,16064,16074,16084,16094,16104,16114,16124,16134,16144,16154,16164,16174,16184,16194,16204,16214,16224,16234,16244,16257,16304,16354,16404,16454]

    lframe=[14350,14400,14450,14500,14550,14600,14650,14700,14750,14800,14850,14900,14950,15000,15050,15100,15150,15200,15250,15300,15350,15400,15450,15500,15550,15600,15650,15700,15750,15781,15863,15873,15883,15893,15903,15913,15923,15933,15943,15953,15963,15973,15983,15993,16003,16013,16023,16033,16043,16053,16063,16073,16083,16093,16103,16113,16123,16133,16143,16153,16163,16173,16183,16193,16203,16213,16223,16233,16243,16253,16303,16353,16403,16453,16503]

    #db.insertSet(fframe,lframe,lframe,sourceDir,tpe,loadPref,switch,cur)
    checkSet(fframe,lframe,lframe,sourceDir,tpe,loadPref,switch,cur)

    #---------------------------------------------------------------------

    sourceDir="/Volumes/M5/Aug19Data/Day4"
    day=4
    tpe=0

    
    fframe=[16831,16885,16935,17118,17123,17128,17133,17140,17145,17150,17155,17160,17165,17170,17175,17180,17185,17190,17195,17200,17205,17210,17215,17220,17225,17230,17235,17240,17245,17250,17255,17265,17270,17275,17292,17297,17302,17307,17312,17317,17322,17327,17337,17347,17357,17362,17367,17372,17377,17382,17387,17392,17397,17402,17407,17412,17417,17422,17427,17432,17437,17442,17447,17452,17457,17462,17467,17472,17478,17483,17488,17493,17530,17535,17540,17545,17550,17555,17560,17565,17570,17575,17580,17585]


    lframe=[16884,16934,16984,17122,17127,17132,17137,17144,17149,17154,17159,17164,17169,17174,17179,17184,17189,17194,17199,17204,17209,17214,17219,17224,17229,17234,17239,17244,17249,17254,17259,17269,17274,17279,17296,17301,17306,17311,17316,17321,17326,17331,17346,17356,17361,17366,17371,17376,17381,17386,17391,17396,17401,17406,17411,17416,17421,17426,17431,17436,17441,17446,17451,17456,17461,17466,17471,17477,17482,17487,17492,17497,17534,17539,17544,17549,17554,17559,17564,17569,17574,17579,17584,17589]

    
    #db.insertSet(fframe,lframe,lframe,sourceDir,tpe,loadPref,switch,cur)
    checkSet(fframe,lframe,lframe,sourceDir,tpe,loadPref,switch,cur)

    #---------------------------------------------------------------------

    sourceDir="/Volumes/M5/Aug19Data/Day5"
    day=5
    tpe=0

    fframe=[1760400,1760500,1760600,1760700,1762200,1762300,1762400,1762500,1762600,1762700,1762800,1762900,1763000,1763100,1763200,1763300,1763400,1764000,1764100]

    lframe=[1760449,1760549,1760649,1760749,1762249,1762349,1762449,1762549,1762649,1762749,1762849,1762949,1763049,1763149,1763249,1763349,1763449,1764029,1764129]


    
    #fframe=[1763400]
    #lframe=[1763440]
    mv1=[0]
    mv2=[49]
    #tpe=1
    #db.insertSet(fframe,lframe,lframe,sourceDir,tpe,loadPref,switch,cur)
    checkSet(fframe,lframe,lframe,sourceDir,tpe,loadPref,switch,cur)

    #---------------------------------------------------------------------

    
def checkSet(plist1,plist2,plist3,sourceDir,tpe,loadPref,switch,cur):
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

        
        if(frameIDs[0] < 1760000):
            for frameID in frameIDs:
                frameID=frameID*100

                if(switch=="visit"):
                    comm="select * from pfs_visit where pfs_visit_id = """+str(int(frameID // 100))
                    data = sqlio.read_sql_query(comm, con)
                    if(data.shape != (1,2)):
                        print(frameID,data.shape)
                elif(switch=="phys"):
                    comm="select * from mcs_exposure where mcs_frame_id = """+str(int(frameID))
                    data = sqlio.read_sql_query(comm, con)
                    if(data.shape != (1,17)):
                        print(frameID,data.shape)
                elif(switch=="cent"):
                    comm="select * from mcs_data where mcs_frame_id = """+str(int(frameID))
                    data = sqlio.read_sql_query(comm, con)
                    cents,sz=db.readCentroids(frameID,loadPref)
                    if((data.shape[0],data.shape[1]-1) != sz):
                        print(frameID,[data.shape[0],data.shape[1]-1],data.shape,sz)
                elif(switch=="fps"):
                    comm="select * from cobra_status where mcs_frame_id = """+str(int(frameID))
                    data = sqlio.read_sql_query(comm, con)
                    if(data.shape != (1574,9)):
                        print(frameID,data.shape)
        else:
            if(switch=="visit"):
                comm="select * from pfs_visit where pfs_visit_id = """+str(int(frameIDs[0] // 100))
                data = sqlio.read_sql_query(comm, con)
                if(data.shape != (1,2)):
                    print(frameIDs[0],data.shape)
            if(switch=="phys"):
                comm="select * from mcs_exposure where mcs_frame_id = """+str(int(frameIDs[0]))
                data = sqlio.read_sql_query(comm, con)
                if(data.shape != (1,17)):
                    print(frameIDs[0],data.shape)
            if(switch=="cent"):
                for frameID in frameIDs:
                    comm="select * from mcs_data where mcs_frame_id = """+str(int(frameID))
                    data = sqlio.read_sql_query(comm, con)
                    cents,sz=db.readCentroids(frameID,loadPref)
                    if((data.shape[0],data.shape[1]-1) != sz):
                        print(frameID,[data.shape[0],data.shape[1]-1],data.shape,sz)
            if(switch=="fps"):
                for frameID in frameIDs:
                    comm="select * from cobra_status where mcs_frame_id = """+str(int(frameID))
                    data = sqlio.read_sql_query(comm, con)
                    if(data.shape != (1574,9)):
                        print(frameID,data.shape)



con = psycopg2.connect(user='pfs', host = 'db-ics', database="opdb", password='2394f4s3d', port=5432)
#con = psycopg2.connect(user="karr", host="localhost", database="opdb")
con.autocommit = True

#create a cursor object, by which we interact with the database
cur = con.cursor()

runAll("visit",cur)
runAll("phys",cur)
runAll("cent",cur)
runAll("fps",cur)

