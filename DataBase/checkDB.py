"""   

Routines to check integrity of database population 

Basically, it reads in the input files used to create the database, loads
the values from the database, and compares tthem.

"""


import pandas as pd
import pandas.io.sql as sqlio
import psycopg2
import dbStuff as db
import numpy as np

def checkPhys(con,frameID,loadPref):

    """ telescope data """

    #read from the DB
    
    comm="""select * from mcs_exposure where mcs_frame_id = """+str(int(frameID))
    data = sqlio.read_sql_query(comm, con)

    #read from file (which was used to populate database)
    
    fname=loadPref+"/ldump"+str(frameID)+".npy"
    phys=db.readPhys(frameID,loadPref)

    #different cases for comparing different types - float, int, timestamp and string
    for key in phys:
        if(type(data[key][0]) == np.float64):
            if(not np.isclose(np.float(phys[key]),data[key][0])):
                print("Problem: ",frameID)
        elif(type(data[key][0]) == np.int64):
            if(np.int(phys[key]) != data[key][0]):
                print("Problem: ",frameID,key)
        elif(type(data[key][0]) == pd.Timestamp):
            if(pd.Timestamp(phys[key]) != data[key][0]):
                print("Problem: ",frameID)
        else:
            if(phys[key] != data[key][0]):
                print("Problem: ",frameID,key)
 
def checkVisit(cur,frameID,loadPref):

    pass

def checkCentroids(con,frameID,loadPref):

    """ centroid data """
    
    #read from database
    comm="""select * from mcs_data where mcs_frame_id ="""+str(int(frameID))+""" order by spot_id"""
    data = sqlio.read_sql_query(comm, con)

    #read from file  (which was used to populate database)
    cents,sz=db.readCentroids(frameID,loadPref)

    #compare - here all the entries are numbers
    for key in cents:
        for i in range(sz[0]):
            if(not np.isclose(np.float(cents[key][i]),data[key][i])):
                print(key,i,cents[key][i],data[key][i])
        

def checkCobra(con,frameID,loadPref):

    """ compare matched fibres """

    #read from database
    comm="""select * from cobra_status where mcs_frame_id ="""+str(int(frameID))+""" order by cobra_id"""
    data = sqlio.read_sql_query(comm, con)

    #read from file  (which was used to populate database)    
    fps,sz=db.readFPS(frameID,loadPref)

    #all the entries here are numbers
    for key in fps:
        for i in range(sz[0]):

            #turn None into NaN for comparison
            if(fps[key][i] == None):
                fps[key][i] = np.nan

            #equal_nan so that un matched cobras are compared correctly
            if(not np.isclose(np.float(fps[key][i]),data[key][i],equal_nan=True)):
                print(key,i,fps[key][i],data[key][i])


con = psycopg2.connect(user='pfs', host = 'db-ics', database="opdb", password='2394f4s3d', port=5432)
#con = psycopg2.connect(user="karr", host="localhost", database="opdb")
con.autocommit = True
cur = con.cursor()

#test a couple of sets for each d
frameIDs=[845000,1124800,1131100,1421900,1430100,1645400,1683100,1758500,1760400,1764129]
loadPref="dbSet"

for frameID in frameIDs:
    print(frameID)
    checkPhys(con,frameID,loadPref)
    checkCentroids(con,frameID,loadPref)
    checkCobra(con,frameID,loadPref)
