"""
routines to calculate statistics regarding non-converging cobras
"""

from opdb import opdb
import numpy as np

def readConvergenceRun(db, visitId):

    """
    pull the data for a single convergence run from the database, 
    from the pfs_visit_id, ordered by cobra id and then iteration

    Input:
      db: database connection
      visitId: pfs_visit_id value

    Returns; 
      df: results in pandas dataframe

    """

    sql = f'SELECT cobra_match.pfs_visit_id, cobra_match.iteration, cobra_match.cobra_id, cobra_match.pfi_center_x_mm, cobra_match.pfi_center_y_mm, cobra_target.pfi_target_x_mm, cobra_target.pfi_target_y_mm FROM cobra_match JOIN cobra_target ON (cobra_match.pfs_visit_id = cobra_target.pfs_visit_id AND cobra_match.iteration = cobra_target.iteration AND cobra_match.cobra_id = cobra_target.cobra_id) WHERE cobra_match.pfs_visit_id = {visitId} ORDER BY cobra_match.cobra_id, cobra_match.pfs_visit_id, cobra_match.iteration'

    df = db.fetch_query(sql)
    return df

def convStats(df, tol, nomIt):

    """
    for a single convergence run, given the tolerance value used to calculate convergence, 
    find the iteration of convergence for all cobras. Note that tol must match the value
    used in the convergence run. 

    Input:
      df: dataframe returned by readConvergenceRun
      tol: distance from target (in mm) regraded as converged
      nomIt: iteration which is considered a successful convergence

    returns:
      convIt: a 2394 element array, set to -1 for broken cobras, iteration of convergence, or maxIt+1 for non converged
      convSuc: a 2394 element array, set to -1 for broken, 0 for nonconverged and 1 for converged, at nomIt

    """

    convIt = np.zeros((2394))

    # extract iteration, cobra numbers, and distance from target in mm
    iNum = df['iteration'].values
    cNum = df['cobra_id'].values
    dd = np.sqrt((df['pfi_center_x_mm'].values-df['pfi_target_x_mm'].values)**2+(df['pfi_center_y_mm'].values-df['pfi_target_y_mm'].values)**2)

    
    # get number of iterations. This method will work regardless of 0 or 1 start
    nIter = len(np.unique(iNum))

    # reshape the arrays
    iN = iNum.reshape((int(len(iNum)/nIter), nIter))
    cN = cNum.reshape((int(len(iNum)/nIter), nIter))
    dist = dd.reshape((int(len(iNum)/nIter), nIter))
    
    # with the reshaped arrays, find the first iteration where the cobra has converged
    convIt = np.argmax(dist < tol, axis = 1)
    
    # set unconverged values to nIter+1
    ind = np.where(convIt == 0)
    convIt[ind] = nIter+1

    # put into a 2394 length array
    
    # the values will be -1 if the cobra/fibre is broken, N for the convergence iteration, or nIter+1 if the 
    # cobra doesn't converge
    
    convIt = np.zeros((2394))-1
    convIt[cN[:, 0].ravel()-1] = convIt
    
    ind = np.where(np.all([convIt > 0, convIt <= nomIt]))
    ind1 = np.where(convIt < 0)
    
    convSuc = np.zeros((2394))
    convSuc[ind] = 1
    convSuc[ind1] = -1


    return convIt, convSuc

def tallyConvegence(db, visitIds, nomIts, tols):

    """
    keep a running tally of cobra success stats for a series of convergence runs

    input:
      db: database connection
      visitIds: list of visitIds
      nomIts: list of iteration of convergence each run (same length as visitIds)
      tols: list of distance for convergence for each run  (same length as visitIds)

    returns:
      failStat: number of failed convergence runs for each cobra
      nRun: number of runs where the cobra wasn't considered broken
      nRuns: total number of cobras

      divide failStat / nRun for non zero nRun to get fraction of runs that failed.

    """

    failStat = np.zeros((2394))
    nRun = np.zeros((2394))
    nRuns = len(visitIds)

    
    # cycle through the list of visits
    for visitIt, nomit, tol in zip(visitIds, nomIts, tols):

        #get values
        db = readConvergenceRun(db, visitId)
        #calculate convergence stats for single run
        convIt, convSuc = convStats(df, tol, nomIt)

        #increment nRun if cobra not broken. 
        ind = np.where(convSuc >= 0):
        nRun[ind] += 1
        #increment failState if failed
        ind = np.where(convSuc == 0):
        failStat[ind] += 1

        
    return failStat, nRun, nRuns
        


