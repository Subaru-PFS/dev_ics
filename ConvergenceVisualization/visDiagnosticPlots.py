"""

Class to visualize convergence results

"""


import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import glob
import astropy.io.fits as pyfits
import numpy as np
import matplotlib.pyplot as plt
import os
import subprocess
import sys
import math
import pathlib
from matplotlib.colors import to_hex

from bokeh.io import output_notebook, show, export_png,export_svgs
from selenium import webdriver
from webdriver_manager.chrome import ChromeDriverManager


from bokeh.plotting import figure, show, output_file
import bokeh.palettes
from bokeh.layouts import column,gridplot
from bokeh.models import HoverTool, ColumnDataSource, LinearColorMapper
from bokeh.models.glyphs import Text

from bokeh.transform import linear_cmap
from bokeh.palettes import Category20

from ics.cobraCharmer import pfiDesign
from ics.cobraCharmer import func
import fnmatch
from ics.fpsActor import fpsFunction as fpstool
import pandas as pd
from opdb import opdb
import logging

from pfs.utils.butler import Butler
import pfs.utils.coordinates.transform as transformUtils

from ics.cobraCharmer import pfi 
from ics.cobraCharmer import func




#logging.basicConfig(format="%(asctime)s.%(msecs)03d %(levelno)s %(name)-10s %(message)s",
#                    datefmt="%Y-%m-%dT%H:%M:%S")
#logger = logging.getLogger('visDianosticPlot')
#logger.setLevel(logging.INFO)


class VisDiagnosticPlot(object):

    def __init__(self, xmlFile, dotFile, hostname = 'opdb', port= '5432', dbname = 'opdb', username = 'pfs'):


        # initialize convergence data to None
        
        self.convergeData = None
        self.visitId = None
        self.nIter = None
        
        
        #load pfi design file
        des = pfiDesign.PFIDesign(xmlFile)
        self.calibModel = des
        cobras = []

        #get goodIdx 
        for i in des.findAllCobras():
            c = func.Cobra(des.moduleIds[i],
                           des.positionerIds[i])
            cobras.append(c)
        allCobras = np.array(cobras)
        nCobras = len(allCobras)
        
        goodNums = [i+1 for i,c in enumerate(allCobras) if
                des.cobraIsGood(c.cobraNum, c.module)]
        badNums = [e for e in range(1, nCobras+1) if e not in goodNums]


        self.goodIdx = np.array(goodNums, dtype='i4') - 1
        self.badIdx = np.array(badNums, dtype='i4') - 1

        #load dots
        self.dots=pd.read_csv(dotFile)

        #initialize database connection
        self.db=opdb.OpDB(hostname = hostname, port = port,dbname = dbname,
                                username = username)

        #load fiducials
        # Read fiducial and spot geometry
        butler = Butler(configRoot=os.path.join(os.environ["PFS_INSTDATA_DIR"], "data"))

        self.fids = butler.get('fiducials')

        # for referencing cobras in angle measurements
        self.cobraDefs()
        self.pfi = pfi.PFI(doLoadModel=False,
                      doConnect=False)
        
    def _sequentialColour(self):

        """
        Get an evenly spaced set of colours for motion plots, based on number of iterations

        """

        # colours from red -> purple
        
        cmap = plt.get_cmap('rainbow')
        cols = cmap(np.linspace(0, 1, self.nIter+1))
        self.colSeq = [] 

        # turn into hex values for plots
        for i in range(self.nIter):
            self.colSeq.append(to_hex(cols[i]))
            
    def loadConvergence(self,visitId):

        """
        load data to plot the results of a convergence run.
        This does a join on cobra_target and cobra_match to get both target and actual positions.
        This loads the results at a given iteration
        """
        self.visitId = visitId

        # make sql select statuement. 
        sql = f'select cobra_target.pfs_visit_id, cobra_target.iteration, cobra_target.cobra_id, cobra_target.pfi_nominal_x_mm, cobra_target.pfi_nominal_y_mm, cobra_target.pfi_target_x_mm, cobra_target.pfi_target_y_mm, cobra_match.pfi_center_x_mm, cobra_match.pfi_center_y_mm  from cobra_target full join cobra_match  on cobra_target.pfs_visit_id = cobra_match.pfs_visit_id and cobra_target.iteration = cobra_match.iteration  and cobra_target.cobra_id = cobra_match.cobra_id where cobra_match.pfs_visit_id = {self.visitId} '

        # get data
        self.convergeData = self.db.fetch_query(sql)

        # save the maximum iteration, as this is used a lot
        self.nIter = self.convergeData['iteration'].values.max()

        # set the colour sequenced based on nIter
        self._sequentialColour()
        
                           
    def visConvergenceHist(self,range=(0,0.08), bins=20, saveFile = False):

        """
        Plot histograms of distance from the target for a convergence sequence

        """

        if(isinstance(self.convergeData, type(None))):
            print("No data loaded")
            return

        
        fig,ax = plt.subplots()

        ax.set_aspect('auto')

        for iterVal in np.arange(0,self.nIter):

            # filter for the iteration
            filterInd = self.convergeData['iteration'] == iterVal
            cD = self.convergeData[filterInd]

            #calculate distance from targets at this iteration
        
            dist=np.sqrt((cD['pfi_center_x_mm'].values-cD['pfi_target_x_mm'].values)**2+
                         (cD['pfi_center_y_mm'].values-cD['pfi_target_y_mm'].values)**2)

            # the histogram
            n, bins, patches = ax.hist(dist,range=range, bins=bins, alpha=0.7,
                    histtype='step',linewidth=3,
                    label=f'{iterVal+1}-th Iteration')
        plt.legend(loc='upper right')
        ax.set_title("Distance to Target")
        ax.set_xlabel("Distance (mm)")
        ax.set_ylabel("N")

        if(saveFile == True):
            fName = f'{self.visitId:d}_{iterVal:d}.pdf'
            plt.savefig(fName)
        
        return fig,ax

    def visConvergenceBool(self,iterVal=None, threshNon = 0.01, threshBad = 0.08, saveFile = False, **kwargs):
        """ 
        Plot convergence/non convergence from the target position at the specified iteration.

        Three colours are used, one for convergence, one for close to convergence, one for non converged. 
        
        iterVal starts from 0, if set to None will default to maximum value
        """

     
        fig,ax=plt.subplots()

        # check if convergence data has been loaded yet
        
        if(isinstance(self.convergeData, type(None))):
            print("No data loaded")
            return

        # check iterVal, set to max if undefined
        
        if(iterVal == None):
            iterVal = self.nIter

        # filter the dataframe for the iteration value
        
        filterInd = self.convergeData['iteration'] == iterVal
        cD = self.convergeData[filterInd]

        #calculate distance from targets at this iteration
        
        dist=np.sqrt((cD['pfi_center_x_mm'].values-cD['pfi_target_x_mm'].values)**2+
            (cD['pfi_center_y_mm'].values-cD['pfi_target_y_mm'].values)**2)

        ind1 = np.where(dist < threshNon)
        ind2 = np.where(np.all([dist >= threshNon, dist < threshBad],axis=0))
        ind3 = np.where(dist >= threshBad)

        
        # do a scatter plot
        sc=ax.scatter(self.calibModel.centers.real[self.goodIdx[ind1]],self.calibModel.centers.imag[self.goodIdx[ind1]],
                      c='purple',marker='o', s=20, **kwargs)
        sc=ax.scatter(self.calibModel.centers.real[self.goodIdx[ind2]],self.calibModel.centers.imag[self.goodIdx[ind2]],
                      c='green',marker='o', s=20, **kwargs)
        sc=ax.scatter(self.calibModel.centers.real[self.goodIdx[ind3]],self.calibModel.centers.imag[self.goodIdx[ind3]],
                      c='yellow',marker='o', s=20, **kwargs)

        plt.colorbar(sc)
            
        if(saveFile == True):
            fName = f'{self.visitId:d}_{iterVal:d}.pdf'
            plt.savefig(fName)

        return fig,ax
     
    def visConvergenceMap(self,iterVal=None, saveFile = False, **kwargs):

        """ 
        Plot distance from the target position at the specified iteration.
        
        iterVal starts from 0, if set to None will default to maximum value
        """
        
        fig,ax=plt.subplots()

        # check if convergence data has been loaded yet
        
        if(isinstance(self.convergeData, type(None))):
            print("No data loaded")
            return

        # check iterVal, set to max if undefined
        
        if(iterVal == None):
            iterVal = self.nIter

        # filter the dataframe for the iteration value
        
        filterInd = self.convergeData['iteration'] == iterVal
        cD = self.convergeData[filterInd]

        #calculate distance from targets at this iteration
        
        dist=np.sqrt((cD['pfi_center_x_mm'].values-cD['pfi_target_x_mm'].values)**2+
            (cD['pfi_center_y_mm'].values-cD['pfi_target_y_mm'].values)**2)


        # do a scatter plot
        sc=ax.scatter(self.calibModel.centers.real[self.goodIdx],self.calibModel.centers.imag[self.goodIdx],
                      c=dist,marker='o', s=20, **kwargs)

        plt.colorbar(sc)

        # some labels
        ax.set_xlabel("X (mm)")
        ax.set_ylabel("Y (mm)")
    
        tString = f'pfsVisitId = {self.visitId:d}, iteration = {iterVal:d}'
        ax.set_title(tString)

        ax.set_aspect('equal')

        if(saveFile == True):
            fName = f'{self.visitId:d}_{iterVal:d}.pdf'
            plt.savefig(fName)

        
        return fig,ax

    def cobraDefs(self):

        """
        Get cobra definitions for angle calculations
        """
        
        cobras=[]
        for i in self.calibModel.findAllCobras():
            c = func.Cobra(self.calibModel.moduleIds[i],
                           self.calibModel.positionerIds[i])
            cobras.append(c)
            self.allCobras = np.array(cobras)

    def singleCobraMotion(self,cobraNum, saveFile = False, **kwargs):
        """
        diagnostic plots for a single cobra.
        """
        
        fig,axes=plt.subplots(2,2)

        # check for loaded data
        if(isinstance(self.convergeData, type(None))):
            print("No data loaded")
            return

        if(cobraNum-1 not in self.goodIdx):
            print("Not a valid cobra number")
            return 
        
        # extract motion of a single cobra; target and measured

        cobraInd = cobraNum - 1
        filterInd = self.convergeData['cobra_id'] == cobraNum
        cM = self.convergeData[filterInd]

        xM = cM['pfi_center_x_mm'].values
        yM = cM['pfi_center_y_mm'].values
        xT = cM['pfi_target_x_mm'].values
        yT = cM['pfi_target_y_mm'].values


        # turn (x,y) into (theta,phi)
        
        pM, tM = self.getAnglesOld(xM, yM, cobraNum)
        pT, tT = self.getAnglesOld(xT, yT, cobraNum)
        dR=np.sqrt((xT-xM)**2+(yT-yM)**2)
    
        xC = self.calibModel.centers[cobraInd].real
        yC = self.calibModel.centers[cobraInd].imag
        aL = self.calibModel.L1[cobraInd] + self.calibModel.L2[cobraInd]
        
        self.movPlot(axes[0,0],xM,yM,xT,yT,xC,yC,aL)
        self.convPlot(axes[0,1],"t",tM-tT)
        self.convPlot(axes[1,0],"p",pM-pT)
        self.convPlot(axes[1,1],'r',dR)

        tstring = f'Cobra = {cobraNum:d}; VisitId = {self.visitId:d}'
        print(tstring)
        fig.suptitle(tstring)
        plt.tight_layout()

        if(saveFile == True):
            fName = f'{self.visitId:d}_{cobraNum:d}.pdf'
            plt.savefig(fName)

        
        return fig,ax

    def getAnglesOld(self, xM, yM, cNum):

        positions = xM+yM*1j
        cIdx=cNum-1
        relativePositions = positions - self.calibModel.centers[cIdx]
        distance = np.abs(relativePositions)
        L1 = self.calibModel.L1[cIdx]
        L2 = self.calibModel.L2[cIdx]
        distanceSq = distance ** 2
        L1Sq = L1 ** 2
        L2Sq = L2 ** 2
        phiIn = self.calibModel.phiIn[cIdx] + np.pi
        phiOut = self.calibModel.phiOut[cIdx] + np.pi
        tht0 = self.calibModel.tht0[cIdx]
        tht1 = self.calibModel.tht1[cIdx]
        phi = np.full((len(xM), 2), np.nan)
        tht = np.full((len(xM), 2), np.nan)

        for i in range(len(xM)):
            ang1 = np.arccos((L1Sq + L2Sq - distanceSq[i]) / (2 * L1 * L2))
            ang2 = np.arccos((L1Sq + distanceSq[i] - L2Sq) / (2 * L1 * distance[i]))
            
            phi[i][0] = ang1 - phiIn
            tht[i][0] = (np.angle(relativePositions[i]) + ang2 - tht0) % (2 * np.pi)

            # check if there are additional solutions
            if ang1 <= np.pi/2 and ang1 > 0:
                # phiIn < 0
                phi[i][1] = -ang1 - phiIn
                tht[i][1] = (np.angle(relativePositions[i]) - ang2 - tht0) % (2 * np.pi)
                # check if tht is within two theta hard stops
            elif ang1 > np.pi/2 and ang1 < np.pi:
                phi[i][1] = 2 * np.pi - ang1 - phiIn
                tht[i][1] = (np.angle(relativePositions[i]) - ang2 - tht0) % (2 * np.pi)

        return phi[:,0],tht[:,0]


    def sequencePlot(self, centrePos = True, hardStop = False, blackDots = False, badCob = True, patrolRegion = True, ff = True, saveFile = False, **kwargs):

        """
        Plot a sequence of moves over the whole array.

        The iterations go in sequence from red -> purple 

        """

        
        fig,ax = plt.subplots()

        # plot spots for each iteration in a given colour
        for iterVal in range(self.nIter):
            filterInd = self.convergeData['iteration'] == iterVal
            cD = self.convergeData[filterInd]

            xM = cD['pfi_center_x_mm'].values
            yM = cD['pfi_center_y_mm'].values

            sc = ax.scatter(xM, yM, c=self.colSeq[iterVal], marker = 'o', s=10, **kwargs)
        ax.set_aspect('equal')
        ax.set_xlabel("X (mm)")
        ax.set_ylabel("Y (mm)")

        tString = f'pfsVisitId = {self.visitId:d}'
        ax.set_title(tString)
                
        # various optional overlays
        if(centrePos == True):
            ax = self.overlayCentres(ax)
        if(hardStop == True):
            ax = self.overlayHardStop(ax)
        if(blackDots == True):
            ax = self.overlayBlackDots(ax)
        if(badCob == True):
            ax = self.overlayBadCob(ax)
        if(patrolRegion == True):
            ax = self.overlayPatrolRegion(ax)
        if(patrolRegion == True):
            ax = self.overlayPatrolRegion(ax)
        if(ff == True):
            ax =self.overlayFF(ax)

        if(saveFile == True):
            fName = f'{self.visitId:d}.pdf'
            plt.savefig(fName)
            
        return fig,ax

    def overlayFF(self, ax):

        """
        overlay positions of fiducial fibres

        """
                
        ax.scatter(self.fids['x_mm'],self.fids['y_mm'],marker="d")

        return ax
                
    def overlayPatrolRegion(self, ax, cobraNum = None):

        """
        overlay cobra patrol regions on given axis.

        if cobraNum == None do for all good cobras, else for cobraNum

        """

        if(cobraNum == None):
            ind = self.goodIdx
        else:
            ind = [cobraNum - 1]

        armLength = self.calibModel.L1 + self.calibModel.L2

        
        for i in ind:
        
            circle=plt.Circle((self.calibModel.centers[i].real, self.calibModel.centers[i].imag),armLength[i],fill=False,color='black')
            a=ax.add_artist(circle)

        return ax
    
    def overlayCentres(self,ax,cobraNum = None):

        """
        overlay centre positions of cobra patrol regions on given axis.

        if cobraNum == None do for all good cobras, else for cobraNum

        """

        if(cobraNum == None):
            ind = self.goodIdx
        else:
            ind = [cobraNum - 1]

        ax.scatter(self.calibModel.centers[ind].real, self.calibModel.centers[ind].imag, c='black', marker = 'o', s=25)

        return ax
    
    def overlayBadCob(self, ax,cobraNum= None):

        """
        overlay centre positions of cobra patrol regions on given axis.

        if cobraNum == None do for all good cobras, else for cobraNum

        """
        
        if(cobraNum == None):
            ind = self.goodIdx
        else:
            ind = [cobraNum - 1]

        ax.scatter(self.calibModel.centers[ind].real, self.calibModel.centers[ind].imag, c='black', marker = '*', s=25)

        return ax
    
    def overlayHardStop(self,ax,cobraNum= None):

        """
        overlay theta hard stops

        if cobraNum == None do for all good cobras, else for cobraNum
        """
        

        if(cobraNum == None):
            ind = self.goodIdx
        else:
            ind = [cobraNum - 1]
        
        length=self.calibModel.L1+self.calibModel.L1 

        self._addLine(ax,self.calibModel.centers[ind],length[ind],self.calibModel.tht0[ind],color='orange',
                            linewidth=0.5,linestyle='--')

        self._addLine(ax,self.calibModel.centers[ind],length[ind],self.calibModel.tht1[ind],color='black',
                            linewidth=0.5,linestyle='-.')

        return ax

    def _addLine(self, ax, centers, length, angle, cobraNum= None, **kwargs):

        x = length*np.cos(angle)
        y = length*np.sin(angle)
        
        for i in range(len(centers)):
                       
            ax.plot([centers.real[i],  centers.real[i]+x[i]],
                        [centers.imag[i],centers.imag[i]+y[i]],**kwargs)
        
        return ax
    def overlayBlackDots(self,ax,cobraNum= None):

        """
        overlay black dots

        if cobraNum == None do for all good cobras, else for cobraNum
        """

        if(cobraNum == None):
            ind = self.goodIdx
        else:
            ind = [cobraNum - 1]

    
        for i in ind:
            e = plt.Circle((self.dots['x'].values[i], self.dots['y'].values[i]), self.dots['r'].values[i], 
                        color='grey', fill=True, alpha=0.5)
            ax.add_artist(e)
     
        return ax
        
    def movPlot(self,ax,xM,yM,xT,yT,xC,yC,aL):
    
        """
        plot the 2D motion of a single cobra over a convergence run. Generally called by badCobraDiagram.

        input: 
        ax: axis on which to plot the diagram
        xC, yC, aL: centres and armlengths for cobra patrol regions
        xM, yM: measured positions (same units as above)
        xT, yT: target positions (same units as above)
    
        output: plot to the provided axis

        """

        #sequence of colours for spots, goes basically red -> purple in chromatic order
    
        
        #plot the motions in sequence
        for i in range(self.nIter):
            ax.scatter(xM[i],yM[i],color=self.colSeq[i])

        #draw patrol region and center
        circle=plt.Circle((xC,yC),aL,fill=False,color='black')
        a=ax.add_artist(circle)
        a=ax.scatter(xC,yC,color='black')

        #target - black adn white x so it shows over background and spots
        a=ax.scatter(xT,yT,c='black',marker="+",s=80)
        a=ax.scatter(xT,yT,c='white',marker="x",s=80)
        
        #adjust limits
        a=ax.set_xlim((xC-aL*1.3,xC+aL*1.3))
        a=ax.set_ylim((yC-aL*1.3,yC+aL*1.3))
        a=ax.set_aspect('equal')

    def convPlot(self,ax,var,diff):
        
        """
        plot linear plots of delta(Values) for convergence.
        
        Inputs
          ax: axis to plot on
          var: variable to plot (p, t or r for phi, theta or distance from target)
          diff: detlta value (measured - target)
        
        """

        # plot delta vs iteration, use same colours as movePlot
        ax.plot(np.arange(1,self.nIter+2),diff,color='black')
        for i in range(0,self.nIter):
            ax.plot([i+1],diff[i],color=self.colSeq[i],marker="d")
            
            # a line at zero
            ax.axhline(y=0,color='black',linestyle='--')
            
            #labels
            ax.set_xlabel("Iteration")
            if(var=="p"):
                ax.set_ylabel("d(Phi)")
            if(var=="t"):
                ax.set_ylabel("d(Theta)")
            if(var=="r"):
                ax.set_ylabel("d(R)")        
        
def getAngles(self,xM,yM,cNum):

    """
    Replace this with proper code from cobraCharmer
    """
    

    positions = xM+yM*1j
    cIdx=cNum-1
    relativePositions = positions - self.calibModel.centers[cIdx]
    distance = np.abs(relativePositions)
    L1 = self.calibModel.L1[cIdx]
    L2 = self.calibModel.L2[cIdx]
    distanceSq = distance ** 2
    L1Sq = L1 ** 2
    L2Sq = L2 ** 2
    phiIn = self.calibModel.phiIn[cIdx] + np.pi
    phiOut = self.calibModel.phiOut[cIdx] + np.pi
    tht0 = self.calibModel.tht0[cIdx]
    tht1 = self.calibModel.tht1[cIdx]
    phi = np.full((len(xM), 2), np.nan)
    tht = np.full((len(xM), 2), np.nan)

    for i in range(len(xM)):
        ang1 = np.arccos((L1Sq + L2Sq - distanceSq[i]) / (2 * L1 * L2))
        ang2 = np.arccos((L1Sq + distanceSq[i] - L2Sq) / (2 * L1 * distance[i]))

        phi[i][0] = ang1 - phiIn
        tht[i][0] = (np.angle(relativePositions[i]) + ang2 - tht0) % (2 * np.pi)

        # check if there are additional solutions
        if ang1 <= np.pi/2 and ang1 > 0:
            # phiIn < 0
            phi[i][1] = -ang1 - phiIn
            tht[i][1] = (np.angle(relativePositions[i]) - ang2 - tht0) % (2 * np.pi)
            # check if tht is within two theta hard stops
        elif ang1 > np.pi/2 and ang1 < np.pi:
            phi[i][1] = 2 * np.pi - ang1 - phiIn
            tht[i][1] = (np.angle(relativePositions[i]) - ang2 - tht0) % (2 * np.pi)

    return phi[:,0],tht[:,0]
        
