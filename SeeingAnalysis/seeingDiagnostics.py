"""
produce seeing/image quality plots with database retrieval
"""
import numpy as np
import matplotlib.pylab as plt

def getSeeingRunFromDB(db,mcsFrameIds):

    sql=f'select * from mcs_data full join cobra_match on mcs_data.mcs_frame_id=cobra_match.mcs_frame_id and mcs_data.spot_id = cobra_match.spot_id where cobra_match.mcs_frame_id between {mcsFrameIds[0]} and {mcsFrameIds[1]} order by cobra_match.mcs_frame_id , cobra_match.cobra_id'

    df = db.fetch_query(sql)

    return df

def calcSeeing(df):

    #extract the variables
    fx=df['mcs_second_moment_x_pix'].values
    fy=df['mcs_second_moment_y_pix'].values
    x=df['pfi_center_x_mm'].values
    y=df['pfi_center_y_mm'].values
    p=df['peakvalue'].values
    b=df['bgvalue'].values
    itNum=df['iteration'].values
    cobNum=df['cobra_id'].values
    itMax=df['iteration'].values.max()
    nIter=itMax+1

    #this massages them into a numpy array of the right form
    xArray=np.zeros((2394,nIter))
    yArray=np.zeros((2394,nIter))
    fxArray=np.zeros((2394,nIter))
    fyArray=np.zeros((2394,nIter))
    pArray=np.zeros((2394,nIter))
    cArray=np.zeros((2394,nIter))
    
    
    for i in range(1,nIter+1):
        ind=np.where(itNum==i)
        xArray[:,i]=x[ind].ravel()
        yArray[:,i]=y[ind].ravel()
        fxArray[:,i]=fx[ind].ravel()
        fyArray[:,i]=fy[ind].ravel()
        pArray[:,i]=p[ind].ravel()
        bArray[:,i]=b[ind].ravel()
       
        xArray=ma.masked_where(xArray <= 0 ,xArray)
        yArray=ma.masked_where(xArray <= 0 ,yArray)
        fxArray=ma.masked_where(xArray <= 0 ,fxArray)
        fyArray=ma.masked_where(xArray <= 0 ,fyArray)
        backArray=ma.masked_where(xArray <= 0 ,backArray)
        peakArray=ma.masked_where(xArray <= 0 ,peakArray)
        qualArray=ma.masked_where(xArray <= 0 ,qualArray)

        dd = np.zeros(xArray.shape)
        xd = np.zeros(xArray.shape)
        yd = np.zeros(xArray.shape)

        xFirst=xArray[:,0]
        yFirst=yArray[:,0]
        
        for i in range(xArray.shape[0]):
            dd[:,i]=np.sqrt((xArray1[:,i]-xfirst[i])**2+(yArray1[i,:]-yfirst[i])**2)
            xd[:,i]=np.sqrt((xArray1[:,i]-xfirst[i])**2)
            yd[:,i]=np.sqrt((yArray1[:,i]-yfirst[i])**2)

        #adjust the masks if needed
    dd=ma.masked_where(((dd <=0) | (xArray.mask == True)), dd)
    xd=ma.masked_where(((dd <=0) | (xArray.mask == True)), xd)
    yd=ma.masked_where(((dd <=0) | (xArray.mask == True)), yd)

    #get rms of the values
    nzero=np.count_nonzero(xArray,axis=1)
    rmsVal=dd.std(axis=1)
    rmsX=xd.std(axis=1)
    rmsY=yd.std(axis=1)
    xAv,yAv,fxAv,fyAv,peakAv = getAverages(frameIDs,xArray, yArray, fxArray, fyArray, peakArray)

    return xAv,yAv,fxAv,fyAv,rmsVal,rmsX,rmsY,nzero,peakAv
    
def getAverages(frameIDs,xArray, yArray, fxArray, fyArray, peakArray):

    xAv=xArray.mean(axis=1)
    yAv=yArray.mean(axis=1)
    fxAv=fxArray.mean(axis=1)
    fyAv=fyArray.mean(axis=1)
    peakAv=peakArray.mean(axis=1)

    return xAv,yAv,fxAv,fyAv,peakAv


def pairPlot(xs,ys,val1,val2,plotrange,titl,outPref, suffix, valtype,   units,   nbins,stitle):

    """

    plot a pair of plots, one showing the map, the other the histogram.

    Input: 
       xs,ys - coordinates of points
       val1 - value for the map, same dimensions as xs,ys
       val2 - value for the histogram
       plotrange - [min,max] for polotting range on both plots
       titl - title
       outPref - outPref for output files
       suffix - suffix for output files
       valtype - type of variable plotted (e.g., FWHM (x))
       units - unites (eg pixels)
       nbins - number of bins for histogram
       stitle - optional subbbtitle.

    """

    #set up plot
    fig,axes=plt.subplots(1,2)
    fig.set_figheight(4)
    fig.set_figwidth(10)

    if(plotRange==None):

        mm=val1.mean()
        st=val1.std()
        plotRange=[mm-2*st,mm+2*st]


    #calculate the bins
    binsize=(plotRange[1]-plotRange[0])/nbins
    bins=np.arange(plotRange[0],plotRange[1]+binsize,binsize)

    #histogram - compressed deals correctly with masked values
    hi=axes[0].hist(val2.compressed(),bins=bins)
    #labels
    axes[0].set_xlabel(valtype)
    axes[0].set_ylabel("N")
    #scatter plot. size of points optimized for file/notebooks
    sc=axes[1].scatter(xs,ys,c=val1,marker="o",cmap='Purples',lw=0,s=20,vmin=plotRange[0],vmax=plotRange[1])
    fig.colorbar(sc)
    #label the axes
    axes[1].set_xlabel("X ("+units+")")
    axes[1].set_ylabel("Y ("+units+")")

    plt.suptitle(valtype+stitle)

    if(outFile != None):
        plt.savefig(outPref+suffix+".png")
 

def runSeeing(db,mcsFrameIds,stitle,plotRangeRMS=None,plotRangeX=None,plotRangeY=None,plotRangeP=None,outFile=None):

    """
    call the sequence
    """
    
    df=getSeeingRunFromDB(db,mcsFrameIds)
    
    xAv,yAv,fxAv,fyAv,rmsVal,rmsX,rmsY,nzero,peakAv=calcSeeing(df)

    nbins=30
    dr.pairPlot(xAv.ravel(),yAv.ravel(),rmsVal.ravel(),rmsVal.ravel(),plotRangeRMS,"RMS",outPref,"_rms"+source,"RMS","mm",nbins,stitle)
    dr.pairPlot(xAv.ravel(),yAv.ravel(),rmsVal.ravel(),rmsVal.ravel(),plotRangeRMS,"Spot Size (x)",outPref,"_fx"+source,"Spot Size (x)","mm",nbins,stitle)
    dr.pairPlot(xAv.ravel(),yAv.ravel(),rmsVal.ravel(),rmsVal.ravel(),plotRangeRMS,"Spot Size (y)",outPref,"_fy"+source,"Spot Size (y)","mm",nbins,stitle)
    dr.pairPlot(xAv.ravel(),yAv.ravel(),rmsVal.ravel(),rmsVal.ravel(),plotRangeRMS,"Peak Value",outPref,"_peak"+source,"Peak Value","mm",nbins,stitle)

