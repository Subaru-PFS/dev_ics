
import numpy as np
import numpy.ma as ma
from scipy.stats import sigmaclip
import matplotlib.pylab as plt
import visRoutines as vis


def scatterPlot(xs,    ys,    val1,plotrange,titl,      prefix, suffix, valtype,   units,inter,rx,ry,swap,stitle=""):

    #set up plot
    fig,axes=plt.subplots()

    if(swap==1):
        xx=ys
        yy=xs
    else:
        xx=xs
        yy=ys
    
    #scatter plot. size of points optimized for file/notebooks
    sc=axes.scatter(xx,yy,c=val1,marker="o",cmap='Purples',lw=0,vmin=plotrange[0],vmax=plotrange[1])
    fig.colorbar(sc)
    #label the axes
    axes.set_xlabel("X ("+units+")")
    axes.set_ylabel("Y ("+units+")")

    plt.suptitle(valtype)

    ylim=axes.get_ylim()
    xlim=axes.get_xlim()

    if(rx==1):
        axes.set_xlim((xlim[1],xlim[0]))
    if(ry==1):
        axes.set_ylim((ylim[1],ylim[0]))
    axes.axis('equal')

    #show and save
    if(inter == 1):
        plt.show()
    plt.savefig(prefix+suffix+".png")


def pairPlotRev(xs,    ys,    val1,val2,           plotrange,titl,      prefix, suffix, valtype,   units,   nbins,inter,rx,ry,stitle=""):

    """

    plot a pair of plots, one showing the map, the other the histogram.

    Input: 
       xs,ys - coordinates of points
       val1 - value for the map, same dimensions as xs,ys
       val2 - value for the histogram
       plotrange - [min,max] for polotting range on both plots
       titl - title
       prefix - prefix for output files
       suffix - suffix for output files
       valtype - type of variable plotted (e.g., FWHM (x))
       units - unites (eg pixels)
       nbins - number of bins for histogram
       inter - interactive flag
       stitle - optional subbbtitle.

    REturns

       plots to png files, if inter=1 to screen

    """

    #set up plot
    fig,axes=plt.subplots(1,2)
    fig.set_figheight(4)
    fig.set_figwidth(10)

    #scatter plot. size of points optimized for file/notebooks
    sc=axes[0].scatter(xs,ys,c=val1,marker="o",cmap='Purples',lw=0,s=20,vmin=plotrange[0],vmax=plotrange[1])
    #label the axes
    axes[0].set_xlabel("X ("+units+")")
    axes[0].set_ylabel("Y ("+units+")")

    #calculate the bins
    binsize=(plotrange[1]-plotrange[0])/nbins
    bins=np.arange(plotrange[0],plotrange[1]+binsize,binsize)

    #histogram
    hi=axes[1].hist(val2,bins=bins)

    #labels
    axes[1].set_xlabel(valtype)
    axes[1].set_ylabel("N")

    plt.suptitle(valtype)

    ylim=axes[0].get_ylim()
    xlim=axes[0].get_xlim()

    if(rx==1):
        axes[0].set_xlim((xlim[1],xlim[0]))
    if(ry==1):
        axes[0].set_ylim((ylim[1],ylim[0]))

    
    #show and save
    if(inter == 1):
        plt.show()
    plt.savefig(prefix+suffix+".png")

    
def pairPlotCB(xs,    ys,    val1,val2,           plotrange,titl,      prefix, suffix, valtype,   units,   nbins,inter,stitle=""):

    """

    plot a pair of plots, one showing the map, the other the histogram.

    Input: 
       xs,ys - coordinates of points
       val1 - value for the map, same dimensions as xs,ys
       val2 - value for the histogram
       plotrange - [min,max] for polotting range on both plots
       titl - title
       prefix - prefix for output files
       suffix - suffix for output files
       valtype - type of variable plotted (e.g., FWHM (x))
       units - unites (eg pixels)
       nbins - number of bins for histogram
       inter - interactive flag
       stitle - optional subbbtitle.

    REturns

       plots to png files, if inter=1 to screen

    """

    #set up plot
    fig,axes=plt.subplots(1,2)
    fig.set_figheight(4)
    fig.set_figwidth(10)

    #scatter plot. size of points optimized for file/notebooks
    sc=axes[0].scatter(xs,ys,c=val1,marker="o",cmap='Purples',lw=0,s=20,vmin=plotrange[0],vmax=plotrange[1])
    #label the axes
    axes[0].set_xlabel("X ("+units+")")
    axes[0].set_ylabel("Y ("+units+")")
    fig.colorbar(sc)

    #calculate the bins
    binsize=(plotrange[1]-plotrange[0])/nbins
    bins=np.arange(plotrange[0],plotrange[1]+binsize,binsize)

    #histogram
    hi=axes[1].hist(val2,bins=bins)

    #labels
    axes[1].set_xlabel(valtype)
    axes[1].set_ylabel("N")

    plt.suptitle(valtype)
    
    #show and save
    if(inter == 1):
        plt.show()
    plt.savefig(prefix+suffix+".png")
     
def pairPlot(xs,    ys,    val1,val2,           plotrange,titl,      prefix, suffix, valtype,   units,   nbins,inter,stitle=""):

    """

    plot a pair of plots, one showing the map, the other the histogram.

    Input: 
       xs,ys - coordinates of points
       val1 - value for the map, same dimensions as xs,ys
       val2 - value for the histogram
       plotrange - [min,max] for polotting range on both plots
       titl - title
       prefix - prefix for output files
       suffix - suffix for output files
       valtype - type of variable plotted (e.g., FWHM (x))
       units - unites (eg pixels)
       nbins - number of bins for histogram
       inter - interactive flag
       stitle - optional subbbtitle.

    REturns

       plots to png files, if inter=1 to screen

    """

    #set up plot
    fig,axes=plt.subplots(1,2)
    fig.set_figheight(4)
    fig.set_figwidth(10)

    if(plotrange==None):

        mm=val1.mean()
        st=val1.std()
        plotrange=[mm-2*st,mm+2*st]

    #scatter plot. size of points optimized for file/notebooks
    sc=axes[0].scatter(xs,ys,c=val1,marker="o",cmap='Purples',lw=0,s=20,vmin=plotrange[0],vmax=plotrange[1])
    #label the axes
    axes[0].set_xlabel("X ("+units+")")
    axes[0].set_ylabel("Y ("+units+")")

    #calculate the bins
    binsize=(plotrange[1]-plotrange[0])/nbins
    bins=np.arange(plotrange[0],plotrange[1]+binsize,binsize)

    #histogram - compressed deals correctly with masked values
    hi=axes[1].hist(val2.compressed(),bins=bins)

    #labels
    axes[1].set_xlabel(valtype)
    axes[1].set_ylabel("N")

    plt.suptitle(valtype+stitle)
    
    #show and save
    if(inter == 1):
        plt.show()
    plt.savefig(prefix+suffix+".png")
 

def checkCentroids(xc,yc,cutrange,prefix,inter):

    """

    Quick plot of centroids to check results

    input

    xc,yc: centroid coordinates
    cutrange: limit the region of plotting if needed (for bad data)
    prefix: prefix for plots

    returns: plot, to screen and file

    """

    fig,ax = plt.subplots()

    #scatter plot
    
    ax.scatter(xc,yc)
    ax.set_aspect('equal')
    
    #display and save
    if(inter == 1):
        plt.show()
    plt.savefig(prefix+"_checkpoints.png")

def checkMatched(xx,yy,xs,ys,prefix,inter):

    """

    quick plotting routine for measured centroids and pinhole coordinates

    input: 

    xx,yy mask coordiantes
    xs,ys: spot coordinates
    prefix: prefix for image files

    """

    
    fig,ax = plt.subplots()
    
    #scatter plot: centroids in circles, mask in red dots
    
    ax.scatter(xs,ys)
    ax.scatter(xx,yy,s=20,color='r')

    #save and show
    plt.savefig(prefix+"_checkpoints1.png")
    if(inter == 1):
        plt.show()

def plotVal(xs,ys,val,limit,plotrange,titl,prefix,suffix,units,inter,stitle=""):

    """
    
    routine for scatter plot of a variable

    input

    xs,ys: coordinates

    val: variable to plot
    limit: upper limit to filter out of plots
    plotrange: colour range in [vmin,vmax] format, or None for default
    title: title for plot
    prefix: prefix for output files
    suffix: suffix for output files

    """

    
    #a quick kludge to filter out bad quality points in poorly focussed images
    if(limit > 0):
        ind=np.where((val < limit) & (val > 0) & (val > 0))
    else:
        ind=np.arange(len(val))

    #scatter plot, with or without ragne limit
    
    fig, axes = plt.subplots()
    
    if(plotrange != None):
        sc=axes.scatter(xs[ind],ys[ind],c=val[ind],marker="o",cmap='Purples',lw=0,s=20,vmin=plotrange[0],vmax=plotrange[1])
    
    else:
        sc=axes.scatter(xs[ind],ys[ind],c=val[ind],marker="o",cmap='Purples',lw=0,s=20)

    fig.colorbar(sc)
    plt.title(titl+stitle)
    plt.xlabel("X ("+units+")")
    plt.ylabel("Y ("+units+")")
    if(inter == 1):
        plt.show()
    plt.savefig(prefix+suffix+".png")
 
def checkPlots(files,inter,stitle=""):

    """

    TAkes a list of files and plots the mean and rms of the pixels by frame. 

    input: list of files

    output: plots

    """

    
    #text file for output
    nfiles=len(files)
    

    #set up variables
    av=[]
    rms=[]
    frame=np.arange(nfiles)+1

    #cycle through files
    for file in files:
        print(file)
        #read in image
        image=getImage(file)

        #calculate Stats
        rms.append(image.std())
        av.append(image.mean())

    #plot
    fig, (ax1, ax2) = plt.subplots(2,1)
    ax1.plot(frame,av,marker='d',linestyle="-")
    ax2.plot(frame,rms,marker='d',linestyle="-")
    ax1.set_title("Frame Mean"+stitle)
    ax2.set_title("Frame RMS"+stitle)
    ax2.set_xlabel("Frame")
    ax1.set_ylabel("Mean")
    ax2.set_ylabel("RMS")
    plt.tight_layout()
    if(inter == 1):
        plt.show()

def plotTransByFrame(fxFrameAv,fyFrameAv,peakFrameAv,sxAll,syAll,xdAll,ydAll,rotAll,prefix,inter,stitle=""):

    """

    plot the calculated transformations values and averages by frame number.
    takes output generated by getTransByFrame

    input:
    fxFrameAv,fyFrameAv: average FWHM by frame
    peakFrameAv: peak value by frame

    sxAll,syAll: scale in x and y direction
    xdAll,ydAll: trnaslation in x and y
    rotAll: rotation 

    output: plots

    """

    
    #get number of frames
    frames=np.arange(len(fxFrameAv))

    #first set - fwhms and translation (most useful)
    fig,axes=plt.subplots(1,2)
    fig.set_figheight(4)
    fig.set_figwidth(10)

    axes[0].plot(frames,fxFrameAv,marker='d',linestyle="-",color="#1f77b4")
    axes[0].plot(frames,fyFrameAv,marker='s',linestyle="-",color="#ff7f0e")
    axes[0].set_title("FWHM Average by Frame"+stitle)
    axes[0].set_xlabel("Frame #")
    axes[0].set_ylabel("FHWM (pixels)")

    axes[1].plot(frames,xdAll-xdAll.mean(),marker='d',linestyle="-",color="#1f77b4")
    axes[1].plot(frames,ydAll-ydAll.mean(),marker='s',linestyle="-",color="#ff7f0e")
    axes[1].set_title("Translation Average by Frame"+stitle)
    axes[1].set_xlabel("Frame #")
    axes[1].set_ylabel("Translation (pixels)")
    
    plt.savefig(prefix+"_byframe1.png")
    if(inter == 1):
        fig.show()

    #second set - peaks and bakcgrounds
    fig,axes=plt.subplots(1,2)
    fig.set_figheight(4)
    fig.set_figwidth(10)
        
    axes[0].plot(frames,peakFrameAv,marker='d',linestyle="-")
    axes[0].set_title("Peak Average by Frame"+stitle)
    axes[0].set_xlabel("Frame #")
    axes[0].set_ylabel("Peak")

    axes[1].plot(frames,peakFrameAv,marker='d',linestyle="-")
    axes[1].set_title("Back Average by Frame"+stitle)
    axes[1].set_xlabel("Frame #")
    axes[1].set_ylabel("Back")
    
    plt.savefig(prefix+"_byframe2.png")
    if(inter == 1):
        fig.show()

    #third set - scale and rotation
    fig,axes=plt.subplots(1,2)
    fig.set_figheight(4)
    fig.set_figwidth(10)
 
    axes[0].plot(frames,sxAll,marker='d',linestyle="-")
    axes[0].plot(frames,syAll,marker='d',linestyle="-")
    axes[0].set_title("Scale Average by Frame"+stitle)
    axes[0].set_xlabel("Frame #")
    axes[0].set_ylabel("Scale")
    
    axes[1].plot(frames,rotAll,marker='d',linestyle="-")
    axes[1].set_title("Rotation Average by Frame"+stitle)
    axes[1].set_xlabel("Frame #")
    axes[1].set_ylabel("Rotation (radians)")

    plt.savefig(prefix+"_byframe3.png")
    if(inter == 1):
         fig.show()

    #fourth set - nper frame
         
def plotImageStats(image,prefix,inter,stitle=""):

    """

    plot histogram of an image
    
    """

    
    back = sigmaclip(image, sigma=2, iters=2)
    backImage=back.mean()
    rmsImage=back.std()

    logbins = np.geomspace(image.min(), image.max(), 50)
    
    fig,ax = plt.subplots()
    ax.hist(image.flatten(),bins=logbins,histtype="step")
    plt.title("Histogram of Region of Interest"+stitle)
    plt.xlabel("Flux Value")
    plt.ylabel("N")
    plt.yscale("log")
    plt.savefig(prefix+"_stats.png")

    if(inter == 1):
        fig.show()

    return backImage,rmsImage

def plotValHist(val,plotrange,titl,prefix,suffix,valtype,inter,nbins,stitle=""):

    #routine to plot a histogram of a variable

    fig,ax = plt.subplots()

    binsize=(plotrange[1]-plotrange[0])/nbins
    bins=np.arange(plotrange[0],plotrange[1]+binsize,binsize)
    ax.hist(val,bins=bins)
    plt.title(titl+stitle)
    plt.xlabel(valtype)
    plt.ylabel("N")
    plt.savefig(prefix+"_"+suffix+".png")

    if(inter == 1):
        fig.show()


def makeReg(x,y,outfile):

    """

    Dump a series of points to ds9 region file

    """

    ff=open(outfile,"w")
    for i in range(len(x)):
        print("circle point ",x[i],y[i],file=ff)
    ff.close()

def movieCutout(files,xc,yc,sz,mmin,mmax,prefix):

    """
    
    generate a series of cutouts around an xy point

    """
    
    fig,ax=plt.subplots()

    for i in range(len(files)):
        image=getImage(files[i])
        
        ax.imshow(image[xc-sz:xc+sz,yc-sz:yc+sz],vmin=mmin,vmax=mmax)
        plt.savefig(prefix+str(i).zfill(2)+".png")


    
def quiverPlot(x,y,dx,dy):

    """

    plot a distortion map in quiver and colour format

    Input:

      x,y - positions
      dx,dy - difference from expected
     
    returns: plot to file nad screen
    

    """

    
    fig,ax=plt.subplots(1,2)
    fig.set_figheight(4)
    fig.set_figwidth(10)

    ax[0].quiver(x,y,dx,dy)
    ax[0].set_xlabel("X (pixels)")
    ax[0].set_ylabel("Y (pixels)")

    ##map shows the total deviation
    
    dist=np.sqrt((dx)**2+(dy)**2)
    sc=ax[1].scatter(x,y,c=dist)
    ax[1].set_xlabel("X (pixels)")
    fig.colorbar(sc,ax=ax[1])

    fig.suptitle("Distortion Map")
    
    fig.show()
     
def quiverPlot1(x,y,dx,dy,frameID1,maxVal):

    """

    plot a distortion map in quiver and colour format

    Input:

      x,y - positions
      dx,dy - difference from expected
     
    returns: plot to file nad screen
    

    """

    
    fig,ax=plt.subplots(1,2)
    fig.set_figheight(4)
    fig.set_figwidth(10)

    ax[0].quiver(x,y,dx,dy)
    ax[0].set_xlabel("X (pixels)")
    ax[0].set_ylabel("Y (pixels)")

    ##map shows the total deviation
    
    dist=np.sqrt((dx)**2+(dy)**2)
    sc=ax[1].scatter(x,y,c=dist,vmin=dist.min(),vmax=maxVal)
    ax[1].set_xlabel("X (pixels)")
    fig.colorbar(sc,ax=ax[1])

    fig.suptitle("Distortion Map ("+str(int(frameID1))+")")
    
    
    fig.show()
    plt.savefig("quiver_"+str(int(frameID1))+".png")
    
def diagPlot(image,mCentroids,dx,dy):

    """

    plot the image wtih the distortion overplotted in images. 

    This needs a (working) functino to get the fibreID by clicking on it a point 

    Input: 
       iimage - image
       mCentroids - array with matched centroids
       dx,dy - distortion

    """

    
    
    fig,ax=plt.subplots()
    x=mCentroids[:,1]
    y=mCentroids[:,2]
    ax.imshow(image,origin='lower')
    ax.set_xlim([x.min(),x.max()])
    ax.set_ylim([y.min(),y.max()])
    ax.quiver(x,y,dx,dy,color='white',headlength=0, headaxislength=0)

    #this shoudl read the position, but doesn't work yet.
    def onM(event):
        xpos=event.xdata
        ypos=event.ydata
        print(xpos,ypos)
        fig.canvas.draw()
    
    fig.canvas.mpl_connect('button_press_event', onM)
    plt.show()

    
def checkThreshold(image,xrange,yrange):

    """

    quick polot to show the image and overplot the region for the threshold calculation

    """
    
    fig,ax=plt.subplots()
    ax.imshow(image,origin="lower")
    ax.scatter([yrange[0],yrange[0],yrange[1],yrange[1]],[xrange[0],xrange[1],xrange[0],xrange[1]])
    fig.show()

def quiverChange(frameIDs,centroidFile):
    centroids=np.loadtxt(centroidFile)
    ind=np.where(centroids[:,0]==frameIDs[0])
    xfirst=centroids[ind,2].ravel()
    yfirst=centroids[ind,3].ravel()

    
    #match all the frames to the first frame
    tol=20
    xArray,yArray,fxArray,fyArray,backArray,peakArray,qualArray=vis.matchAllPoints(centroids,xfirst,yfirst,tol,frameIDs)
    
    #get transformations by frame
    xdAll,ydAll,sxAll,syAll,rotAll,fxFrameAv,fyFrameAv,peakFrameAv,transAll = vis.getTransByFrame(xArray,yArray,fxArray,fyArray,peakArray,xfirst,yfirst)
    #sAll=(sxAll+syAll)/2
    xAv,yAv,fxAv,fyAv,peakAv,backAv,rmsVal,nMatch,xArray1,yArray1,dd,rmsX,rmsY,xd,yd = vis.getRMSStats(xArray,yArray,fxArray,fyArray,peakArray,backArray,xdAll,ydAll,sxAll,syAll,rotAll,xfirst,yfirst)

    nFrames=len(frameIDs)

    plt.ioff()
    xd=[]
    yd=[]
    for i in range(nFrames-1):
        xd.append(xArray[:,i]-xArray[:,i+1])
        yd.append(yArray[:,i]-yArray[:,i+1])

    xd=np.array(xd)
    yd=np.array(yd)


    for i in range(nFrames-1):
    #for i in range(1):
        fig,ax=plt.subplots(1,2)
        fig.set_figheight(4)
        fig.set_figwidth(10)

        dd=np.sqrt(xd[i]*xd[i]+yd[i]*yd[i])
        ind=np.where(dd < 3)

        if(i==0):
            
            Q=ax[0].quiver(xArray[ind,0],yArray[ind,0],xd[i,ind],yd[i,ind],scale=None,scale_units=None)
            ax[0].set_xlabel("X (pixels)")
            ax[0].set_ylabel("Y (pixels)")


            ax[0].quiverkey(Q,0.9, 0.98, 2,"2 pixels")

            dist=np.sqrt((xd[i,ind])**2+(yd[i,ind])**2)
            cmin=dist.min()
            cmax=dist.max()
            sc=ax[1].scatter(xArray[ind,0],yArray[ind,0],c=dist,vmin=cmin,vmax=cmax)
            ax[1].set_xlabel("X (pixels)")
            fig.colorbar(sc,ax=ax[1])

            
        else:
            ax[0].quiver(xArray[ind,0],yArray[ind,0],xd[i,ind],yd[i,ind],scale=Q.scale,scale_units=None)
            ax[0].set_xlabel("X (pixels)")
            ax[0].set_ylabel("Y (pixels)")
            ax[0].quiverkey(Q,0.9, 0.98, 2,"2 pixels")

            dist=np.sqrt((xd[i,ind])**2+(yd[i,ind])**2)
            sc=ax[1].scatter(xArray1[ind,0],yArray[ind,0],c=dist,vmin=cmin,vmax=cmax)
            ax[1].set_xlabel("X (pixels)")
            fig.colorbar(sc,ax=ax[1])

        plt.savefig("quiv1_"+str(int(frameIDs[0]))+"_"+str(int(i)).zfill(2)+".png")

    ##############################################



def quiverChange1(frameIDs,centroidFile):

    centroids=np.loadtxt(centroidFile)
    ind=np.where(centroids[:,0]==frameIDs[0])
    xfirst=centroids[ind,2].ravel()
    yfirst=centroids[ind,3].ravel()

    
    #match all the frames to the first frame
    tol=20
    xArray,yArray,fxArray,fyArray,backArray,peakArray,qualArray=vis.matchAllPoints(centroids,xfirst,yfirst,tol,frameIDs)
    
    #get transformations by frame
    xdAll,ydAll,sxAll,syAll,rotAll,fxFrameAv,fyFrameAv,peakFrameAv,transAll = vis.getTransByFrame(xArray,yArray,fxArray,fyArray,peakArray,xfirst,yfirst)
    #sAll=(sxAll+syAll)/2
    xAv,yAv,fxAv,fyAv,peakAv,backAv,rmsVal,nMatch,xArray1,yArray1,dd,rmsX,rmsY,xd,yd = vis.getRMSStats(xArray,yArray,fxArray,fyArray,peakArray,backArray,xdAll,ydAll,sxAll,syAll,rotAll,xfirst,yfirst)

    nFrames=len(frameIDs)
     
    for i in range(nFrames-1):
    #for i in range(1):

        transform,xd,yd,sx,sy,rotation=vis.getTransform(xArray[:,i],yArray[:,i],xArray[:,i+1],yArray[:,i+1],1)
        xArray1, yArray1 = vis.transformPointsNew(xArray[:,i],yArray[:,i],xd,yd,rotation,sx,sy,)
        
        xd=xArray1-xArray[:,i+1]
        yd=yArray1-yArray[:,i+1]
        dd=np.sqrt(xd*xd+yd*yd)
        ind=np.where(dd < 3)

        fig,ax=plt.subplots(1,2)
        fig.set_figheight(4)
        fig.set_figwidth(10)


        if(i==0):
            
            #Q=ax[0].quiver(xArray[ind,0],yArray[ind,0],xd[ind],yd[ind],scale=None,scale_units=None)
            Q=ax[0].quiver(xArray[ind,0],yArray[ind,0],xd[ind],yd[ind],scale=20)
            ax[0].set_xlabel("X (pixels)")
            ax[0].set_ylabel("Y (pixels)")


            ax[0].quiverkey(Q,0.9, 0.98, 2,"2 pixels")

            dist=np.sqrt((xd[ind])**2+(yd[ind])**2)
            cmin=0
            cmax=0.3

            aa=xArray[ind,0].flatten()
            bb=yArray[ind,0].flatten()
            sc=ax[1].scatter(aa,bb,c=dist,vmin=cmin,vmax=cmax)
            ax[1].set_xlabel("X (pixels)")
            fig.colorbar(sc,ax=ax[1])

            
        else:
            ax[0].quiver(xArray[ind,0],yArray[ind,0],xd[ind],yd[ind],scale=20)
            #Q=ax[0].quiver(xArray[ind,0],yArray[ind,0],xd[ind],yd[ind],scale=None,scale_units=None)
            ax[0].set_xlabel("X (pixels)")
            ax[0].set_ylabel("Y (pixels)")
            #ax[0].quiverkey(Q,0.9, 0.98, 2,"2 pixels")

            dist=np.sqrt((xd[ind])**2+(yd[ind])**2)
            aa=xArray[ind,0].flatten()
            bb=yArray[ind,0].flatten()
            sc=ax[1].scatter(aa,bb,c=dist,vmin=cmin,vmax=cmax)
            ax[1].set_xlabel("X (pixels)")
            fig.colorbar(sc,ax=ax[1])

        plt.savefig("quiv2_"+str(int(frameIDs[0]))+"_"+str(int(i)).zfill(2)+".png")
 
