#!/usr/bin/python

##########################################################
#
#    Graphical presentation of various analysis results
#    
# 
##########################################################

#from matplotlib import rc
import matplotlib.pyplot as plt
from math import ceil,log
import numpy as np
import pandas as pd
import matplotlib.patches as mpatches


#from matplotlib import rc
#rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
## for Palatino and other serif fonts use:
#rc('font',**{'family':'serif','serif':['Palatino']})
#rc('text', usetex=True)

    
def prettyCompare(xdata,ydata,yrange=[0.,1.],ylabel="No Label",datalabels=False,
                  title="No Title",caption=False,floatCaption=None,
                  logSpace=False,corePlot=False,plotType="line",
                  adjust=None,textSpace=0.,size=(12,9),xTickremove=None,
                  yTickremove=None, marginLeft=None,marginRight=None,
                  marginBottom=None,marginTop=None,capAdj=0.,
                  xLabelAdj=0.,yLabelAdj=0., xlabel=r"Normalized minor radius ($\rho$)",
                  toSN=False,legend=False):

    if legend:textSpace=0.
    prettyData=[xdata]
    for a in ydata:
        prettyData.append(a)
    pandasFrame=pd.DataFrame.from_items(prettyData)

#    These are the "Tableau 20" colors as RGB.    

#    tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14), (255, 187, 120),    
#                 (44, 160, 44), (152, 223, 138), (214, 39, 40), (255, 152, 150),    
#                 (148, 103, 189), (197, 176, 213), (140, 86, 75), (196, 156, 148),    
#                 (227, 119, 194), (247, 182, 210), (127, 127, 127), (199, 199, 199),    
#                 (188, 189, 34), (219, 219, 141), (23, 190, 207), (158, 218, 229)]  

    tableau20 = [(255,0,0), (0,0,255), (0,255,0), (255,0,255),    
                 (128,0,0), (128,0,128), (0,0,128), (0,128,128)]              
    markerList=["*","x","s","o"]
    for i in range(len(tableau20)):    
        r, g, b = tableau20[i]    
        tableau20[i] = (r / 255., g / 255., b / 255.)      

    fig=plt.figure(figsize=size)
    
#   Remove the plot frame lines. They are unnecessary chartjunk.    
    ax = plt.subplot(111)    
    ax.spines["top"].set_visible(False)    
    ax.spines["bottom"].set_visible(False)    
    ax.spines["right"].set_visible(False)    
    ax.spines["left"].set_visible(False)    

    plt.subplots_adjust(bottom=marginBottom,left=marginLeft,top=marginTop,right=marginRight)
#   Ensure that the axis ticks only show up on the bottom and left of the plot.    
#   Ticks on the right and top of the plot are generally unnecessary chartjunk.    
    ax.get_xaxis().tick_bottom()    
    ax.get_yaxis().tick_left()   
      
#   Limit the range of the plot to only where the data is.    
#   Avoid unnecessary whitespace.    
    plt.ylim(yrange[0],yrange[1])    
    if corePlot==True:
        plt.xlim(0,0.85+textSpace)     # core
        xrange=[0,.85]
    else:
        plt.xlim(.9,1+textSpace)    # edge
        xrange=[.9,1.]

#   Make sure your axis ticks are large enough to be easily read.    
#   You don't want your viewers squinting to read your plot. 
    if logSpace==False:
        if toSN!=False:
            if yTickremove!=None:
                plt.yticks(np.linspace(yrange[0],yrange[1],num=6)[:(-1*yTickremove)], [str(round(x,2)) for x in np.linspace(yrange[0]/(1.*10**toSN),yrange[1]/(1.*10**toSN),num=6)], fontsize=30)    
            else:
#                print np.linspace(yrange[0]/(10**toSN),yrange[1]/(1.*10**toSN),num=6)
                plt.yticks(np.linspace(yrange[0],yrange[1],num=6), [str(round(x,2)) for x in np.linspace(yrange[0]/(1.*10**toSN),yrange[1]/(1.*10**toSN),num=6)], fontsize=26)    
        elif yTickremove!=None:
            plt.yticks(np.linspace(yrange[0],yrange[1],num=6)[:(-1*yTickremove)], [str(round(x,2)) for x in np.linspace(yrange[0],yrange[1],num=6)], fontsize=30)
        else:
            plt.yticks(np.linspace(yrange[0],yrange[1],num=6), [str(round(x,2)) for x in np.linspace(yrange[0],yrange[1],num=6)], fontsize=20)
            
    else:
        if yTickremove!=None:
            plt.yticks([log(a) for a in np.logspace(yrange[0],yrange[1],num=6)][:(-1*yTickremove)],[str(round(x,2)) for x in [log(a) for a in np.logspace(yrange[0],yrange[1],num=6)]], fontsize=26)
        else:
            plt.yticks([log(a) for a in np.logspace(yrange[0],yrange[1],num=6)],[str(round(x,2)) for x in [log(a) for a in np.logspace(yrange[0],yrange[1],num=6)]], fontsize=26)
    
    plt.tick_params(axis='x',pad=18)



    if xTickremove!=None:
        plt.xticks(np.linspace(xrange[0],xrange[1],num=6)[:(-1*xTickremove)], [str(round(x,2)) for x in np.linspace(xrange[0],xrange[1],num=6)], fontsize=26)    
    else:
        plt.xticks(np.linspace(xrange[0],xrange[1],num=6), [str(round(x,2)) for x in np.linspace(xrange[0],xrange[1],num=6)], fontsize=26)    

        
#   Provide tick lines across the plot to help your viewers trace along    
#   the axis ticks. Make sure that the lines are light and small so they    
#   don't obscure the primary data lines. 
    if logSpace==False:        
        for y in np.linspace(yrange[0],yrange[1],num=6):    
            plt.plot(np.linspace(xrange[0],xrange[1],num=6), [y] * len(np.linspace(xrange[0],xrange[1],num=6)), "--", lw=0.5, color="black", alpha=0.3)     
    else:
        for y in [log(a) for a in np.logspace(yrange[0],yrange[1],num=6)]:    
            plt.plot([log(a) for a in np.logspace(yrange[0],yrange[1],num=6)], [y] * len([log(a) for a in np.logspace(yrange[0],yrange[1],num=6)]), "--", lw=0.5, color="black", alpha=0.3)     
        
#   Remove the tick marks; they are unnecessary with the tick lines we just plotted.    
    plt.tick_params(axis="both", which="both", bottom=False, top=False,
                    labelbottom=True, left=False, right=False, labelleft=True)
    ystep=(yrange[1]-yrange[0])/25.
    if legend:
        plotList=[]
        patchList=[]
        for rank,column in enumerate(datalabels):
    #       Build the plots for plt.legend()
            x,=plt.plot(pandasFrame.rhor.values,
                         pandasFrame[column.replace("\n", " ")].values,
                        lw=2, color=tableau20[rank],marker=markerList[rank],markersize=10)
            plotList.append(x)
        plt.legend(plotList, datalabels, loc=2, fontsize=20,framealpha=1.)
    else:
        for rank, column in enumerate(datalabels):

    #       Plot each line separately with its own color, using the Tableau 20
    #       color set in order.
            if plotType=="line":
                plt.plot(pandasFrame.rhor.values,
                         pandasFrame[column.replace("\n", " ")].values,
                        lw=2.5, color=tableau20[rank])
                y_pos=pandasFrame[column.replace("\n", " ")].values[-1]
            elif plotType=="scat":
                plt.scatter(pandasFrame.rhor.values,
                        pandasFrame[column.replace("\n", " ")].values,
                        color=tableau20[rank])
                y_pos=pandasFrame[column.replace("\n", " ")].values[-1]
                if legend:
                    ax.legend(column, datalabels)
            else:
                raise ("No line type requested")
            try:
                if rank in list(adjust.keys()):
                    y_pos += adjust[rank]*ystep
            except:
                pass

        # Again, make sure that all labels are large enough to be easily read
        # by the viewer.
            #plt.text(xrange[1]+0.0025, y_pos, column, fontsize=16, color=tableau20[rank]) # Default
            plt.text(xrange[1] + 0.0025, y_pos, column, fontsize=26, color=tableau20[rank])

    # matplotlib's title() call centers the title on the plot, but not the graph,    
    # so I used the text() call to customize where the title goes.    
      
    # Make the title big enough so it spans the entire plot, but don't make it    
    # so big that it requires two lines to show.    
      
    # Note that if the title is descriptive enough, it is unnecessary to include    
    # axis labels; they are self-evident, in this plot's case.

#    plt.text((xrange[1]-xrange[0])/2.+xrange[0]+0.015, (yrange[1]-yrange[0])/20.+yrange[1], title, fontsize=26, ha="center")      
    plt.text((xrange[1]-xrange[0])/2.+xrange[0]+.5*textSpace, (yrange[1]-yrange[0])/20.+yrange[1], title, fontsize=30, ha="center")

    if caption!=False:
        for num,cap in enumerate(caption):
            plt.text(xrange[0],yrange[0]-(yrange[1]-yrange[0])/5.+capAdj*(yrange[1]-yrange[0])/10.-(num)*(yrange[1]-yrange[0])/22.5, cap, fontsize=16,ha="left")
    if floatCaption!=None:
        plt.text(floatCaption[1][0],floatCaption[1][1],floatCaption[0],fontsize=22,ha="left",va="top")
        
    if not legend:
        for num,label in enumerate(ylabel):
            plt.text(xrange[0]-(xrange[1]-xrange[0])/10.+yLabelAdj*(xrange[1]-xrange[0])/10.,(yrange[1]-yrange[0])/2.+yrange[0]-num*(yrange[1]-yrange[0])/10.,label,fontsize=36,ha="center") # default


    #plt.text((xrange[1] - xrange[0]) / 2. + xrange[0],yrange[0] - (yrange[1] - yrange[0]) / 9. + xLabelAdj * (yrange[1] - yrange[0] / 10.), xlabel, fontsize=16, ha="center") #Default
    plt.text((xrange[1] - xrange[0]) / 2. + xrange[0],yrange[0] - (yrange[1] - yrange[0]) / 9. + xLabelAdj * (yrange[1] - yrange[0] / 10.), xlabel, fontsize=20, ha="center")