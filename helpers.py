#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May 18 12:07:06 2017

@author: max
"""
import numpy as np
from math import sqrt,acos,degrees,sin,cos,atan2,pi
from shapely.geometry import LineString,Point

# PLOTTING HELPERS
def plot_line(ax, ob):
    x, y = ob.xy
    ax.plot(x, y, linewidth=3, solid_capstyle='round', zorder=2)

def plot_lines(ax,line1,line2,line3,line4):
    plot_line(ax, line1)
    plot_line(ax, line2)
    plot_line(ax, line3)
    plot_line(ax, line4)
        
def plot_mesh(ax,mesh):
    for row in mesh:
        xgrid = ([row[1],row[3],row[5],row[7],row[1]])
        ygrid = ([row[2],row[4],row[6],row[8],row[2]])  
        ax.plot(xgrid,ygrid,'black')
        
def plot_verts(ax,x,y,rmeshnum,pmeshnum):
    for j in range(0,pmeshnum):
        for i in range(0,rmeshnum):
            ax.plot(x[i,j],y[i,j],"o")
            
def plot_xpt(ax,xpt):
    ax.plot(xpt[0],xpt[1],"o")
    
#GEOMETRY HELPERS
def getangle3ptsdeg(p1,p2,p3):
    a = sqrt((p1[0]-p2[0])**2+(p1[1]-p2[1])**2)
    b = sqrt((p2[0]-p3[0])**2+(p2[1]-p3[1])**2)
    c = sqrt((p1[0]-p3[0])**2+(p1[1]-p3[1])**2) 
    #print ("a,b,c = ",a,b,c)
    #print ("line2 = ",(c**2 - a**2 - b**2)/(-2*a*b))
    theta = degrees(acos((c**2 - a**2 - b**2)/(-2*a*b))) #returns degree in radians
    return theta

def makeline(pt,r,theta):
    #print (pt,r,theta)
    pt2 = [r*cos(theta) + pt[0],r*sin(theta) + pt[1]]
    line = LineString([(pt[0],pt[1]), (pt2[0],pt2[1])])
    return pt2,line

def getangle(p1,p2):
    if isinstance(p1,Point) and isinstance(p2,Point):
        p1 = [p1.coords.xy[0][0],p1.coords.xy[1][0]]
        p2 = [p2.coords.xy[0][0],p2.coords.xy[1][0]]
    p1 = np.asarray(p1)
    p1 = np.reshape(p1,(-1,2))
    p2 = np.asarray(p2)
    p2 = np.reshape(p2,(-1,2))
    theta = np.arctan2(p1[:,1]-p2[:,1],p1[:,0]-p2[:,0])
    theta_mod = np.where(theta<0,theta+pi,theta) #makes it so the angle is always measured counterclockwise from the horizontal
    return theta
    
def getpoint(p1,angle,distance,shapely):
    if isinstance(p1,Point):
        p1 = [p1.coords.xy[0][0],p1.coords.xy[1][0]]
    p1 = np.asarray(p1)
    angle = np.asarray(angle)
    distance = np.asarray(distance)

    p1 = np.reshape(p1,(-1,2))

    px = p1[:,0]+distance*np.cos(angle)
    py = p1[:,1]+distance*np.sin(angle)    
    if shapely=='s':
        point = Point(px,py)
    else:
        point = np.column_stack((px,py))
    return point
    
def getdist(p1,p2):
    x = p1[0]-p2[0]
    y = p1[1]-p2[1]
    dist = sqrt(x**2 + y**2)
    return dist
    
def getpc(pa,pb,lb,theta_a):
    xa,ya,xb,yb = pa[0],pa[1],pb[0],pb[1]
    la = sqrt((xa-xb)**2+(ya-yb)**2)
    xc = (lb/la)*((xa - xb)*cos(theta_a) - (ya - yb)*sin(theta_a)) + xb
    yc = (lb/la)*((ya - yb)*cos(theta_a) + (xa - xb)*sin(theta_a)) + yb
    pc = [xc,yc]
    return pc
    
def getangle3pts(p1,p2,p3):
    a = sqrt((p1[0]-p2[0])**2+(p1[1]-p2[1])**2)
    b = sqrt((p2[0]-p3[0])**2+(p2[1]-p3[1])**2)
    c = sqrt((p1[0]-p3[0])**2+(p1[1]-p3[1])**2) 
    theta = acos(round((c**2 - a**2 - b**2)/(-2*a*b),5)) #returns degree in radians
    return theta
    
def cut(line, distance):
    # Cuts a line in two at a distance from its starting point
    if distance <= 0.0 or distance >= 1.0:
        return [LineString(line)]
    coords = list(line.coords)
    for i, p in enumerate(coords):
        pd = line.project(Point(p),normalized=True)
        if pd == distance:
            return [
                LineString(coords[:i+1]),
                LineString(coords[i:])]
        if pd > distance:
            cp = line.interpolate(distance,normalized=True)
            return [
                LineString(coords[:i] + [(cp.x, cp.y)]),
                LineString([(cp.x, cp.y)] + coords[i:])]
            
def PolyArea(x,y):
    return 0.5*np.abs(np.dot(x,np.roll(y,1))-np.dot(y,np.roll(x,1)))

#def FindAngle(x1,y1,x2,y2,x3,y3):
#    adotb = (lim_exp[:,0]-np.roll(lim_exp[:,0],1))*(lim_exp[:,0]-np.roll(lim_exp[:,0],-1)) + (lim_exp[:,1]-np.roll(lim_exp[:,1],1))*(lim_exp[:,1]-np.roll(lim_exp[:,1],-1))
#    mag_a = np.sqrt((lim_exp[:,0]-np.roll(lim_exp[:,0],1))**2+(lim_exp[:,1]-np.roll(lim_exp[:,1],1))**2)
#    mag_b = np.sqrt((lim_exp[:,0]-np.roll(lim_exp[:,0],-1))**2+(lim_exp[:,1]-np.roll(lim_exp[:,1],-1))**2)
#    return 0.5*np.abs(np.dot(x,np.roll(y,1))-np.dot(y,np.roll(x,1)))
            
#Physics helpers
def SOL_width(ID1_width,OD4_width,soltheta,ib_mp_width,ob_mp_width):
    index2 = np.where(soltheta==0)[0][0] + 1
    index3 = np.where(soltheta==pi)[0][0] + 1
    
    region1 = np.delete(np.linspace(OD4_width,ob_mp_width,index2),1)
    region23 = np.delete(np.linspace(ob_mp_width,ib_mp_width,index3-index2+1),1)
    region4 = np.linspace(ib_mp_width,ID1_width,len(soltheta)-index3+1)
    
    widths = np.concatenate((region1, region23, region4), axis=0)
    return widths

def isinline(pt,line):
    pt_s = Point(pt[0],pt[1])
    dist = line.distance(pt_s)
    if dist < 1E-6:
        return True
    else:
        return False