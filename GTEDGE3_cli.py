#!/usr/bin/python2.7

#########################################
#
#   GTEDGE V3
#  
#  Main CLI script
#
#########################################

import matplotlib.pyplot as plt
#import lib.beams as beams
#import lib.funcs as funcs
import graphs as graphs
from gt3 import gt3, gt3Prep
#from lib.funcs.dataGen import newScatPlot,smooth,multiScatPlot
#from scipy.interpolate import UnivariateSpline
from numpy import interp, pi
from math import sqrt
import inspect, sys, argparse,warnings,datetime,os
import numpy as np
from scipy import interpolate

##########################################################
#
#    Module to catologue file names, graphs, etc.
#
#   When a new shot is added, include its details here
#   where necessary (e.g., vpolID for bad pertrubation data
#
##########################################################


def CatalogueCall(shotid,timeid,runid):
    fileCat={}
    fileCat["GTEDGEsupp"]="Inputs/GTEDGEsupp_"+str(shotid)+"_"+str(timeid)+".csv"
    fileCat["GT3Consts"]="Inputs/GT3Consts_"+str(shotid)+"_"+str(timeid)+".csv"
    fileCat["GT3Profs"]="Inputs/GT3Profs_"+str(shotid)+"_"+str(timeid)+".csv"
    fileCat["GT3NBI"]="Inputs/GT3NBI_"+str(shotid)+"_"+str(timeid)+".csv"

    fileCat["Erspl"]="p"+str(shotid)+"_"+str(timeid)+"_"+str(runid)+"_ersplrhob_sc.txt"
    fileCat["toplasma"]=str(shotid)+"_"+str(timeid)+"_toplasma"
    fileCat["bpol"]=str(shotid)+"_"+str(timeid)+"_"+str(runid)+"_bpol.txt"
    fileCat["btor"]=str(shotid)+"_"+str(timeid)+"_"+str(runid)+"_btor.txt"
    fileCat["R"]=str(shotid)+"_"+str(timeid)+"_"+str(runid)+"_Z.txt"
    fileCat["BDRY"]=str(shotid)+"_"+str(timeid)+"_"+str(runid)+"_BDRY.txt"
    fileCat["EPOTEN"]=str(shotid)+"_"+str(timeid)+"_"+str(runid)+"_EPOTEN.txt"
    fileCat["lim"]=str(shotid)+"_"+str(timeid)+"_"+str(runid)+"_lim.txt"
    fileCat["psirz"]=str(shotid)+"_"+str(timeid)+"_"+str(runid)+"_psirz.txt"
    fileCat["Z"]=str(shotid)+"_"+str(timeid)+"_"+str(runid)+"_Z.txt"

    vTorPoints={(166606,1950,"j8099"):[198],
            (144977,3000,"j3000"):[193],
            (118890,1515,"r90"):[181],
            (118890,1560,"r90"):[191]}

    try:
        fileCat["vtorID"]=vTorPoints[(shotid,timeid,runid)]
    except:
        pass
    return fileCat

def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

def calc_fsa(x, R, Z):
    R1 = R[:, :-1]
    R2 = np.roll(R[:, :-1], -1, axis=1)
    Z1 = Z[:, :-1]
    Z2 = np.roll(Z[:, :-1], -1, axis=1)
    x1 = x[:, :-1]
    x2 = np.roll(x[:, :-1], -1, axis=1)

    dl = np.sqrt((R2 - R1)**2 + (Z2 - Z1)**2)

    R_av = (R1 + R2)/2

    dA = dl * (2 * pi * R_av)

    x_av = (x1 + x2)/2

    fsa = np.sum(x_av * dA, axis=1) / np.sum(dA, axis=1)
    fsa[0] = x[0,0]
    return fsa

def argFunc(a):
    argDic={}
    try:
        argDic['shotid']=int(a.shotid)
    except:
        raise NameError("shotid is not formatted correctly")
    try:
        argDic['timeid']=int(a.timeid)
    except:
        raise NameError("timeid is not formatted correctly")
    try:
        argDic['runid']=str(a.runid)
    except:
        raise NameError("runid is not formatted corectly")
    argDic['nbRun']=a.nbRun
    argDic['IOL']=a.IOL
    argDic['quiet']=a.q
    argDic['gt3Method']=a.gt3Method
    argDic['reNeut']=a.reNeut
    
    return argDic

def customwarn(message, category, filename, lineno, file=None, line=None):
    sys.stderr.write(warnings.formatwarning(message, category, filename, lineno))
    
def maxDebug(coreData):
    
    maxiandic={}
    maxiandic['rhor']=coreData.rhor
    maxiandic['fiol']=coreData.fiol
    maxiandic['eiol']=coreData.eiol
    maxiandic['miol']=coreData.miol
    maxiandic['xti']=coreData.xti
    maxiandic['xni']=coreData.xni
    maxiandic['xte']=coreData.xne
    maxiandic['fiolC']=coreData.fiolC
    maxiandic['eiolC']=coreData.eiolC
    maxiandic['miolC']=coreData.miolC
    maxiandic['xnC']=coreData.xnC
    maxiandic['xre']=coreData.xer
    maxList=maxiandic.keys()
    
    funcs.dataGen.csvDump("fioldebug.csv",maxiandic)
    

#def run(shotid,timeid,runid,nbRun=True,IOL=True):
def run(shotargs):
    ###########################################################################
    #
    #    New GTEDGE 3 with full core data and interpretation
    #
    #   Data provided as 201-point CSV files with the following naming conv.:
    #
    #   Input/GT3Profs_shotid_timeid.csv                Profiles w/ Core
    #   Input/GT3Consts_shotid_timeid.csv               Constants
    #   Input/nudraginputs_shotid_timeid_2v2.csv        GTEDGE profiles
    #   Input/GT3NBIConsts_144977_3000                  NBI-Relevent Constants
    #
    ############################################################################
    shotid=shotargs['shotid']
    runid=shotargs['runid']
    timeid=shotargs['timeid']
    IOL=shotargs['IOL']
    nbRun=shotargs['nbRun']
    quiet=shotargs['quiet']
    gt3Method=shotargs['gt3Method']
    reNeut=shotargs['reNeut']
    
    errLog=open("GT3err.%s.%s.log" % (shotid,timeid),"w+")
    errLog.write("\n")
    errLog.write("Time: %s \n" % str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
    errLog.write("%s.%s.%s  nbRun=%s IOL=%s \n" % (str(shotid),str(timeid),str(runid),str(nbRun),str(IOL)))
    
    sys.stderr=errLog
    warnings.simplefilter('always', UserWarning)
    warnings.simplefilter('always', RuntimeWarning)
    warnings.showwarning=customwarn
    


#    fileNames=funcs.dataCat.CatalogueCall(shotid,timeid,runid)

    
    GTedgeProfsFile=fileNames['GTEDGEsupp']
    constsFile=fileNames['GT3Consts']
    coreProfsFile=fileNames['GT3Profs']
    nbiFile=fileNames['GT3NBI']
    
#   GTEDGE Data
    
    GTEDGEdata=funcs.dataGen.dataGrab(profs=GTedgeProfsFile,const=constsFile)
    coreDataDic=funcs.dataGen.dataGrab(profs=coreProfsFile,const=constsFile)
    beamDic=funcs.dataGen.dataGrab(const=nbiFile)
    beamParams=funcs.dataGen.dataAug(beamDic)

    ###########################################################################
    #
    #   Augment coreData to allow things such as
    #   coreData.xni = []
    #   coreData.xti = []
    #   coreData.aminor = #.#
    #
    #   coreData will include 1) the coreDataDic dictionary, which includes
    #   all the mathematical functions, and 2) the above augmenting data
    #
    #   This removes the need to do things like coreData['key'].values()
    #
    ###########################################################################
    
    coreData=funcs.dataGen.dataAug(coreDataDic)
    coreData.rhor=[0.+x*(1./(len(coreData.xte)-1)) for x in range(len(coreData.xte))]    
    #   Calculate carbon density from experimental data
    #   Set zbar2 = 6. if data not provided
    try:
        coreData.zbar2
    except:
        coreData.zbar2=[6.]*201
    coreData.xnC=funcs.physFuncs.xnCcalc(coreData.fz1,coreData.zbar2)                    # Possibly not legit    
    #      Set q=4. if data not provided
    #   Note: This is "qedge" legacy from GTEDGE and should just be q
    #   in reality
    try:
        coreData.q
    except:
        coreData.q=[4.]*201
    
    coreData.xnuc=funcs.physFuncs.xnucCalc(coreData)

    ###################################################################
    #	NeutPy
    #
    #   Runs at low resolution for the plasma to get neutrals before
    #   running at full resolution
    #
    ###################################################################
    
    print "shotid=%s   runid=%s    timeid=%s   IOL Correction=%s" % (str(shotid),str(runid),str(timeid),str(IOL))
    maxfile='togt3_d3d_'+str(shotid)+'_'+str(timeid)
    funcs.maxPrep.maxPrep(shotid,timeid,runid,maxfile,coreData,beamParams,reNeut)

    #beams.run(nbiFile,coreData.xte,coreData.xti,coreData.xne,coreData.xni,coreData.xnC,coreData.zbar2,coreData.rmajor,coreData.elong,fastiolzip,run=False)
    os.chdir("MaxPlasma") 

    #   Delete neutrals data if rerun flag (reNeut) is toggled True    
    if reNeut==True:
        try:
            os.remove('gt3_%s_%s_neut.dat' % (str(shotid),(timeid)))
        except:
            pass



    ###########################################################
    #
    #    Debug
    #
    ###########################################################
#    funcs.dataGen.newScatPlot(coreData.rhor,coreData.xte,ylabel="xte")

    myPlasma=gt3(maxfile) 
    maxMethodDic={
#        'coreonly':myPlasma.coreonly,
        'coreandiol':myPlasma.coreandiol,
        'coreiolnbi':myPlasma.coreiolnbi,
        'coreiolnbineuts':myPlasma.coreiolnbineuts,
        'corenbi':myPlasma.coreandnbi,
        'corenbineuts':myPlasma.corenbineuts
        }
    if "neuts"in gt3Method:
        if reNeut==True:
            maxMethodDic[gt3Method](ntrl_switch=2)
        else:
            maxMethodDic[gt3Method](ntrl_switch=1)
    else:
        maxMethodDic[gt3Method]()
    

    os.chdir("..")

    # [:1] is simply taking the theta = 1 distribution 

    #Are neutrals available?
     
    try:
        
#        coreData.xnuati=map(lambda x: x*2.,myPlasma.ntrl.nn_xnuioni_1D[:,1])
#        coreData.Lz = myPlasma.imp.Lz
#        coreData.gt3=myPlasma
        ionRateInp = interpolate.interp1d(myPlasma.core.rho[:, 0], myPlasma.core.izn_rate.tot[:, 0])
        coreData.ionRate = np.array(map(lambda x: ionRateInp(x), coreData.rhor))

#        xnuioniInp = interpolate.interp1d(myPlasma.core.rho[:, 0],calc_fsa(myPlasma.core.sv.ion.st * myPlasma.core.n.n.t, myPlasma.core.R, myPlasma.core.Z))

        # Try to use theta = 8 to minimize neutrals density
        xnuioniInp = interpolate.interp1d(myPlasma.core.rho[:, 0], myPlasma.core.sv.ion.st[:, 8]*myPlasma.core.n.n.t[:, 8]
                                          )
        sigvionInp = interpolate.interp1d(myPlasma.core.rho[:, 0], myPlasma.core.sv.ion.st[:, 0])

        coreData.xnuioni=np.array(map(lambda x: xnuioniInp(x),coreData.rhor))
        coreData.sigvion = np.array(map(lambda x: sigvionInp(x), coreData.rhor))
        coreData.sigvata = myPlasma.core.sv.el.st + myPlasma.core.sv.cx.st
#        xnuatiInp = interpolate.interp1d(myPlasma.core.rho[:, 0], calc_fsa(coreData.sigvata * myPlasma.core.n.n.s,myPlasma.core.R ,myPlasma.core.Z))

        # Try to use theta = 8 to minimize neutrals density
        xnuatiInp = interpolate.interp1d(myPlasma.core.rho[:, 0], (myPlasma.core.sv.cx.st + myPlasma.core.sv.el.st)[:,8]*myPlasma.core.n.n.s[:,8])
        coreData.xnuati = np.array(map(lambda x: xnuatiInp(x), coreData.rhor))

#        coreData.xnuioni=myPlasma.ntrl.nn_xnuioni[:,1] # IN RZ

        radcoolInp = interpolate.interp1d(myPlasma.core.rho[:, 0],calc_fsa(myPlasma.core.n.e * myPlasma.core.n.C * np.nan_to_num(myPlasma.core.Lz_thermal[0]) +
                                                                           myPlasma.core.n.e * myPlasma.core.n.C * np.nan_to_num(myPlasma.core.Lz_slow[0]),
                                                                           myPlasma.core.R, myPlasma.core.Z))
        coreData.radcool=np.array(map(lambda x: radcoolInp(x), coreData.rhor))
        # Kill off core neutrals because BS

        coreData.xnuati[:180]=np.array([0.]*180)
        coreData.xnuioni[:180]=np.array([0.]*180)


        if "neuts" not in gt3Method:
            coreData.xnuati,coreData.xnuioni=np.zeros(len(coreData.xne)),np.zeros(len(coreData.xne))
    except:
        print "No neutrals data or crashed"
        coreData.xnuioni,coreData.xnuati,coreData.radcool=np.array([0.]*len(coreData.xne)), np.array([0.]*len(coreData.xne)), np.array([0.]*len(coreData.xne))
    ###################################################################
    #
    #   Display summary of shot being run
    #   Will error if shot/run IDs do not match between Max's code
    #   and this code
    #   Also takes in flag to say if NBeams is running
    #
    ###################################################################
    funcs.dataGen.runSummary(funcs.dataGen.shotID(coreProfsFile,fileNames['Erspl']),quiet=quiet,IOLflag=IOL)
    
    ###################################################################
    #
    #   IOL Calculations
    #
    ###################################################################
#    maxrhor=[0.+x*(1./(myPlasma.inp.rhopts-1)) for x in range(myPlasma.inp.rhopts)]
    maxrhor=np.array(myPlasma.core.rho[:,0])
        
#    coreData.bthet=myPlasma
#    coreData.bphi=myPlasma
    
    ep = np.array([a*coreData.aminor/coreData.rmajor for a in coreData.rhor])
    coreData.bthet=np.array([a * abs(coreData.bphi) / coreData.q95 for a in ep])
    
    if IOL==False:
        coreData.fiol,coreData.miol,coreData.eiol=np.zeros((3,myPlasma.inp.rhopts))
        coreData.fiolT,coreData.miolT,coreData.eiolT=np.zeros((3,myPlasma.inp.rhopts))
        coreData.fiolC,coreData.miolC,coreData.eiolC=np.zeros((3,myPlasma.inp.rhopts))

    else:
        coreData.fiol=interp(coreData.rhor,maxrhor,myPlasma.iol.forb_d_therm[:,0])
#        coreData.fiol=coreData.fiol.tolist()
        coreData.miol=interp(coreData.rhor,maxrhor,myPlasma.iol.morb_d_therm[:,0])
#        coreData.miol=coreData.miol.tolist()
        coreData.eiol=interp(coreData.rhor,maxrhor,myPlasma.iol.eorb_d_therm[:,0])

        coreData.fiolT=interp(coreData.rhor,maxrhor,myPlasma.iol.forb_t_therm[:,0])
#        coreData.fiolT=coreData.fiol.tolist()
        coreData.miolT=interp(coreData.rhor,maxrhor,myPlasma.iol.morb_t_therm[:,0])
#        coreData.miolT=coreData.miol.tolist()
        coreData.eiolT=interp(coreData.rhor,maxrhor,myPlasma.iol.eorb_t_therm[:,0])
#        coreData.eiolT=coreData.eiol.tolist()

#       Hydrogen IOL if ever needed

        # coreData.fiolH = interp(coreData.rhor, maxrhor, myPlasma.iol.forb_h_therm[:, 0])
        # #        coreData.fiolT=coreData.fiol.tolist()
        # coreData.miolH = interp(coreData.rhor, maxrhor, myPlasma.iol.morb_h_therm[:, 0])
        # #        coreData.miolT=coreData.miol.tolist()
        # coreData.eiolH = interp(coreData.rhor, maxrhor, myPlasma.iol.eorb_h_therm[:, 0])
        # #        coreData.eiolT=coreData.eiol.tolist()

        coreData.fiolC=interp(coreData.rhor,maxrhor,myPlasma.iol.forb_c_therm[:,0])
#        coreData.fiolC=coreData.fiol.tolist()
        coreData.miolC=interp(coreData.rhor,maxrhor,myPlasma.iol.morb_c_therm[:,0])
#        coreData.miolC=coreData.miol.tolist()
        coreData.eiolC=interp(coreData.rhor,maxrhor,myPlasma.iol.eorb_c_therm[:,0])

        # Trim IOL

        coreData.fiol[:160]=[0.]*160
        coreData.miol[:160]=[0.]*160
        coreData.eiol[:160]=[0.]*160

        coreData.fiolT[:160]=[0.]*160
        coreData.miolT[:160]=[0.]*160
        coreData.eiolT[:160]=[0.]*160

        coreData.fiolC[:160]=[0.]*160
        coreData.miolC[:160]=[0.]*160
        coreData.eiolC[:160]=[0.]*160


    

    
    ###################################################################
    #
    #   any IOL nan converted to 0
    #    
    ####################################################################        
    
    coreData.fiol=np.nan_to_num(coreData.fiol)
    coreData.eiol=np.nan_to_num(coreData.eiol)
    coreData.miol=np.nan_to_num(coreData.miol)
    
																					
    ###################################################################
    #
    #   NBeams call, set data to nbiDep attribute in coreData
    #   Data obtained as NBIresults are as follows:
    #
    #   NBIresults[0]=nbiDep list (0-2)
    #   NBIresults[1]=snbi distribution
    #   NBIresults[2]=qnbi distribution
    #   NBIresults[3]=qnbe distribution 
    #   NBIresults[4]=beam input toroidal momentum deterium
    #   NBIresults[5]=beam input toroidal momentum carbon
    ###################################################################    
    
    #   Need FIOL data from Max's plasma
    
    #NBIresults=beams.run(nbiFile,coreData.xte,coreData.xti,coreData.xne,coreData.xni,coreData.xnC,coreData.zbar2,coreData.rmajor,coreData.elong,fastiolzip,run=nbRun)

    beamDic=funcs.dataGen.dataGrab(const=nbiFile)

    #NBIresults=beams.beamSources(coreData, beams.beamDataClass(beamDic), myPlasma.nbi.dep_prof1[:, 1],myPlasma.nbi.dep_prof2[:, 1],myPlasma.nbi.dep_prof3[:, 1], fastiolzip)

#    raw_input("test")
    NBIresults=beams.beamSourcesgt3(coreData,myPlasma)

    coreData.Snbi=np.array(NBIresults[1])
    coreData.qnbi=np.array(NBIresults[2])
    coreData.qnbe=np.array(NBIresults[3])
    coreData.xmomtor1=np.array(NBIresults[4])
    coreData.xmomtor2=np.array(NBIresults[5])
    
    ###################################################################
    #
    #   CXR data to be generated here. Currently set zbar2=6.0
    #   
    #
    ###################################################################
    

    coreData.gamC=np.zeros(len(coreData.rhor))

    funcs.coreflux.fluxCore(coreData,iolFlag=True)
#    gammahat,qhatHeati,qhatHeate=funcs.coreflux.fluxCore(coreData,iolFlag=True)
    

    
#    funcs.dataGen.dictAppend(coreData,coreData.rhor,(('gamma',gamma),('qHeati',qHeati),('qHeate',qHeate),('gamC',coreData.gamC)))

    coreData.y11=np.array(funcs.y11.y11(coreData))
    coreData.y22=np.array(funcs.y22.y22(coreData))
    
    ###################################################################
    #
    #    Intrinsic Rotation
    #
    ###################################################################
    
    rlossiol=0.5

    coreData.intrin=np.array([2./sqrt(pi)*morbl*sqrt(2*coreData.xk*Temp/coreData.xmas1) for morbl,Temp in zip(coreData.miol,coreData.xti)])
#    coreData.intrinC=[2./sqrt(pi)*morbl*sqrt(2*coreData.xk*Temp/coreData.xmas2) for morbl,Temp in zip(coreData.miolC,coreData.xti)]
    coreData.intrinC=np.zeros(len(coreData.xne))
    coreData.vtorC=np.array([a-b for a,b in zip(coreData.vtorC,coreData.intrinC)])

 
#    coreData.vtorChat=[a-b for a,b in zip(coreData.vtorC,coreData.yy2)]    
    
    ##################################################################
    #
    #   Calculate deuterium velocities with Carbon data or exp. data
    #
    ##################################################################
    
    nudragResults=funcs.nuDragMIOL.nuDragMIOL(6,coreData,y11=coreData.y11,y22=coreData.y22)
  
    if ((hasattr(coreData,'vtorD')==False) and (hasattr(coreData,'vtord')==False)):
        coreData.vtorD=nudragResults[2]
        coreData.vtorD=[a-b for a,b in zip(coreData.vtorD,coreData.intrin)]
    else:
        coreData.vtorD=[a-b for a,b in zip(coreData.vtorD,coreData.intrin)]
 
    ##################################################################
    #
    #   Remove bad points from crazy perturbation theory calculation
    #
    ##################################################################    
    
    try:
        for a in fileNames['vtorID']:
            funcs.dataGen.pointRemove(coreData.vtorD,a) 
            funcs.dataGen.pointRemove(coreData.vtorDhat,a)
    except:
        print "Failed to remove points"
        pass   
    coreData.nudrag=nudragResults[0]


#   TODO: Assuming Vpolc = 0.4*V_pold until better answer figured out
    
    coreData.vpolD=[a/.4 for a in coreData.vpolC]
    
    
    
#    chi=funcs.chi.chiClass(coreData)
    chi=funcs.chi.chiorbClass(coreData)  

    coreData.diff=funcs.diff.diffCalcs(coreData)
###############################################################################
#
#    Chi inference from expermental data
#
#   Chicalc(data,conv15=False,conv25=False,visc=False,pressure=False)
#
#   Setting an argument to True subtracts away that argument from qcond
#
###############################################################################
    
    coreData.chi1=chi.chiiCalc(coreData) 
    coreData.chi2=chi.chiiCalc(coreData) 
    coreData.chi3=chi.chiiCalc(coreData,conv25=True) 
    coreData.chi4=chi.chiiCalc(coreData,conv25=True,pressure=True) 
    coreData.chi5=chi.chiiCalc(coreData,conv25=True,pressure=True,visc=True)
    coreData.chie=chi.chieCalc(coreData)
    coreData.gameltemp=chi.gameltemp

    coreData.q1=chi.qtotal 
    coreData.q2=chi.qtotal 
    coreData.q3=chi.conv25 
    coreData.q4=chi.heatin 
    coreData.q5=chi.heatvisc 
    
    chiList=[coreData.chi1,coreData.chi2,coreData.chi3,coreData.chi4,coreData.chi5]
    qList=[coreData.q1,coreData.q2,coreData.q3,coreData.q4,coreData.q5]#    
    
    coreData.chiGraphsDic={}
    coreData.diffGraphsDic={}
    coreData.nuGraphsDic={}
    
    for name,obj in inspect.getmembers(graphs.chiGraphs):
        if inspect.isclass(obj):
            coreData.chiGraphsDic[obj.shotid]=obj
            
    for name,obj in inspect.getmembers(graphs.diffGraphs):
        if inspect.isclass(obj):
            coreData.diffGraphsDic[obj.shotid]=obj
            
    for name,obj in inspect.getmembers(graphs.nuGraphs):
        if inspect.isclass(obj):
            coreData.nuGraphsDic[obj.shotid]=obj
            




    #   Dump FIOL debug data
#    maxDebug(coreData)
#
    funcs.dataGen.variableFlushtoCSV("log.%s.%s.%s.csv" % (shotid,timeid,gt3Method),coreData)
#    plt.show(block=True)
    errLog.close()
    return coreData

def runGT3(shotargs):
    ###########################################################################
    #
    #   GTEDGE 3 with full core data and interpretation
    #   Uses gt3 as backend for background plasma calculations
    #
    #   Data provided as 201-point CSV files with the following naming conv.:
    #
    #   Input/GT3Profs_shotid_timeid.csv                Profiles w/ Core
    #   Input/GT3Consts_shotid_timeid.csv               Constants
    #   Input/nudraginputs_shotid_timeid_2v2.csv        GTEDGE profiles
    #   Input/GT3NBIConsts_144977_3000                  NBI-Relevent Constants
    #
    ############################################################################
    shotid = shotargs['shotid']
    runid = shotargs['runid']
    timeid = shotargs['timeid']
    IOL = shotargs['IOL']
    nbRun = shotargs['nbRun']
    quiet = shotargs['quiet']
    gt3Method = shotargs['gt3Method']
    reNeut = shotargs['reNeut']
    neutrals = shotargs['neutrals']
    debug = shotargs['debug']

    errLog = open("GT3err.%s.%s.log" % (shotid, timeid), "w+")
    errLog.write("\n")
    errLog.write("Time: %s \n" % str(datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S")))
    errLog.write("%s.%s.%s  nbRun=%s IOL=%s \n" % (str(shotid), str(timeid), str(runid), str(nbRun), str(IOL)))

    sys.stderr = errLog
    warnings.simplefilter('always', UserWarning)
    warnings.simplefilter('always', RuntimeWarning)
    warnings.showwarning = customwarn

    #fileNames = funcs.dataCat.CatalogueCall(shotid, timeid, runid)

    print "shotid=%s   runid=%s    timeid=%s   IOL Correction=%s" % (str(shotid),str(runid),str(timeid),str(IOL))
    maxfile='togt3_d3d_'+str(shotid)+'_'+str(timeid)
    gt3Prep(shotid,timeid,runid,maxfile,quiet = quiet, genFiles=False)


    if reNeut==True:
        try:
            os.remove(os.path.join('inputs','%s_%s','gt3_%s_%s_neut.dat') % (str(shotid),(timeid), str(shotid),(timeid)))
        except:
            pass

    myPlasma = gt3(shotlabel=maxfile, mode=gt3Method, iolFlag = IOL, neutFlag = neutrals, debugRT = debug)


    return myPlasma




###############################################################################
#
#
#    MAIN PROGRAM as CLI Interface
#
#	usage: GTEDGE3_cli.py [-h] [-nbRun [NBRUN]] [-IOL [IOL]] shotid timeid runid
#
#	positional arguments:
#	shotid          DIII-D shot id
#	timeid          DIII-D time id
#	runid           DIII-D run id
#
#	optional arguments:
#	-h, --help      show this help message and exit
#	-nbRun [NBRUN]  Run NBeams (default: True)
#	-IOL [IOL]      Correct for IOL (default: True)
#
#
###############################################################################   
    
if __name__== "__main__":


    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("shotid",help="DIII-D shot id")
    parser.add_argument("timeid",help="DIII-D time id")
    parser.add_argument("runid",help="DIII-D run id")
    parser.add_argument("-nbRun",type=str2bool,nargs='?',#const=True,
                        default=True,help="Run NBeams")
    parser.add_argument("-IOL",type=str2bool,nargs='?',#const=True,
                        default=True,help="Correct for IOL")
    parser.add_argument("-reNeut",type=str2bool,nargs='?',#const=True,
                        default=True,help="Rerun the neutrals")
    praser.add_argument("-neutrals", type=str2bool, nargs='?',
                        default=True, help="Calculates neutrals")
    parser.add_argument("-gt3Method",type=str2bool,nargs='?',#const="brndiolneuts",
                        default="brndiolneuts",help="What GT3 method to run")
    parser.add_argument("-debug", type=str2bool, nargs='?',
                        default="False", help="Produce Radial Transport debug profiles")
    parser.add_argument("-q",action='store_true',help="Disable informational popus")
    

    try:
        args=parser.parse_args()
    except SystemExit as err:
        if err.code==2:
            parser.print_help()
        sys.exit(0)
    run(argFunc(args))

