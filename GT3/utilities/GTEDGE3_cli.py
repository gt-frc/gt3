#!/usr/bin/python2.7

#########################################
#
#   GTEDGE V3
#  
#  Main CLI script
#
#########################################


from GT3.gt3 import gt3
from GT3.GT3Prep import gt3Prep
import sys, argparse, warnings, datetime, os
import deprecation



def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

@deprecation.deprecated(deprecated_in="0.0.3", details="GT3 no longer utilizes this CLI interface")

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

@deprecation.deprecated(deprecated_in="0.0.3", details="GT3 no longer utilizes this CLI interface")
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


    print("shotid=%s   runid=%s    timeid=%s   IOL Correction=%s" % (str(shotid),str(runid),str(timeid),str(IOL)))
    maxfile='togt3_d3d_'+str(shotid)+'_'+str(timeid)
    gt3Prep(shotid, timeid, runid, maxfile, quiet = quiet, genFiles=False)


    if reNeut==True:
        try:
            os.remove(os.path.join('inputs','%s_%s','gt3_%s_%s_neut.dat') % (str(shotid),(timeid), str(shotid),(timeid)))
        except:
            pass

    myPlasma = gt3(shotlabel=maxfile, mode=gt3Method, iolFlag = IOL, neutFlag = neutrals)


    return myPlasma