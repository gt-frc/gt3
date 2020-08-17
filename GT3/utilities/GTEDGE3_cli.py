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



def str2bool(v):
    if v.lower() in ('yes', 'true', 't', 'y', '1'):
        return True
    elif v.lower() in ('no', 'false', 'f', 'n', '0'):
        return False
    else:
        raise argparse.ArgumentTypeError('Boolean value expected.')

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


    print "shotid=%s   runid=%s    timeid=%s   IOL Correction=%s" % (str(shotid),str(runid),str(timeid),str(IOL))
    maxfile='togt3_d3d_'+str(shotid)+'_'+str(timeid)
    gt3Prep(shotid, timeid, runid, maxfile, quiet = quiet, genFiles=False)


    if reNeut==True:
        try:
            os.remove(os.path.join('inputs','%s_%s','gt3_%s_%s_neut.dat') % (str(shotid),(timeid), str(shotid),(timeid)))
        except:
            pass

    myPlasma = gt3(shotlabel=maxfile, mode=gt3Method, iolFlag = IOL, neutFlag = neutrals)


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

