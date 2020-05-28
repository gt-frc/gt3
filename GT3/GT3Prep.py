#!/usr/bin/env python2
# -*- coding: utf-8 -*-
import os
from Functions.GetVals import getVals
from Processors.WriteFile import writeFile

def gt3Prep(s, t, r, m, quiet=False, genFiles=False):
    """
    This function allows you to prep GT3 in an interactive manner.

    Your file structure should look something like this if you want to use this template, or you could create your own

├── inputs
│   ├── 118888_1525
│   ├── 118888_1570
│   │   ├── gt3_118888_1570_bpol.dat
│   │   ├── gt3_118888_1570_btor.dat
│   │   ├── gt3_118888_1570_er.dat
│   │   ├── gt3_118888_1570_exlni.dat
│   │   ├── gt3_118888_1570_exlte.dat
│   │   ├── gt3_118888_1570_exlti.dat
│   │   ├── gt3_118888_1570_fracz.dat
│   │   ├── gt3_118888_1570_fz1.dat
│   │   ├── gt3_118888_1570_nD.dat
│   │   ├── gt3_118888_1570_ne.dat
│   │   ├── gt3_118888_1570_neut.dat
│   │   ├── gt3_118888_1570_psirz.dat
│   │   ├── gt3_118888_1570_q.dat
│   │   ├── gt3_118888_1570_R.dat
│   │   ├── gt3_118888_1570_Te.dat
│   │   ├── gt3_118888_1570_Ti.dat
│   │   ├── gt3_118888_1570_vpolC.dat
│   │   ├── gt3_118888_1570_vpolD.dat
│   │   ├── gt3_118888_1570_vtorC.dat
│   │   ├── gt3_118888_1570_vtorD.dat
│   │   ├── gt3_118888_1570_zbar2.dat
│   │   ├── gt3_118888_1570_Z.dat
│   │   └── gt3_diiid_wall.dat
│   ├── 144977_3000
│   ├── 164436_3720
│   ├── 164988_1915
│   ├── 175826_2010
│   ├── 92980_3000
│   ├── 92980_3600
│   ├── togt3_d3d_118888_1570
│   ├── togt3_d3d_118890_1560
│   ├── togt3_d3d_164436_3720



    :param s: shot id
    :param t: time id
    :param r: run id (unsure if necessary)
    :param m: filename for input file to be created
    :param quiet: If True, does not provide interactive requests. Will use input file that it gt3 finds. If no data file,
                  system will halt
    :param genFiles: Flag to tell code in the future to transform DIII-D data into format we want for gt3
    :return: Your mom
    """

    varDict = {}

    #    direct="MaxPlasma\\inputs\%s\%s" % (str(s),str(t))
    direct = os.path.join("inputs", "%s_%s") % (str(s), str(t))
    if not os.path.exists(direct):
        print """Directory for shot %s.%s does not exist at inputs/%s_%s
                 Directory used: %s
                 Shutting down now.""" % (str(s), str(t), str(s), str(t), str(os.getcwd() + direct))
        raise Exception
    fileName = m
    fpath = os.path.join("inputs", fileName)
    if not os.path.exists(os.path.join("inputs", fileName)):
        if quiet:
            print """Silent mode active: No input file found.

                     Shutting down now."""
            raise Exception
        data = getVals(s, t, fileName)
        writeFile(s, t, fpath, data, reNeut=False)
    else:
        print """Input file for shot %s.%s exists""" % (str(s), str(t))
        while True:
            if quiet:
                print """Silent mode active: Using current input file for shot %s.%s""" % (str(s), str(t))
                break
            ch = raw_input("Use current input file? (Y/N) ")
            if (str(ch) == 'Y' or str(ch) == 'y'):
                break
            elif (str(ch) == 'N' or str(ch) == 'n'):
                os.remove(os.path.join("inputs", fileName))
                data = getVals(s, t, fileName)
                writeFile(s, t, fpath, data, reNeut=False)
                break
            else:
                print "Invalid selection \n"

        #####################################################################################################
        #
        #   This section will probably be usable for when we do not have data ready for gt3
        #
        #   TODO: Implement this somehow
        #
        #####################################################################################################


    if genFiles:
        fileList = {
            'data.xne': ('gt3_%s_%s_ne.dat' % (str(s), str(t)), data.xne),
            'data.xni': ('gt3_%s_%s_ni.dat' % (str(s), str(t)), data.xni),
            'data.xte': ('gt3_%s_%s_Te.dat' % (str(s), str(t)), data.xte),
            'data.xti': ('gt3_%s_%s_Ti.dat' % (str(s), str(t)), data.xti),
            'data.xer': ('gt3_%s_%s_er.dat' % (str(s), str(t)), data.xer),
            'data.fz1': ('gt3_%s_%s_fz1.dat' % (str(s), str(t)), data.fz1),
            'data.fracz': ('gt3_%s_%s_fracz.dat' % (str(s), str(t)), data.fracz),
            'data.exlti': ('gt3_%s_%s_exlti.dat' % (str(s), str(t)), data.exlti),
            'data.exlte': ('gt3_%s_%s_exlte.dat' % (str(s), str(t)), data.exlte),
            'data.exlni': ('gt3_%s_%s_exlni.dat' % (str(s), str(t)), data.exlni),
            'data.vpolC': ('gt3_%s_%s_vpolC.dat' % (str(s), str(t)), data.vpolC),
            'data.vtorC': ('gt3_%s_%s_vtorC.dat' % (str(s), str(t)), data.vtorC),
            'data.q': ('gt3_%s_%s_q.dat' % (str(s), str(t)), data.q),
            'data.zbar2': ('gt3_%s_%s_zbar2.dat   ' % (str(s), str(t)), data.zbar2)}
        try:
            os.chdir("inputs")
        except:
            raise Exception("gt3 input directory not found")

        #####################################################################################################
        #
        #   This section will probably be usable for when we do not have data ready for gt3
        #
        #   TODO: Implement this somehow
        #
        #####################################################################################################

        unitCon = {'data.xti': 1E-3,
                   'data.xte': 1E-3,
                   'data.xer': 1E-3}
        for var in fileList.keys():
            if var in unitCon.keys():
                with open(fileList[var][0], "w") as openFile:
                    for a, b in zip(data.rhor, fileList[var][1]):
                        openFile.write("%s    %s \n" % (str(a), str(b * unitCon[var])))
            else:
                with open(fileList[var][0], "w") as openFile:
                    for a, b in zip(data.rhor, fileList[var][1]):
                        openFile.write("%s    %s \n" % (str(a), str(b)))

        if hasattr(data, 'vpolD'):
            with open('gt3_%s_%s_vpolD.dat' % (str(s), str(t)), "w") as openFile:
                for a, b in zip(data.rhor, data.vpolD):
                    openFile.write("%s   %s \n" % (a, b))
        else:
            with open('gt3_%s_%s_vpolD.dat' % (str(s), str(t)), "w") as openFile:
                for a, b in zip(data.rhor, [0.] * len(data.rhor)):
                    openFile.write("%s   %s \n" % (a, b))

        if hasattr(data, 'vtorD'):
            with open('gt3_%s_%s_vtorD.dat' % (str(s), str(t)), "w") as openFile:
                for a, b in zip(data.rhor, data.vtorD):
                    openFile.write("%s   %s \n" % (a, b))
        else:
            with open('gt3_%s_%s_vtorD.dat' % (str(s), str(t)), "w") as openFile:
                for a, b in zip(data.rhor, [0.] * len(data.rhor)):
                    openFile.write("%s   %s \n" % (a, b))
