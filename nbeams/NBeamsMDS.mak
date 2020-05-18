# Microsoft Developer Studio Generated NMAKE File, Format Version 4.00
# ** DO NOT EDIT **

# TARGTYPE "Win32 (x86) Console Application" 0x0103

!IF "$(CFG)" == ""
CFG=NBeamsMDS - Win32 Debug
!MESSAGE No configuration specified.  Defaulting to NBeamsMDS - Win32 Debug.
!ENDIF 

!IF "$(CFG)" != "NBeamsMDS - Win32 Release" && "$(CFG)" !=\
 "NBeamsMDS - Win32 Debug"
!MESSAGE Invalid configuration "$(CFG)" specified.
!MESSAGE You can specify a configuration when running NMAKE on this makefile
!MESSAGE by defining the macro CFG on the command line.  For example:
!MESSAGE 
!MESSAGE NMAKE /f "NBeamsMDS.mak" CFG="NBeamsMDS - Win32 Debug"
!MESSAGE 
!MESSAGE Possible choices for configuration are:
!MESSAGE 
!MESSAGE "NBeamsMDS - Win32 Release" (based on\
 "Win32 (x86) Console Application")
!MESSAGE "NBeamsMDS - Win32 Debug" (based on "Win32 (x86) Console Application")
!MESSAGE 
!ERROR An invalid configuration is specified.
!ENDIF 

!IF "$(OS)" == "Windows_NT"
NULL=
!ELSE 
NULL=nul
!ENDIF 
################################################################################
# Begin Project
# PROP Target_Last_Scanned "NBeamsMDS - Win32 Debug"
F90=fl32.exe
RSC=rc.exe

!IF  "$(CFG)" == "NBeamsMDS - Win32 Release"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 0
# PROP BASE Output_Dir "Release"
# PROP BASE Intermediate_Dir "Release"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 0
# PROP Output_Dir "Release"
# PROP Intermediate_Dir "Release"
# PROP Target_Dir ""
OUTDIR=.\Release
INTDIR=.\Release

ALL : "$(OUTDIR)\NBeamsMDS.exe"

CLEAN : 
	-@erase ".\Release\NBeamsMDS.exe"
	-@erase ".\Release\erf.obj"
	-@erase ".\Release\zinteg.obj"
	-@erase ".\Release\nbdriver.obj"
	-@erase ".\Release\qsimp.obj"
	-@erase ".\Release\hunt.obj"
	-@erase ".\Release\sigfit.obj"
	-@erase ".\Release\calcbeams.obj"
	-@erase ".\Release\rinteg.obj"
	-@erase ".\Release\ssum.obj"
	-@erase ".\Release\fastions.obj"
	-@erase ".\Release\hofr.obj"
	-@erase ".\Release\beamconsts.obj"
	-@erase ".\Release\sinteg.obj"
	-@erase ".\Release\eisplit.obj"
	-@erase ".\Release\sig_olson.obj"
	-@erase ".\Release\getrho.obj"
	-@erase ".\Release\frate.obj"
	-@erase ".\Release\gausswts.obj"
	-@erase ".\Release\beams_mfp.obj"
	-@erase ".\Release\coulomb.obj"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

# ADD BASE F90 /Ox /I "Release/" /c /nologo
# ADD F90 /Ox /I "Release/" /c /nologo
F90_PROJ=/Ox /I "Release/" /c /nologo /Fo"Release/" 
F90_OBJS=.\Release/
# ADD BASE RSC /l 0x409 /d "NDEBUG"
# ADD RSC /l 0x409 /d "NDEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
BSC32_FLAGS=/nologo /o"$(OUTDIR)/NBeamsMDS.bsc" 
BSC32_SBRS=
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /subsystem:console /machine:I386
# ADD LINK32 kernel32.lib /nologo /subsystem:console /machine:I386
LINK32_FLAGS=kernel32.lib /nologo /subsystem:console /incremental:no\
 /pdb:"$(OUTDIR)/NBeamsMDS.pdb" /machine:I386 /out:"$(OUTDIR)/NBeamsMDS.exe" 
LINK32_OBJS= \
	"$(INTDIR)/erf.obj" \
	"$(INTDIR)/zinteg.obj" \
	"$(INTDIR)/nbdriver.obj" \
	"$(INTDIR)/qsimp.obj" \
	"$(INTDIR)/hunt.obj" \
	"$(INTDIR)/sigfit.obj" \
	"$(INTDIR)/calcbeams.obj" \
	"$(INTDIR)/rinteg.obj" \
	"$(INTDIR)/ssum.obj" \
	"$(INTDIR)/fastions.obj" \
	"$(INTDIR)/hofr.obj" \
	"$(INTDIR)/beamconsts.obj" \
	"$(INTDIR)/sinteg.obj" \
	"$(INTDIR)/eisplit.obj" \
	"$(INTDIR)/sig_olson.obj" \
	"$(INTDIR)/getrho.obj" \
	"$(INTDIR)/frate.obj" \
	"$(INTDIR)/gausswts.obj" \
	"$(INTDIR)/beams_mfp.obj" \
	"$(INTDIR)/coulomb.obj"

"$(OUTDIR)\NBeamsMDS.exe" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

!ELSEIF  "$(CFG)" == "NBeamsMDS - Win32 Debug"

# PROP BASE Use_MFC 0
# PROP BASE Use_Debug_Libraries 1
# PROP BASE Output_Dir "Debug"
# PROP BASE Intermediate_Dir "Debug"
# PROP BASE Target_Dir ""
# PROP Use_MFC 0
# PROP Use_Debug_Libraries 1
# PROP Output_Dir "Debug"
# PROP Intermediate_Dir "Debug"
# PROP Target_Dir ""
OUTDIR=.\Debug
INTDIR=.\Debug

ALL : "$(OUTDIR)\NBeamsMDS.exe"

CLEAN : 
	-@erase ".\Debug\NBeamsMDS.exe"
	-@erase ".\Debug\beams_mfp.obj"
	-@erase ".\Debug\gausswts.obj"
	-@erase ".\Debug\erf.obj"
	-@erase ".\Debug\sinteg.obj"
	-@erase ".\Debug\nbdriver.obj"
	-@erase ".\Debug\getrho.obj"
	-@erase ".\Debug\qsimp.obj"
	-@erase ".\Debug\frate.obj"
	-@erase ".\Debug\calcbeams.obj"
	-@erase ".\Debug\sigfit.obj"
	-@erase ".\Debug\eisplit.obj"
	-@erase ".\Debug\fastions.obj"
	-@erase ".\Debug\hunt.obj"
	-@erase ".\Debug\beamconsts.obj"
	-@erase ".\Debug\zinteg.obj"
	-@erase ".\Debug\sig_olson.obj"
	-@erase ".\Debug\ssum.obj"
	-@erase ".\Debug\hofr.obj"
	-@erase ".\Debug\coulomb.obj"
	-@erase ".\Debug\rinteg.obj"
	-@erase ".\Debug\NBeamsMDS.ilk"
	-@erase ".\Debug\NBeamsMDS.pdb"

"$(OUTDIR)" :
    if not exist "$(OUTDIR)/$(NULL)" mkdir "$(OUTDIR)"

# ADD BASE F90 /Zi /I "Debug/" /c /nologo
# ADD F90 /Zi /I "Debug/" /c /nologo
F90_PROJ=/Zi /I "Debug/" /c /nologo /Fo"Debug/" /Fd"Debug/NBeamsMDS.pdb" 
F90_OBJS=.\Debug/
# ADD BASE RSC /l 0x409 /d "_DEBUG"
# ADD RSC /l 0x409 /d "_DEBUG"
BSC32=bscmake.exe
# ADD BASE BSC32 /nologo
# ADD BSC32 /nologo
BSC32_FLAGS=/nologo /o"$(OUTDIR)/NBeamsMDS.bsc" 
BSC32_SBRS=
LINK32=link.exe
# ADD BASE LINK32 kernel32.lib /nologo /subsystem:console /debug /machine:I386
# ADD LINK32 kernel32.lib /nologo /subsystem:console /debug /machine:I386
LINK32_FLAGS=kernel32.lib /nologo /subsystem:console /incremental:yes\
 /pdb:"$(OUTDIR)/NBeamsMDS.pdb" /debug /machine:I386\
 /out:"$(OUTDIR)/NBeamsMDS.exe" 
LINK32_OBJS= \
	"$(INTDIR)/beams_mfp.obj" \
	"$(INTDIR)/gausswts.obj" \
	"$(INTDIR)/erf.obj" \
	"$(INTDIR)/sinteg.obj" \
	"$(INTDIR)/nbdriver.obj" \
	"$(INTDIR)/getrho.obj" \
	"$(INTDIR)/qsimp.obj" \
	"$(INTDIR)/frate.obj" \
	"$(INTDIR)/calcbeams.obj" \
	"$(INTDIR)/sigfit.obj" \
	"$(INTDIR)/eisplit.obj" \
	"$(INTDIR)/fastions.obj" \
	"$(INTDIR)/hunt.obj" \
	"$(INTDIR)/beamconsts.obj" \
	"$(INTDIR)/zinteg.obj" \
	"$(INTDIR)/sig_olson.obj" \
	"$(INTDIR)/ssum.obj" \
	"$(INTDIR)/hofr.obj" \
	"$(INTDIR)/coulomb.obj" \
	"$(INTDIR)/rinteg.obj"

"$(OUTDIR)\NBeamsMDS.exe" : "$(OUTDIR)" $(DEF_FILE) $(LINK32_OBJS)
    $(LINK32) @<<
  $(LINK32_FLAGS) $(LINK32_OBJS)
<<

!ENDIF 

.for{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

.f{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

.f90{$(F90_OBJS)}.obj:
   $(F90) $(F90_PROJ) $<  

################################################################################
# Begin Target

# Name "NBeamsMDS - Win32 Release"
# Name "NBeamsMDS - Win32 Debug"

!IF  "$(CFG)" == "NBeamsMDS - Win32 Release"

!ELSEIF  "$(CFG)" == "NBeamsMDS - Win32 Debug"

!ENDIF 

################################################################################
# Begin Source File

SOURCE=.\zinteg.for
DEP_F90_ZINTE=\
	".\nbparams.inc"\
	".\BeamsLoc.inc"\
	".\nbconsts.inc"\
	

"$(INTDIR)\zinteg.obj" : $(SOURCE) $(DEP_F90_ZINTE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\ssum.for

"$(INTDIR)\ssum.obj" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\sinteg.for
DEP_F90_SINTE=\
	".\nbparams.inc"\
	".\BeamsLoc.inc"\
	".\nbplasma.inc"\
	

"$(INTDIR)\sinteg.obj" : $(SOURCE) $(DEP_F90_SINTE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\sigfit.for

"$(INTDIR)\sigfit.obj" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\sig_olson.for

"$(INTDIR)\sig_olson.obj" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\rinteg.for
DEP_F90_RINTE=\
	".\nbparams.inc"\
	".\BeamsLoc.inc"\
	".\nbplasma.inc"\
	

"$(INTDIR)\rinteg.obj" : $(SOURCE) $(DEP_F90_RINTE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\qsimp.for

"$(INTDIR)\qsimp.obj" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\nbdriver.for

"$(INTDIR)\nbdriver.obj" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\hunt.for

"$(INTDIR)\hunt.obj" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\hofr.for
DEP_F90_HOFR_=\
	".\nbparams.inc"\
	".\nbplasma.inc"\
	".\BeamsLoc.inc"\
	".\nbconsts.inc"\
	

"$(INTDIR)\hofr.obj" : $(SOURCE) $(DEP_F90_HOFR_) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\getrho.for
DEP_F90_GETRH=\
	".\nbparams.inc"\
	".\nbplasma.inc"\
	".\BeamsLoc.inc"\
	

"$(INTDIR)\getrho.obj" : $(SOURCE) $(DEP_F90_GETRH) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\gausswts.for

"$(INTDIR)\gausswts.obj" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\frate.for

"$(INTDIR)\frate.obj" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\fastions.for
DEP_F90_FASTI=\
	".\nbparams.inc"\
	".\BeamsLoc.inc"\
	".\nbplasma.inc"\
	".\nbconsts.inc"\
	

"$(INTDIR)\fastions.obj" : $(SOURCE) $(DEP_F90_FASTI) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\erf.for

"$(INTDIR)\erf.obj" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\eisplit.for
DEP_F90_EISPL=\
	".\nbparams.inc"\
	".\nbconsts.inc"\
	".\nbplasma.inc"\
	

"$(INTDIR)\eisplit.obj" : $(SOURCE) $(DEP_F90_EISPL) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\coulomb.for

"$(INTDIR)\coulomb.obj" : $(SOURCE) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\calcbeams.for
DEP_F90_CALCB=\
	".\nbparams.inc"\
	".\BeamsLoc.inc"\
	".\nbplasma.inc"\
	".\nbconsts.inc"\
	

"$(INTDIR)\calcbeams.obj" : $(SOURCE) $(DEP_F90_CALCB) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\beams_mfp.for
DEP_F90_BEAMS=\
	".\nbparams.inc"\
	".\BeamsLoc.inc"\
	".\nbplasma.inc"\
	

"$(INTDIR)\beams_mfp.obj" : $(SOURCE) $(DEP_F90_BEAMS) "$(INTDIR)"


# End Source File
################################################################################
# Begin Source File

SOURCE=.\beamconsts.for
DEP_F90_BEAMC=\
	".\nbconsts.inc"\
	

"$(INTDIR)\beamconsts.obj" : $(SOURCE) $(DEP_F90_BEAMC) "$(INTDIR)"


# End Source File
# End Target
# End Project
################################################################################
