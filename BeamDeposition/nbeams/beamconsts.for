      block data nbconsts

c/    This blockdata subprogram initializes several constants that are
c/    used in the NB module
      
      implicit none
      include 'nbconsts.inc'
      
      data pi /3.1415926536/, muZero /1.256637E-06/, 
     .     eCharge /1.6022e-19/
      
      end
