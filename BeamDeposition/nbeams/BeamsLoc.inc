c/    File BeamsLoc.inc with variables local to the NB module
c/    Created by John Mandrekas for the stand-alone NB module,
c/    07-07-98
c/
c/    NB deposition related variables:
c/ 
      real ab
      real pitchangl(mxrho,3,maxBeams), fbpwr(3,maxBeams)

      common /nbvars/ ab, pitchangl, fbpwr

c/    Variables related to the Gaussian integration:                            
                                                        
      real wr, ws, wz, sr(nr), ss(ns), sz(nz)                                
      common /gauss/ wr, ws, wz, sr, ss, sz                         
                                                                                
c/    Miscellaneous variables related to the RZ integration: 
                   
      integer iflux, ibeam, maxit, nstart                                              
      real sgn, rb, zb, rflux
      real  alongr(3), alongz(3), alongksi(3), avksi(3), capd0(3),
     .      capd1(3)                  
      common /integr/ iflux, ibeam, maxit, nstart, sgn, rb, zb, alongr,
     .    alongz, rflux, alongksi, avksi, capd0, capd1                            
                                                                                
c/    Local geometry and MHD variables                                          
      real rhoi, rhoSurf, elongr, elonga, elong0, delong, shft, dshft                    
     .                                                           
      common /local/  rhoi, rhoSurf, elongr, elonga, elong0, delong, 
     .   shft, dshft                               
                                                                                
c/    Cross section information:                                                
      real omfp(mxrho,3,maxBeams), denz(4)                                                              
      common /boley/ omfp, denz                                                       
                                                                                
c/    Beamline related private variables:                                                   
c/    These variables are the local equivalents of the global                   
c/    beamline parameters.                                                      
                                                                            
      integer shape, ptype
      real ebkev, rtan, rbeam, width, height, gaussR, gaussZ,
     .     zpos                                                                 
      common /privbeam/ ebkev, rtan, rbeam, width, height, gaussR,
     .     gaussZ, zpos, shape, ptype                                                 

c/    Current drive and fast ion related private variables: 
      
      real fNBion(mxrho,3,maxBeams), fNBelec(mxrho,3,maxBeams), 
     .   fDDHe3ion(mxrho), fDDHe3elec(mxrho), fDDpion(mxrho), 
     .   fDDpelec(mxrho), fDDtion(mxrho), fDDtelec(mxrho),
     .   fDTHe4ion(mxrho), fDTHe4elec(mxrho) 
      
      common /nbfast/ fNBion, fNBelec, fDDHe3ion, fDDHe3elec, fDDpion, 
     .   fDDpelec, fDDtion, fDDtelec, fDTHe4ion, fDTHe4elec
