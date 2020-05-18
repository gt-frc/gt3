c///////////////////////////////////////////////////////////////////////
c/ 
c/ NAME
c/    hunt
c/
c/ DESCRIPTION
c/    Routine for search with correlated values:
c/    Given an array xx of length n, and given a value x, this
c/    subroutine returns a value jlo, such that x is between
c/    xx(jlo) and xx(jlo+1). The array xx must be monotonic, either
c/    increasing or decreasing. On input, jlo is taken as the initial
c/    guess for jlo on output.
c/    Written by John Mandrekas, GIT, 02-02-92, for the SuperCode.
c/    Adapted from "Numerical Recipes", by W. Press, et al., pg 91
c/
c/ PARAMETERS
c/    xx: array of length n
c/    n : lenth of array
c/    x : value of x whose position we seek
c/
c/ RETURNS
c/    jlo: on output, x is between xx(jlo) and xx(jlo+1)
c/
c/ NOTE:
c/    This is a "hidden" routine, i.e., it does not appear in Beams.mod
c/
c///////////////////////////////////////////////////////////////////////

      subroutine hunt(xx, n, x, jlo)

      integer inc, jhi, jlo, jm, n
      real xx(n), x
      logical ascnd

      ascnd = xx(n).gt.xx(1)
      if (jlo.le.0.or.jlo.gt.n) then
        jlo = 0
        jhi = n + 1
        go to 3
      endif
      inc = 1
      if(x.ge.xx(jlo).eqv.ascnd) then    ! Hunt up
1       jhi = jlo + inc
        if (jhi.gt.n) then
          jhi = n + 1
        else if(x.ge.xx(jhi).eqv.ascnd) then
          jlo = jhi
          inc = inc + inc
          go to 1
        endif                   ! Done hunting, value bracketed
      else                      ! Hunt down
        jhi = jlo
2       jlo = jhi - inc
        if(jlo.lt.1) then
          jlo = 0
        else if (x.lt.xx(jlo).eqv.ascnd) then
          jhi = jlo
          inc = inc + inc
          go to 2
        endif                   ! Done hunting, value bracketed
      endif

c-----Begin final bisection phase:
3     if (jhi-jlo.eq.1) return
      jm = (jhi + jlo)/2
      if (x.gt.xx(jm).eqv.ascnd) then
        jlo = jm
      else
        jhi = jm
      endif
      go to 3 
      
      end

