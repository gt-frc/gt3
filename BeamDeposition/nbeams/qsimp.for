      subroutine qsimp(func,a,b,s)

c-----This subroutine returns as -s- the integral of the function -func-
c-----from -a- to -b-. The variable -eps- determines the desired fractional
c-----accuracy, and the maximum number of steps is 2**(jmax-1). Integration
c-----is performed by Simpson's rule. The routine has been taken from the 
c-----book "Numerical recipes" , by W. Press et al.
c-----Adapted for the SuperCode, John Mandrekas, GIT, 11-15-91

      implicit none

      external func

      integer j, jmax
      real a, b, eps, func, os, ost, s, st
      data eps/1.e-4/,jmax/100/

      ost=-1.e30
      os= -1.e30
      do  j = 1, jmax
         call trapzd(func,a,b,st,j)
         s=(4.*st-ost)/3.0
         if(abs(s-os).lt.eps*abs(os)) return
         os=s
         ost=st
      enddo

      return
      end


      subroutine trapzd(func,a,b,s,n)

      implicit none
      integer it, j, n
      real a, b, func, del, sum, s, tnm, x

      if (n.eq.1) then
         s=0.5*(b-a)*(func(a)+func(b))
         it=1
      else
         tnm=it
         del=(b-a)/tnm
         x=a+0.5*del
         sum=0.0
         do  j = 1, it
            sum=sum+func(x)
            x=x+del
         enddo
         s=0.5*(s+(b-a)*sum/tnm)
         it=2*it
      endif

      return
      end
