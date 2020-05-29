
      real function erf(x)

c-----This function returns the error function ERF(x) with fractional
c-----error everywhere less than 1.2*10-7.
c-----Adapted from the book "Numerical Recipes"

      real x, z, t, erfcc
      if (x.lt.9.0) then
         z=abs(x)      
         t=1./(1.+0.5*z)
         erfcc=t*exp(-z*z-1.26551223+t*(1.00002368+t*(.37409196+
     *     t*(.09678418+t*(-.18628806+t*(.27886807+t*(-1.13520398+
     *     t*(1.48851587+t*(-.82215223+t*.17087277)))))))))
         if (x.lt.0.) erfcc=2.-erfcc
      endif
      if (x.ge.9.0) erfcc = 0.0
      erf = 1.0 - erfcc
      
      return
      end


