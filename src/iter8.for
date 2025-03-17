      subroutine iter8 (errfct, x, xmin, xmax, ierr, trace)
C
C
      implicit double precision (a-h, o-z)

      integer i, ierr
      character(LEN=6) trace

      data tol, rtol /1.E-7,1.e-4/

      xlow = xmin
      xhigh = xmax
      ierr = 1

      if(xlow. gt. xhigh) stop ' Limits contradiction in ITER'

      if(xlow. eq. xhigh) then
        x = xlow
        ierr = 0
        return
      end if

      call errfct (xlow, fxlow, deriv)

      call errfct (xhigh, fxhigh, deriv)

      if (fxlow * fxhigh .gt. 0.) then

C     root is not within interval [xlow,xhigh]
C     return endpoint closest to the root

        if (abs(fxlow) .lt. abs(fxhigh)) then
          x = xlow
        else
          x = xhigh
        end if
        call errfct (x, fx, deriv)
        ierr = 0
        return
      else
        ierr = 1
      end if

      x1 = x

      call errfct (x1, fx, deriv)

C     does [xlow,x1] bracket the root?

      if (fx * fxlow .lt. 0.) xhigh = x1

C                                                                iteration loop
C------------------------------------------------------------------------------

      do i = 1, 50

        if (deriv .ne. 0. .and. i .le. 10) then
C                                                       compute newton estimate
          x2 = x1 - fx / deriv

          if (x2 .le. xlow .or. x2 .ge. xhigh) then
C                                                                        bisect
            x2 = (xlow + xhigh) / 2.

          end if

        else
C                                                                        bisect
          x2 = (xlow + xhigh) / 2.

        end if
C                                            get function value of new estimate
        ofx = fx
  10    call errfct (x2, fx, deriv)
C
  900 format(' itr:', i2,3g12.4)

        test1 = abs(fx)
C                                             scale the test criteria if x2 < 1
        if (abs(x2) .gt. 1000.*tol)   test1 = abs (fx / x2)
C!!  modify denominator in case x2 is zero:  12/02
        xdenom = max(abs(x1),abs(x2))
        if(xdenom .le. 1.e-9) then
C          write(77,'(a6,i2,6g12.3)') trace,i,x,x1,x2,xmin,xmax
C          stop ' error stop in iter '
          test2 = 2.*rtol
        else
          test2 = abs ((x2 - x1) / xdenom)
        end if
        dfx = abs(ofx - fx) 
C                                                              convergence test
        if (test1 .lt. tol .or. test2 .lt. rtol .or. (dfx .lt. tol
     &  .and. abs(x2) .lt. tol)) then

          x = x2
C          ierr = 0
          ierr = -i
          return

        end if
C                                                              bracket the root
        if (fx * fxhigh .gt. 0.) then
          xhigh = x2
        else
          xlow = x2
        end if

        x1 = x2

      end do
C                                                            end iteration loop
C------------------------------------------------------------------------------

      return

      end
