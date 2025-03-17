      subroutine fmt10 (var, ref, string, n)


      integer n, j

      real :: var, valu, mult, ref

      character string*10


      valu = var

C                                          compute # of decimal places based on
C                                            7 sig. fig.'s from reference value
      if (ref .ge. 1.0) then

        j = 6 - int (alog10 (ref))

      else

        j = 7

      end if

      if (valu .ge. 1.E9 .or. -valu .ge. 1.E8) then
C                                                      use exponential notation
        write (string, 185) valu

      else if (j .ge. 7) then

        write (string, 105) valu

      else if (j .eq. 6) then

        write (string, 115) valu

      else if (j .eq. 5) then

        write (string, 125) valu

      else if (j .eq. 4) then

        write (string, 135) valu

      else if (j .eq. 3) then

        write (string, 145) valu

      else if (j .eq. 2) then

        write (string, 155) valu

      else if (j .eq. 1) then

        write (string, 165) valu

      else if (j .le. 0) then

        mult = 1.

        if (var .ge. 1.E8) then
C                                                  remove nonsignificant digits
          mult = 100.

        else if (var .ge. 1.E7) then

          mult = 10.

        end if

        valu = mult * float (int (var / mult))

        write (string, 175) valu

      end if

      write (string, fmt) valu
C                                     determine position of first nonblank char
      do n = 1, 10

        if (string(n:n) .ne. ' ') go to 20

      end do

20    return

105   format (f10.7)
115   format (f10.6)
125   format (f10.5)
135   format (f10.4)
145   format (f10.3)
155   format (f10.2)
165   format (f10.1)
175   format (f10.0)
185   format (e10.4)

      end
