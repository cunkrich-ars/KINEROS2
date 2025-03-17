C------------------------------------------------------------------------------

      entry sedbf_smith (qup, a2, q2, topw)

C                                       baseflow: compute initial concentration
C------------------------------------------------------------------------------

C     Upper boundary condition

      if( a2(1,1) .gt. 1.E-8) then
        do n = 1, nps
          sd(n,1) = 0.
        end do
        do j = 1, nu
           do n = 1, nps
             call get (upc(j,n,1), 1, cj)
             sd(n,1) = sd(n,1) + qup(j,1) * cj
           end do
        end do
        cs(n,1,1) = sd(n,1) / q2(1,1)
      else
        do n = 1, nps
          cs(n,1,1) = 0.
        end do
      end if

      cl(n,1,1) = cs(n,1,1)

C     Transport capacity

      do i = 1, nk
        if (a2(i,1) .gt. 1.E-8) then
          do n = 1, nps
            v = q2(i,1) / a2(i,1)
            d = a2(i,1) / topw(i,1)
            cml(n,i,1) = trnscap(d, v, slope(i,1), dpd(n), rho(n), grav)
          end do
        else
          cml(n,i,1) = 0.
        end if
      end do

!     Loop until concentrations stabilize

711   continue

C     Loop through spatial nodes...

      do i = 2, nk

        im = i - 1
        sumdel = 0.

        dcpx = cml(nlp(1),i,1) - cl(nlp(1),i,1)

        if (a2(i,1) .gt. 1.E-8) then
          v = q2(i,1) / a2(i,1)
          d = a2(i,1) / topw(i,1)
        else
          v = 0.
          d = 0.
        end if

C       Loop thru particle sizes...

        do n = 1, nps

          if (a2(i,1) .gt. 1.E-8) then

            cp = cl(n,i,1)
            cpm = cl(n,im,1)
            pro = 0.
            upav = 1.

            if (dcpx .gt. 0. .and. n .ge. nlp(1)) then
              pro = pros(n,1) * bcoh(1)
              if (deprt(n,i,1) .le. 1.E-6)then
                upav = 1. - pav(1)
                if (pav(1) .ge. 0.999) upav = 0.
              else if (pav(1) .ge. .999 .and.
     &                 depnet(i,1) .gt. 1.E-5) then
                prt = deprt(n,i,1) / depnet(i,1)
                pro = amin1 (prt, 1.)
                upav = 1.
              end if
            else
              if (cml(n,i,1) .lt. cp) pro = 1.
            end if

C           Compute erosion or deposition term

            cgu = vsl(n) * pro
            dfac = upav * cgu * topw(i,1)
            eterm = dfac * cml(n,i,1)

            rpt = .false.

715         dtrc = eterm

C           Compute concentration

            supc = dtrc
            af = cp * a2(i,1) / delt
            bf = (omega * cs(n,im,1) * q2(im,1) - comega *
     &           (cp * q2(i,1) - cpm * q2(im,1))) / dx
            df = 1. / delt + omega * v / dx
            deno = a2(i,1) * df + dfac
            cs(n,i,1) = (supc + af + bf) / deno

            if (cs(n,i,1) .lt. 0.) cs(n,i,1) = 0.

C           Check for mixture limit - limit erosion to proportion of least erodable class

            eros = eterm - cs(n,i,1) * dfac

            if (pav(1) .lt. 1. .and. n .gt. nlp(1)
     &          .and. .not. rpt .and. erolim .gt. 0.) then
              erop = erolim * pros(n,1) / pros(nlp(1),1)
              if (eros .gt. erop .and. eros .gt. 0.) then
                eterm = erop
                rpt = .true.
                go to 715
              end if
            end if

C           Erosion / deposition depth accounting

            if (n .eq. nlp(1)) erolim = eros

            dela = eros * delt * rho(n)
            deprt(n,i,1) = deprt(n,i,1) - dela
            sumdel = sumdel + dela

          else

            cs(n,i,1) = 0.

          end if

        end do
          depnet(i,1) = depnet(i,1) - sumdel
      end do

      rpt = .false.

      do i = 1, nk
        do n = 1, nps
          if (.not. rpt) then
            if (cs(n,i,1) .gt. 0.) then
              test = (cl(n,i,1) - cs(n,i,1)) / cs(n,i,1)
              if (abs(test) .gt. 0.001) rpt = .true.
            end if
          end if
          cl(n,i,1) = cs(n,i,1)
        end do
      end do

      if (sed_debug) then
C        do i = 1, nk
C          conc(i) = 0.
C          do n = 1, nps
C            conc(i) = conc(i) + cml(n,i,1)
C          end do
C        end do
C        write(88,'(15F10.6)') (conc(i), i = 1, nk)
        do i = 1, nk
          conc(i) = 0.
          do n = 1, nps
            conc(i) = conc(i) + cs(n,i,1)
          end do
        end do
        write(88,'(15F10.6)') (conc(i), i = 1, nk)
        write(88, '(1x)')
      end if

      if ( rpt ) go to 711

      do i = 1, nk
        do n = 1, nps
          cl(n,i,1) = cs(n,i,1)
        end do
      end do

      do n = 1, nps
        call store (outc(n,1), 1, cs(n,nk,1))
      end do

      return

