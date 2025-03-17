C   Description:

C   Combines water and sediment from up to 10 upstream elements

C-------------------------------------------------------------------------------

C  Global parameters:

C     sed            log           .true. = route sediment (module runpars)

C     nps            int            number of particle classes (module runpars)

C     limit          int            number of time steps (module runpars)

C     typname        char          'ADDER' (module elpars)

C-------------------------------------------------------------------------------

C  Input Block (input file labels in parenthesis):

C     id (ID)        int*4         identifier (up to 6 digits)

C     idup (UP)      int*4         identifiers of up to 10 upstream
C                                  contributing elements

C-------------------------------------------------------------------------------

C  Subroutines/functions:

C     getr4         (reader.for)
C     getstr             "
C     geti4              "
C     new           (miscel.for)
C     old                "
C     get                "
C     store              "
C     geta               "
C     stora              "
C     errxit        (errxit.for)
C-------------------------------------------------------------------------------

      subroutine adder ()

      use runpars
      use elpars

      implicit none

      integer ierr, i, j, ms, n, nu, outr, outq(2), nr, nq, nc 
      integer upq(10), idup(10), upr(10), upc(10,5), outc(5)

      real qrain, qrmax, q1, q2, vo, qru, qu, cu, area, coarea, wsi
      real qs_sum(5), qp, tp, c1(5), c2(5), wso(5), twso, qsp, tsp, dso

      character msg*18

      logical first

      character(LEN=2),dimension(5,2) :: attr = reshape((/'S1', 'S2', 
     & 'S3', 'S4', 'S5', 's1', 's2', 's3', 's4', 's5'/),(/5,2/))


      do j = 1, 15
        qbal(j) = 0.
      end do

      ltyp = 8

      call geti4 ('ID', 0, id, ierr)

      if (ierr .eq. 3) then

        call errxit ('ADDER', 'invalid identifier (ID)')

      else if (ierr .ge. 1) then

        call errxit ('ADDER', 'missing identifier (ID)')

      end if

      msg(1:16) = typname(8)
      write (msg(15:18), '(i4)') id

      do ms = 1, 5
        if (msg(ms:ms) .ne. ' ') go to 1
      end do

1     continue

      nu = 0

      do j = 1, 10

        call geti4 ('UP', j, idup(j), ierr)

        if (ierr .eq. 3) call errxit
     &  (msg(ms:18), 'invalid upstream identifier (UP)')

        if (ierr .eq. 0) then

          nu = nu + 1

        else

          go to 11

        end if
      end do

11    continue

      coarea = 0.

      do j = 1, nu

        upr(j) = 0

        call old (idup(j), 'RA', upr(j), nr, ierr)

        call old (idup(j), 'RO', upq(j), nq, ierr)

        if (ierr .gt. 0) call errxit
     &  (msg(ms:18), 'upstream inflow not found')

        call geta (idup(j), area)
        coarea = coarea + area

        if (sed) then

          do n = 1, nps

            call old (idup(j), attr(n,1), upc(j,n), nc, ierr)

            if (ierr .gt. 0) call errxit
     &      (msg(ms:18), 'upstream sediment concentration not found')

          end do

        end if

      end do

      call stora (id, coarea)

      call new (id, 'RA', outr)
      call new (id, 'RO', outq(1))

      if (sed) then

        do n = 1, nps

          call new (id, attr(n,1), outc(n))

          wso(n) = 0.

        end do

      end if

      vo = 0.
      qp = 0.
      tp = 0.
      qsp = 0.
      tsp = 0.
      wsi = 0.
      qrmax = 0.
      first = .true.

      do i = 1, limit - 1

        qrain = 0.
        q2 = 0.

        if (sed) then
          do n = 1, nps
            qs_sum(n) = 0.
          end do
        end if

        do j = 1, nu

          if (upr(j) .gt. 0) then

            call get (upr(j), i, qru)

            qrain = qrain + qru

          end if

          call get (upq(j), i, qu)

          q2 = q2 + qu

          if (sed) then

            do n = 1, nps

              call get (upc(j,n), i, cu)

              qs_sum(n) = qs_sum(n) + cu * qu ! sed. discharge

            end do

          end if

        end do

        call store (outr, i, qrain)

        if (qrain .gt. qrmax) qrmax = qrain

        call store (outq(1), i, q2)

        if (sed) then

          twso = 0.

          do n = 1, nps

            if (q2 .gt. 0.) then
              c2(n) = qs_sum(n)/q2
            else
              c2(n) = 0.
            end if

            call store (outc(n), i, c2(n))

            if (.not. first) then
              dso = (q1 * c1(n) + q2 * c2(n)) * .5 * rho(n)
              twso = twso + dso
              wso(n) = wso(n) + dso * delt
            end if

            c1(n) = c2(n)

          end do

          if (twso .gt. qsp) then
            qsp = twso
            tsp = float (i - 1) * delt
          end if

        end if

        if (first) then
          first = .false.
        else
          vo = vo + delt * 0.5 * (q1 + q2)
        end if

        if (q2 .gt. qp) then
          qp = q2
          tp = float (i - 1) * delt
        end if

        q1 = q2

      end do

      qbal(10) = vo

      call qwrt (id, 1, msg, qp, tp, vo, coarea, outq, 8, outr, qrmax)

      if (sed) then

        do n = 1, nps
          wsi = wsi + wso(n)
        end do

        call swrt ( wso, wsi, 0., 0., qsp, tsp, outc, 0., 0.)

      end if

      rsoil = .false.

      return

      end
