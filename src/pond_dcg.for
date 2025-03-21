C   Code type: fortran subroutine

C   Compiler: microsoft fortran v5.1

C   Programmed by: C. Unkrich

C   Date: 2/96

C   Description:

C     computes the volume and associated outflow of a detention pond with a
C     given volume-discharge relation in tabular form.  inputs can include
C     inflow from channels or pipes, lateral inflow from planes, and
C     rainfall. Infiltration is computed as a constant rate.  Sediment
C     coming in with upstream and lateral inflow is also routed through
C     the pond.  uniform mixing is assumed by this method.

C     For output, volumetric rainfall rate from upstream areas for each time
C     step are added and stored.  Dividing by contributing area gives the
C     integrated rainfall rate over the area contributing to the pond.
C------------------------------------------------------------------------------

C   Arguments:

C     qbal        real          array with global volume balance components;
C                               qbal(1)  = area                  (sq.m/sq.ft),
C                               qbal(2)  = rainfall              (cu.m/cu.ft),
C                               qbal(5)  = initial pond storage  (cu.m/cu.ft),
C                               qbal(7)  = infiltration          (cu.m/cu.ft),
C                               qbal(9)  = storage               (cu.m/cu.ft),
C                               qbal(10) = outflow               (cu.m/cu.ft),

C     delt          real        user-specified time step (sec),

C     limit         int         number of time steps,

C     units         int         1 = metric, 2 = english,

C     sed           log         .true. = route sediment,

C     diag          log         .true. = write diagnostic info to unit 99,

C------------------------------------------------------------------------------

C   Input block (input file labels in parenthesis):

C     id (id)        int*4         identifier (up to 6 digits),

C     idup (up)      int*4         identifiers of up to 10 upstream
C                                  contributing elements,

C     idlat (la)     int*4         identifier(s) of 1 or 2 lateral elements,

C     nd (n)         int*4         number of data pairs in v,q, and sa,

C     v (vol)        real          tabulated values of volume and
C     q(dis)         real          corresponding discharge and
C     sa(surf)       real          surface area [l**3],[l**3/t],[l**2],

C     voli (stor)    real          initial storage volume [l**3],
C
C     fmin (ks)      real          constant infiltration rate [l/t],
C------------------------------------------------------------------------------

C   Subroutines/functions:

C     erpond         (pond.for)
C     iter           (miscel.for)
C     progss         (K2RunW.for)
C------------------------------------------------------------------------------

C   Variables:

C     i         time index;
C     m         index to breakpoints in the rainfall pattern;
C     k         index to data in the outflow rating table;
C     qup       average inflow rate between time-dt, time [l**3/t];
C     cup       concentration of sediment in flow qup
C     qlat      average lateral inflow rate [l**3/t];
C     clat      concentration of sediment in flow qlat
C     v1,v2     volume at time-dt, time [l**3];
C     vthres    threshold volume when outflow begins [l**3];
C     vlog      logarithm of rating table volume entry;
C     q1,q2     discharge at time-dt, time [l**3/t];
C     c1,c2     sed. conc. at time-dt, time
C     qlog      logarithm of rating table discharge entry;
C     sa1,sa2   surface area at time-dt, time [l**2];
C     sath      surface area at outflow threshold [l**2];
C     samax     maximum surface area;
C     salog     logarithm of rating table surface area entry;
C     dt        variable time increment;
C     rf        rainfall intensity [l/t];
C     qp        peak flow during event [l**3/t];
C     tp        time to peak.
C     vi        total inflow volume [l**3];
C     vf        total infiltrated volume [l**3];
C     vs        volume above voli remaining in pond at tfin [l**3];
C     vr        total rainfall volume;
C     vo        total outflow volume [l**3].
C     wsi       total weight of incoming sediment
C     wso(n)    total weight of outgoing sediment class size n
C     wse       total weight of deposited sediment
C     wss       total weight of sediment suspended in vs
C------------------------------------------------------------------------------


      subroutine pond (t, z, np)

      use elpars
      use runpars

      common /pond1/ v(50), vlog(50), q(50), qlog(50), sa(50), dt, fmin,
     &               salog(50), qup, k, v1, q1, sa1, qlat, vthres, sath,
     &               nd, i, samax, qu1, qu2, ql1, ql2, rin

      integer ierr, i, j, m, k, nd, np, outq(2), outr, outc(5), !id,
     &        idup(10), idlat(2), upq(10,2), latq(2), upc(10,5,2),
     &        latc(2,5), upr(10), latr(2), ms, ngage !, nps, limit

      character msg*20, attr*2, source*6
      character(LEN=1) :: cr

      dimension attr(5,2)

      real qsin1(5), qsin2(5), c1(5), c2(5), cu1(5),cu2(5), cl1(5),
     &     cl2(5), clat1(5), clat2(5), cup1(5), cup2(5), 
     &     qu(10,2), ql(2), wso(5)
C
      real, dimension(5000) :: t, z, r
!     &      !, rho(5), qbal(10)

      external erpond

      data msg / '       (POND)'/

      data attr /'S1', 'S2', 'S3', 'S4', 'S5',
     &           's1', 's2', 's3', 's4', 's5'/

	VMIN = 0.

C                                                          initialize variables
C------------------------------------------------------------------------------
      q1 = 0.
      qup1 = 0.
      qu1 = 0.
      qlat1 = 0.
      ql1 = 0.
      vi=0.
      vf=0.
      vs=0.
      vo=0.
      tp = 0.
      qp = 0.
      ipeak =1
      coarea = 0.
      iexceed = 0

      r_cum = 0.

      do j = 1, 10
        qbal(j) = 0.
      end do
      ltyp = 3
      sumtab(ltab)%twola = .false.
C** 4/03
      sumtab(ltab)%thst = 0.
      if (sed) then

        do n = 1, nps
          qsin1(n) = 0.
          qsin2(n) = 0.
          c1(n)  = 0.
          c2(n)  = 0.
          cup1(n) = 0.
          cup2(n) = 0.
          cu1(n) = 0.
          cu2(n) = 0.
          cl1(n) = 0.
          cl2(n) = 0.
          clat1(n)=0.
          clat2(n)=0.
          wso(n) = 0.
        end do

        wss = 0.
        wsi = 0.
        wse = 0.
        qsp = 0.
        tsp = 0.

      end if
C                                                    get basic input parameters
C------------------------------------------------------------------------------

C                                                                    IDENTIFIER
      call geti4 ('ID', 0, id, ierr)

      if (ierr .eq. 3) then
        call errxit ('POND', 'invalid identifier (ID)')
      else if (ierr .ge. 1) then
        call errxit ('POND', 'missing identifier (ID)')
      end if

      msg(1:16) = typname(3)
      write (msg(15:20), '(i6)') id

      do ms = 1, 5
        if (msg(ms:ms) .ne. ' ') go to 1
      end do

1     continue

      OPEN( 1, FILE = 'pond_aux.csv', STATUS = 'UNKNOWN' )
      WRITE(1,'("Time, Vol., SA, Vol. infiltrated, rainfall depth, Q")')

      nu = 0
C                                                            upstream inflow(s)
      do j = 1, 10

        call geti4 ('UP', j, idup(j), ierr)

        if (ierr .eq. 3) call errxit
     &  (msg(ms:13), 'invalid upstream identifier (UP)')

        if (ierr .eq. 0) then
          nu = nu + 1
        else
          go to 11
        end if

      end do

11    nl = 0

      do j = 1, 2
C                                                             lateral inflow(s)
        call geti4 ('LA', j, idlat(j), ierr)

        if (ierr .eq. 3) call errxit
     &  (msg(ms:13), 'invalid lateral element identifier (LA)')

        if (ierr .eq. 0) then
          nl = nl + 1
        else
          go to 12
        end if

      end do

12    continue
C                                                          number of data pairs
      call geti4 ('N', 0, nd, ierr)

      if (ierr .gt. 0) call errxit
     &(msg(ms:13), 'Number of rating table entries (N) not specified')

      if (nd .gt. 49) call errxit
     &(msg(ms:13), 'Too many rating table entries')

C                                        volume, discharge, and surface area --
C                                  compute logs and determine outflow threshold
      do  j = 1, nd

        call getr4 ('V', j, v(j), ierr)

        if (ierr .gt. 0) call errxit
     &  (msg(ms:13), 'volume entry (V) missing or invalid')

        call getr4 ('D', j, q(j), ierr)

        if (ierr .gt. 0) call errxit
     &  (msg(ms:13), 'discharge entry (D) missing or invalid')

        call getr4 ('SU', j, sa(j), ierr)

        if (ierr .gt. 0) call errxit
     &  (msg(ms:13), 'surface area entry (SU) missing or invalid')

        if (v(j) .gt. 0.) vlog(j) = alog (v(j))

        if (q(j) .gt. 0.) then

          qlog(j) = alog (q(j))

          if (q(j-1) .le. 0.) then
            vthres = v(j-1)
            sath = sa(j-1)
          end if

        end if

        if (sa(j) .gt. 0.) salog(j) = alog (sa(j))

      end do
C                                                        initial storage volume
      call getr4 ('ST', 0, voli, ierr)

      if (ierr .gt. 0) voli = 0.

      qbal(5) = voli
C                                                             infiltration rate
      call getr4 ('K', 0, fmin, ierr)

      if (ierr .eq. 0) then

        if (units .eq. 1) then
          fmin = fmin / 1000.
        else
          fmin = fmin / 12.
        end if

        fmin = fmin / 3600.

      else

        fmin = 0.
      end if

      samax = sa(nd)
      VMAX = 1.1 * v(nd)
C                                              get indices to storage locations
C                                              for rainfall, inflow and outflow
C------------------------------------------------------------------------------
      call new (id, 'RA', outr)
C                                                              rainfall storage
      do j = 1, nu

        upr(j) = 0
C                                                    upstream rain contribution
        call old (idup(j), 'RA', upr(j), nr, ierr)

        call old (idup(j), 'RO', upq(j,1), nq, ierr)
C                                                               upstream inflow
        if (ierr .gt. 0) call errxit
     &  (msg(ms:13), 'upstream inflow not found')

        call old (idup(j), 'ro', upq(j,2), nq, ierr)
C                                                     compound channel upstream
        if (ierr .gt. 0) upq(j,2) = -1

        if (sed) then
C                                                                      sediment
          do n = 1, nps

            call old (idup(j), attr(n,1), upc(j,n,1), nc, ierr)

            if (ierr .gt. 0) call errxit
     &      (msg(ms:13), 'upstream sediment concentration not found')

            if (upq(j,1) .lt. 0) then

              call old (idup(j), attr(n,2), upc(j,n,2), nc, ierr)

              if (ierr .gt. 0) call errxit
     &        (msg, 'upstream sediment concentration not found')

            end if

          end do

        end if

      end do

      do j = 1, nl

        latr(j) = 0
C                                                     lateral rain contribution
        call old (idlat(j), 'RA', latr(j), nr, ierr)

        call old (idlat(j), 'RO', latq(j), nq, ierr)
C                                                                lateral inflow
        if (ierr .gt. 0) call errxit
     &  (msg(ms:13), 'lateral inflow not found')

        if (sed) then
C                                                                      sediment
          do n = 1, nps

            call old (idlat(j), attr(n,1), latc(j,n), nc, ierr)

            if (ierr .gt. 0) call errxit
     &      (msg(ms:13), 'lateral sediment concentration not found')

          end do

        end if

      end do
C                                                               outflow storage
      call new (id, 'RO', outq(1))

      if (sed) then
        do n = 1, nps
          call new (id, attr(n,1), outc(n))
        end do
      end if
C                                                      update contributing area
C------------------------------------------------------------------------------
      do j = 1, nu

        call geta (idup(j), area)
        coarea = coarea + area

      end do

      do j = 1, nl

        call geta (idlat(j), area)
        coarea = coarea + area

      end do

C                                                                      rainfall
C------------------------------------------------------------------------------
      if (pond_rain) then
        rarea = samax
        z(np+1) = z(np)
        t(np+1) = t(np) + delt
        do ir = 2, np + 1
          r(ir-1) = (z(ir) - z(ir-1)) / (t(ir) - t(ir-1))
        end do
      else
        rarea = 0.
        qbal(1) = 0. ! xleng * (width(1,1) + width(nk,1)) /2.
        np = 2
        t(1) = 0.
        r(1) = 0.
        t(2) = float (limit) * delt
        r(2) = 0.
      end if

      qbal(1) = rarea
      coarea = coarea + rarea

      call stora (id, coarea)
C                               compute integrated rain rate for each time step
C------------------------------------------------------------------------------
      qrmax = 0.
      k = 1
      qrf = r(1) * samax

      do i = 1, limit

        time = delt * float (i - 1)

        if (time .ge. t(k+1)) then
          k = k + 1
          qrf = r(k) * samax
        end if

        qrain = qrf

        do j = 1, nu
          if (upr(j) .gt. 0) then
            call get (upr(j), i, qr)
            qrain = qrain + qr
          end if
        end do

        do j = 1, nl
          if (latr(j) .gt. 0) then
            call get (latr(j), i, qr)
            qrain = qrain + qr
          end if
        end do

        call store (outr, i, qrain)
        if (qrain .gt. qrmax) qrmax = qrain

      end do

      k = 1
C      m = 1
C      time = 0.
C** needs initial value/ 4/03
C      rin = r(1)
      if(diag) write(99,876)

      rin = 0.
      d2  = 0.

      if (pond_rain) then
        d1 = z(1)
        ir = 1
      end if

C     compute initial surface area and discharge

      v1 = voli
      v2 = voli

      if (v2 .lt. v(2)) then
        k = 1
      else
C      
   10   if (v2 .gt. v(k+1)) then
          k = k + 1
          if (k .eq. nd) then
            k = nd - 1
            go to 15
          end if
          go to 10
        else if (v2 .lt. v(k)) then
          k = k - 1
          go to 10
        end if

      end if

   15 continue

      if (q(k) .gt. 0. .and. v1 .gt. 0.) then
        q1 = q(k) * exp (( (alog (v1) - vlog(k)) /
     &       (vlog(k+1) - vlog(k)) ) * (qlog(k+1) - qlog(k)))
      else
        q1 = q(k) + ((v1 - v(k)) /
     &      (v(k+1) - v(k))) * (q(k+1) - q(k))
            end if
      if (v(k) .gt. 0.) then
         sa1 = sa(k) * exp(((alog(v1) - vlog(k)) /
     &         (vlog(k+1) - vlog(k))) * (salog(k+1) - salog(k)))
      else
        sa1 = sa(k) + ((v1 - v(k)) /
     &        (v(k+1) - v(k))) * (sa(k+1)-sa(k))
      end if

      call store (outq(1), 1, q1)

      WRITE( 1, '(5(F12.2,","),F12.2)' ) 0., v1, sa1, 0., 0., q1

      v_max = q1
C                                 start fixed time step loop (i, delt, step)...
C------------------------------------------------------------------------------

      do i = 2, limit

C        step = float (i - 1) * delt
        dt = delt

        time = float (i - 1) * delt

        if (pond_rain) then
C                                    compute rainfall rate as average over delt
          do while(time .gt. t(ir))
            ir = ir + 1
          end do

          d2 = z(ir-1) +  r(ir-1) * (time - t(ir-1))
          rin = (d2 - d1) / delt
          if (rin .lt. 0.) rin = 0.

          r_cum = r_cum + (d2 - d1)

          d1 = d2

        end if
C                                                            upstream inflow(s)
        qup2 = 0.

        do j = 1, nu

          call get (upq(j,1), i, qu(j,1))

          qup2 = qup2 + qu(j,1)

          if (upq(j,2) .gt. 0) then
C                                                     overbank section upstream
            call get (upq(j,2), i, qu(j,2))

            qup2 = qup2 + qu(j,2)

          end if

        end do

        if (sed) then

          do n = 1, nps

            cup2(n) = 0.

            do j = 1, nu

              call get (upc(j,n,1), i, cj)

              if (qup2 .gt. 0.) cup2(n) = cup2(n) + qu(j,1) * cj / qup2

              if (upq(j,2) .gt. 0) then
C                                                     overbank section upstream
                call get (upc(j,n,2), i, cj)

                if (qup2 .gt. 0.) cup2(n) = cup2(n) + qu(j,2) * cj
     &                                      / qup2
              end if

            end do

          end do

        end if
C                                                                lateral inflow
        qlat2 = 0.

        do j = 1, nl

          call get (latq(j), i, ql(j))

          qlat2 = qlat2 + ql(j)

        end do

        if (sed) then

          do n = 1, nps

            clat2(n) = 0.

            do j = 1, nl

              call get (latc(j,n), i, cj)

              if (qlat2 .gt. 0.) clat2(n) = clat2(n) + ql(j) * cj
     &                                      / qlat2
            end do

          end do

        end if
C                                start variable time step loop (m, dt, time)...
C------------------------------------------------------------------------------

C30      if (time .eq. t(m)) then
C                                                   last time = last breakpoint
C                                              next time = next fixed time step
C          rin = r(m)
C          dt = step - time
C          m = m + 1
C        end if

C        if (time .lt. t(m) .and. time + dt .gt. t(m)) then
C                                                                   next time =
C                                                               next breakpoint
C          dt = t(m) - time
C          time = t(m)
C        else
C          time = step
C        end if
C                                                               compute average
C                                                                  inflow rates
C        qu2 = qup2 - (qup2 - qup1) * (step - time) / delt
C        qup = 0.5 * (qu1 + qu2)
C        ql2 = qlat2 - (qlat2 - qlat1) * (step - time) / delt
C        qlat = 0.5 * (ql1 + ql2)

        qu2 = qup2
        qup = 0.5 * (qu1 + qu2)
        ql2 = qlat2
        qlat = 0.5 * (ql1 + ql2)

        if (sed) then

          do n = 1, nps
C            cu2(n) = cup2(n) - (cup2(n) - cup1(n)) * (step - time)
C     &               / delt
C            cl2(n) = clat2(n) - (clat2(n) - clat1(n)) * (step-time)
C     &               / delt
            cu2(n) = cup2(n)
            cl2(n) = clat2(n)
          end do

        end if

C                                   Solution for v2, q2 based on continuity and
C                                                 interpolation of rating table
C------------------------------------------------------------------------------

        IF( i == 63 )THEN
          n = 1
        END IF

        if (qup + qlat + rin + v2 .gt. 0.) then

          source = ' Pond '

          call iter (erpond, v2, VMIN, VMAX, ierr, source)

          if (ierr .gt. 0)
     &    call errxit(msg(ms:13), 'no convergence (erpond)')

          if (v2 .gt. v(nd)) then

C           pond is overtopped, set outflow = inflow

            v2  = v(nd)
            q2  = qu2 + ql2 + samax * (rin - fmin)
            sa2 = samax
            if (iexceed .eq. 0) iexceed = i
          else
            if (q(k) .gt. 0. .and. v2 .gt. 0.) then
              q2 = q(k) * exp (( (alog (v2) - vlog(k)) /
     &             (vlog(k+1) - vlog(k)) ) * (qlog(k+1) - qlog(k)))
            else
              q2 = q(k) + ((v2 - v(k)) /
     &             (v(k+1) - v(k))) * (q(k+1) - q(k))
            end if
            if (v(k) .gt. 0.) then
              sa2 = sa(k) * exp(((alog(v2) - vlog(k)) /
     &              (vlog(k+1) - vlog(k))) * (salog(k+1) - salog(k)))
            else
              sa2 = sa(k) + ((v2 - v(k)) /
     &              (v(k+1) - v(k))) * (sa(k+1)-sa(k))
            end if
          end if

          if (v2 .gt. v_max) v_max = v2

          if (v2 .lt. 1.E-6) v2 = 0.

        else

          q2 = 0.
          sa2 = sa(1)

        end if

C      WRITE( 8, '(2F10.2)' ) time / 60., v2
C##
        if(diag) then
          write(99,877)time,qlat,qup,q2,v2,sa2
 876  format(' time/60.  qlat     qup      qout       vol       surfA ')
          write(99,'("k ",I4)') k
 877  format(f8.2, 3g11.3, 2f10.2) 
        end if

        if (sed) then
C                                                    calculate sediment balance
          do n = 1, nps
            if (sa2 .gt. 0.) then
              qsin1(n) = qsin2(n)
              qsin2(n) = cu2(n) * qu2 + cl2(n) * ql2
              c2(n) = (qsin2(n) + qsin1(n) - c1(n)* (q1 + sa1 * vsl(n) -
     &                2.* v1 / dt)) / (2.* v2 / dt + q2 + sa2 * vsl(n))
            else
              c2(n) = 0.
            end if
          end do

        end if

C                                                                  flow volumes
        vi = vi + dt * (qup + qlat)
        vo = vo + dt * 0.5 * (q1 + q2)
        sarea = 0.5 * (sa1 + sa2)
        darea = samax - sarea
        dvf = fmin * sarea + amin1(fmin * darea, qlat + rin * darea)

        if (v2 .gt. 0.) then
          vf = vf + dt * dvf
        else
          vf = vf + dt * amin1 (dvf, qup + qlat + rin * samax)
        end if
C                                                     wt of sediment components
        if (sed) then

          dso = 0.
          twso = 0.

          do n = 1, nps
            wsi = wsi + (qsin1(n) + qsin2(n)) * .5 * dt * rho(n)
            dso = (q1 * c1(n) + q2 * c2(n)) * .5 * rho(n)
            wse = wse + (c1(n) * sa1 * vsl(n) + c2(n) * sa2 * vsl(n))
     &            * .5 * dt * rho(n)
            twso = twso + dso
            wso(n) = wso(n) + dso * dt
          end do

        end if
C                                                               reset variables
        q1 = q2
        v1 = v2
        sa1 = sa2
        qu1 = qu2
        ql1 = ql2

        if (sed) then

          do n = 1, nps
            cu1(n) = cu2(n)
            cl1(n) = cl2(n)
            c1(n) = c2(n)
          end do

        end if

C        if (time .ne. step) go to 30
C                                             ...end of variable time step loop
C------------------------------------------------------------------------------

        call store (outq(1), i, q2)
C                                                           check for peak flow
        if (q2 .gt. qp) then
            qp = q2
            tp = time
        end if

        qup1 = qup2
        qlat1 = qlat2

        if (sed) then

          do n = 1, nps

            call store (outc(n), i, c2(n))

            cup1(n) = cup2(n)
            clat1(n) = clat2(n)

          end do

          if (twso .gt. qsp) then
            qsp = twso
            tsp = time
          end if

        end if

        WRITE(1,'(5(F12.2,","),F12.2)') time/60.,v2,sa2,vf,r_cum*conv,q2

      end do
C                                                ...end of fixed time step loop
C------------------------------------------------------------------------------

      CLOSE( 1 )

      if (vf .lt. 0.) vf = 0.

      qbal(2) = d2 * rarea
      qbal(7) = vf
      qbal(9) = v2 - voli
      qbal(10) = vo
C                                                         pass values to writer

      call qwrt (id, 1, msg, qp, tp, vi, coarea, outq, 3, outr, qrmax)

      if (sed) then
C                                                      wt of suspended sediment
        wss = 0.
        do n = 1, nps
            wss = wss + v2 * c2(n) * rho(n)
        end do
C                                                         pass values to writer
        call swrt ( wso, wsi, wse, wss, qsp, tsp, outc, wse, 0.)

      end if

      IF( v_max > v(nd) )THEN
        xpond_n = xpond_n + 1
        xpond_id(xpond_n) = id
        xpond_pct(xpond_n) = (v_max - v(nd)) / v(nd) * 100.
      END IF

      return

      end


C------------------------------------------------------------------------------


      subroutine erpond (v2, fv2, dfv2)

C   Computes the residual fv2 and derivative w/resp to v2 of the error function
C   form of the continuity equation.

      common /pond1/ v(50), vlog(50), q(50), qlog(50), sa(50), dt, fmin,
     &               salog(50), qup, k, v1, q1, sa1, qlat, vthres, sath,
     &               nd, i, samax, qu1, qu2, ql1, ql2, rin


C                                      find index k to interval in rating table
      if (v2 .lt. v(2)) then
        k = 1
      else
C      
   10   if (v2 .gt. v(k+1)) then
          k = k + 1
          if (k .eq. nd) then
            k = nd - 1
            go to 15
          end if
          go to 10
        else if (v2 .lt. v(k)) then
          k = k - 1
          go to 10
        end if

      end if

   15 continue
C                                   compute discharge by interpolation; compute
C                                             derivative with respect to volume
C------------------------------------------------------------------------------

      if (q(k) .gt. 0.) then
        q2 = q(k) * exp (((alog (v2) - vlog(k)) /
     &       (vlog(k+1) - vlog(k))) * (qlog(k+1) - qlog(k)))
        dqdv = (q2 / v2) * ((qlog(k+1) - qlog(k)) /
     &         (vlog(k+1) - vlog(k)))
      else
        q2 = q(k) + ((v2 - v(k)) / (v(k+1) - v(k))) * (q(k+1) - q(k))
        dqdv = (q(k+1) - q(k)) / (v(k+1) - v(k))
      end if
C                                         compute surface area by interpolation
C------------------------------------------------------------------------------

      if (v(k) .gt. 0.) then
        sa2 = sa(k) * exp (((alog (v2) - vlog(k)) /
     &        (vlog(k+1) - vlog(k))) * (salog(k+1) - salog(k)))
      else
        sa2 = sa(k) + ((v2 - v(k)) /
     &        (v(k+1) - v(k))) * (sa(k+1) - sa(k))
      end if
C                                                    compute error function fv2
C------------------------------------------------------------------------------

      if (v1 .lt. vthres .and. v2 .gt. vthres) then

C                                           recompute dt as the fraction of the
C                                    original dt remaining after outflow begins

        qin1 = qu1 + sa1 * (rin - fmin) + ql1 +
     &         amax1 (0., (rin - fmin) * (samax - sa1))
        qin2 = qu2 + sa2 * (rin - fmin) + ql2 +
     &         amax1 (0., (rin - fmin) * (samax - sa2))
        aa = 0.5 * (qin2 - qin1) / dt
        bb = qin1
        cc = -(vthres - v1)

        if (aa .le. 0.) then

          if (bb .eq. 0.) then
C                                                  no inflow - v2 estimate from
C                                                  subroutine iter is incorrect
            go to 20

          else
            dtp = -cc / bb
          end if

        else

            dtp = (-bb + sqrt (bb * bb - 4. * aa * cc)) * 0.5 / aa

        end if

        if (dtp .ge. dt) then
C                                 v2 estimate from subroutine iter is incorrect
            go to 20
        else
            dtn = dt - dtp
        end if

C        sarea = 0.5 * (sath + sa2)
C        fv2 = dtn * (qup + rin * sarea + amax1 (0., qlat + (rin - fmin)*
C     &        (samax - sarea)) - 0.5 * q2 - fmin * sarea) -
C     &        (v2 - vthres)

        sarea = 0.5 * (sa1 + sa2)
        fv2 = dt * (qup + rin * sarea + amax1 (0., qlat + (rin - fmin) *
     &        (samax - sarea)) - fmin * sarea) - 0.5 * dtn * q2 -
     &        (v2 - v1)

        go to 30

      else if (v1 .gt. vthres .and. v2 .lt. vthres) then

C                                           recompute dt as the fraction of the
C                                    original dt remaining after outflow ceases

        qin1 = qu1 + sa1 * (rin - fmin) + ql1 +
     &         amax1 (0., (rin - fmin) * (samax - sa1))
        qin2 = qu2 + sa2 * (rin - fmin) + ql2 +
     &         amax1 (0., (rin - fmin) * (samax - sa2))
        aa = 0.5 * (qin2 - qin1) / dt
        bb = 0.5 * q1 - qin1
        cc = -(v1 - vthres)

        if (aa .le. 0.) then
          dtp = -cc / bb
        else
          dtp = (-bb + sqrt (bb * bb - 4. * aa * cc)) * 0.5 / aa
        end if

        if (dtp .ge. dt) then
C                                 v2 estimate from subroutine iter is incorrect
          go to 20
        else
          dtn = dt - dtp
        end if

C        sarea = 0.5 * (sath + sa2)
C        fv2 = dtn * (qup + rin * sarea + amax1 (0., qlat + (rin - fmin)*
C     &        (samax - sarea)) - 0.5 * q1 - fmin * sarea) -
C     &        (v2 - vthres)

        sarea = 0.5 * (sa1 + sa2)
        fv2 = dt * (qup + rin * sarea + amax1 (0., qlat + (rin - fmin) *
     &        (samax - sarea)) - fmin * sarea) - 0.5 * dtp * q1 -
     &        (v2 - v1)

        go to 30

      else
C                                                               all other cases
          go to 20

      end if

20    sarea = 0.5 * (sa1 + sa2)

      dtn = dt

      fv2 = dtn * (qup + rin * sarea + amax1 (0., qlat + (rin - fmin)*
     &      (samax - sarea)) - 0.5 * (q1 + q2) - fmin * sarea) -
     &      (v2 - v1)
C                                                   compute derivative of error
C                                                         function w/resp to v2
30    dfv2 = -0.5 * dtn * dqdv - 1.

      return

      end
