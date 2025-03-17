MODULE kinsed_def

  USE runpars
  USE elpars

  IMPLICIT NONE
  PUBLIC

  INTEGER :: sed_model ! 1 = Smith, 2 = DRHEM, 3 = DWEPP

  CONTAINS

!---------------------------------------------------------------------------------------------------
!   PUBLIC PROCEDURES
!---------------------------------------------------------------------------------------------------

    SUBROUTINE sed00()

!   Get global sediment parameters

      INTEGER :: ierr,i,j,m,na
      REAL    :: sumvs

      nps = 0
      DO i = 1, 5
        CALL getr4('DI', i, dpd(i), ierr) ! particle diameter (mm or in)
        IF(ierr .gt. 0) EXIT
        CALL getr4('DE', i, rho(i), ierr) ! particle density (g/cc)
        IF(ierr .gt. 0) STOP ' error - missing particle density'
        nps = i
      END DO
      IF(nps .lt. 1) STOP ' error - particle classes not defined'

      sumvs = 0.
      DO i = 1, nps
        dpd(i) = dpd(i) / conv
        vsl(i) = vsetl(rho(i), dpd(i), xnu, grav)
        sumvs  = sumvs + 1. / vsl(i)
      END DO
      
!     Compute relative erodibility weighting

      DO i = 1, nps
        pwt(i)  = 1. / vsl(i) / sumvs
        nord(i) = i
      END DO

      IF(nps .gt. 1) THEN ! order classes by increasing erodibility
        j = 1
        DO WHILE(j + 1 .le. nps)
          m = j
          DO i = j + 1, nps
            If(pwt(i) .gt. 0. .and. pwt(i) .lt. pwt(m)) m = i
          END DO
          IF (m .ne. j) THEN
            CALL switch(rho(j), rho(m))
            CALL switch(vsl(j), vsl(m))
            CALL switch(dpd(j), dpd(m))
            CALL switch(pwt(j), pwt(m))
            na = nord(j)
            nord(j) = nord(m)
            nord(m) = na
          END IF
          j = j + 1
        END DO
      END IF

    END SUBROUTINE sed00

!---------------------------------------------------------------------------------------------------

    SUBROUTINE sed0(msg, inflow, nup, nlat, nx, omeg, delx, slp, qf)

      CHARACTER(LEN=*) :: msg
      INTEGER          :: nup,nlat,nx,inflow(12)
      REAL             :: omeg,delx,qf,slp(20,nchan)

      INTEGER :: ierr
      REAL    :: ki

      CALL getr4('KE', 0, ki, ierr) ! Ks tag used with RHEM
	
      IF(ierr /= 0) THEN ! Smith or DWEPP
        CALL getr4('KR', 0, ki, ierr) ! tag used with DWEPP
        IF(ierr /= 0) THEN ! Smith
          CALL sed0_smith(msg, inflow, nup, nlat, nx, omeg, delx, slp, qf)
          sed_model = 1
        ELSE ! DWEPP
          CALL sed0_wepp(msg, inflow, nup, nx, omeg, delx, slp, qf)
          sed_model = 3
        END IF
      ELSE ! DRHEM
        CALL sed0_rhem(msg, inflow, nup, nx, omeg, delx, slp, qf)
        sed_model = 2
      END IF

    END SUBROUTINE sed0

!---------------------------------------------------------------------------------------------------

    SUBROUTINE kinsed(indx, dt, dte, ql, a2, q2, topw, rf, qup, qt, time)

      INTEGER :: indx
      REAL    :: dt,dte,ql(20,2),a2(20,nchan),q2(20,nchan),topw(20,nchan),rf,qup(10,2),qt(20),time

      IF(sed_model == 1) THEN ! Smith
        CALL kinsed_smith(indx, dt, dte, ql, a2, q2, topw, rf, qup, qt, time)
      ELSE IF(sed_model == 2) THEN ! RHEM
        CALL kinsed_rhem(indx, dt, dte, ql, a2, q2, rf, qup, time)
      ELSE
        CALL kinsed_wepp(indx, dt, dte, ql, a2, q2, rf, qup, time)
      END IF

    END SUBROUTINE kinsed

!---------------------------------------------------------------------------------------------------

    SUBROUTINE sedfin()

      IF(sed_model == 1) THEN ! Smith
        CALL sedfin_smith()
      ELSE IF(sed_model == 2) THEN ! RHEM
        CALL sedfin_rhem()
      ELSE
        CALL sedfin_wepp ()
      END IF

    END SUBROUTINE sedfin

!---------------------------------------------------------------------------------------------------

    SUBROUTINE sedbf(qup, a2, q2, topw)

      REAL :: qup(10,2),a2(20,nchan),q2(20,nchan),topw(20,nchan)

      IF(sed_model == 1) THEN ! Smith
        CALL sedbf_smith(qup, a2, q2, topw)
      END IF

    END SUBROUTINE sedbf

!---------------------------------------------------------------------------------------------------
!   PRIVATE PROCEDURES
!---------------------------------------------------------------------------------------------------

   REAL FUNCTION vsetl(ss, d, vnu, grav)

!  Finds settling velocity for a round particle of diameter d and specific
!  gravity ss, in water of viscosity vnu.

      REAL :: ss,d,vnu,grav

      REAL    :: ca,ck,cb,cc,vs,f,fp,dvs
      INTEGER :: i

      ca   = 24. * vnu / d
      ck   = 4. / 3. * grav * (ss - 1.) * d
      cb   = 3. * SQRT(vnu / d)
      cc   = .34
      vs   = 2.8E05 * d * d

      DO i = 1, 40

        f = ca * vs + cb * vs**1.5 + cc * vs * vs - ck
        fp = ca + 1.5 * cb * sqrt (vs) + 2. * cc * vs
        dvs = f / fp

        IF(ABS (dvs) .gt. 0.0001 * vs) THEN
          vs = vs - dvs
          IF(vs .le. 0.) vs = 0.5 * (vs + dvs)
        ELSE
          vsetl = vs
          RETURN
        END IF

      END DO

      STOP ' error computing Vsetl for given Temp. and size '

   END FUNCTION vsetl

END MODULE kinsed_def
