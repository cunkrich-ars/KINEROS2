C   Code type: Fortran subroutine

C   Compiler: Fortran77

C   Initial code programmed by R. Smith for sediment routine was modified 
C   by N. Bulygina to fit WEPP-like equations for plane and later 
C   by Mariano Hernandez to fit RHEM Unit Stream Power equations

C   Date: 12/2004

C   Date: 09/2009

C   Date: 08/2011

C   Description:

C     A four-point explicit finite-difference solution of the equation of
C     supply and conservation of transported material. The solution is
C     explicit since advantage is taken of the kinematic wave solution for the
C     same time step. This routine is called from PLANE modules.
C     This version treats sediment supply drawn from as
C     many as 5 different particle size classes, thus to deal with a watershed
C     where soils are not uniform for all elements.
C------------------------------------------------------------------------------

C   Entry points:

C     sed0_plane          reads parameters & initializes vars for current element,

C     sedfin_plane        finish sediment balance computations for current element,
C------------------------------------------------------------------------------

C   Arguments:


C     id             int          element identifier,
C     msg            char*(*)     element id string for error messages,
C     indx           int          current time index,

C     inflow(10)     int          inflow(1:10)  = upstream element id's

C     units          int          1 = metric, 2 = english,
C     qup(10)        real         upstream inflow rate(s) at current time,
C     dt             real         current time increment (sec),
C     dte            real         time elapsed since last user-defined time step,
C     ql(30)         real         lateral inflow of water ...
C                                 ... for a plane, the average rate over the time
C                                     interval, m/s or ft/s;
C     qil(30)        real         bed or soil inflow of water (m/s or ft/s),
C     a2(30)         real         area at end of interval divided by 
C                                 rill spacing (ams) or depth of overland flow,
C     q2(30)         real         discharge per unit width or ams at end of time interval, (m2/s)
C     wid(30)        real         width of plane
C     rf             real         current rainfall rate (m/s or ft/s),


C     time           real         current simulation time,

C     nx             int          number of spatial nodes,
C     omeg           real         time weighting factor from finite-diff. eqn.
C                                        of calling routine,
C     delx           real         spatial increment (m or ft),
C     slp(30)        real         bottom slope at each spatial node.
C     units          int          1 = metric, 2 = english,
C     del            real         fixed time step (sec),
C     lim            int          number of (fixed) time steps,
C     nps            int          number of particle classes
C     dpd(5)         real         particle diameters (m or ft),
C     rho(5)         real         particle densities,
C------------------------------------------------------------------------------

C   Input Block (input file labels in parenthesis):

C     temp               real     water temperature (C or F),
C     pros(5)   (FR)     real     fraction of total material represented by
C                                     each particle class 

C     Kib  (KSS)          real     splash and sheet erodibility coefficient (dimensionless)
C     Kcu  (KOM)         real     undisturbed concentrated erodibility coefficient (s^2/m^2)
C     Kcm  (KCM)         real     maximum concentrated erodibility coefficient (s^2/m^2)
C     Adf  (ADF)         real     beta decay factor in the exponential decay function to estimate
C                                 the erodibility coeff. Eq (13) Al-Hamdan et al. (2012).  

C     rsp  (RSP)          real    rill spacing (m or ft),

C     
C------------------------------------------------------------------------------

C   Subroutines/functions:

C     Intrl    (kinsed_plane.for)
C     Dc                "  
C     vsetl        (kinsed.for) 
C     switch            "
C     getr4        (reader.for)
C     swrt         (writer.for)
C     errxit       (K2RunW.for)
C------------------------------------------------------------------------------


      subroutine kinsed_rhem ( indx, dt, dte, ql, a2, q2, rf, 
     &                 	qup, time )

      use runpars
      use elpars
      use multip
C                                                                     arguments
C------------------------------------------------------------------------------

      integer  indx, nx, nc
      real dt, dte, rf, time, delx
      real, dimension(30) ::  a2, q2, slp
      real, dimension(30) ::  ql, sumq2
      real, dimension(10) :: qup
      integer, dimension(10) :: inflow

      character msg*(*)

      real, dimension(30) :: qso
C                                                               local variables
C------------------------------------------------------------------------------
C  qr2- discharge in a rill 
C  ar2- x-sectional area of flow in a rill 
C  h2 - depth of flow at given time and given node    
      real, dimension(30) :: qr2, ar2, h2, eros, velocity, width
      integer :: ierr, i, im, j, jp, n, na, nup, nlat, redist
      integer, dimension(5):: excd
      integer,save :: nk, ilast, nu
C
	real Dc_v, IR_v, coef1,coef2,coef3,sp_nearing,
     &	 eterm, dtrc, distance,
     &     s1, s2, s3, cmax, uu, dz, rs,
     &     wStreamPower, qss, qs, qs01, qs02
      real, save::  omega, comega, dx, qsp, tsp, xl, qfac,
     &              Kcu, Kcm, rfp, rsp, bare, adf,
     &              Kibadj, Ki, alpha, p 
      real,save,dimension(5) :: wso, dwso, dwsol, wsin, dwsup, dwsupl
      real,save,dimension(5) :: csub, pros, csub1, csub2, dcsub
      real,save,dimension(30) :: depnet, slope, a1, q1, qil, Fr, prob
      real,save,dimension(30) :: depnet_neg, depnet_pos
      real,save,dimension(5,30) :: cs, cl 
      real,save,dimension(30):: cst, clt, erosp, widthp
	real cpt, cpmt
	real,save,dimension(30) ::qr1, ar1, h1
      real, dimension(5):: cp, cpm, qsi
C
      integer, save, dimension(5) :: outc
      integer, save, dimension(10,5) :: upc
      
C
      character(LEN=2),dimension(5,1) :: attr = reshape((/'S1', 'S2', 
     &'S3', 'S4', 'S5'/),(/5,1/))
C------------------------------------------------------------------------------
C     Flow path width on the rectangular channel changes with time
C     and slope length 
C------------------------------------------------------------------------------
 
C     Equation (34) in Al-Hamdan et al. 2012.
C     Characteristics of concentrated flow hydraulics for rangeland ecosystems:
C     implications for hydrologic modeling
      width(nk)=(2.46*(q2(nk)*rsp)**0.39)/(sin(atan(slope(nk)))**0.40)
	if (width(nk). gt. rsp) width(nk) = rsp
C
C     Discharge, x-sectional area and depth of flow per rill space
        
      do j=1,nk
	  sumq2(j) = sumq2(j) + q2(j)
	  qr2(j) = q2(j) * rsp
        width(j)=( 2.46*(q2(j)*rsp)**0.39)/(sin(atan(slope(j)))**0.40 )
C          
         if ( width(j) .gt. rsp ) width(j) = rsp
c         compute rill flow depth h2(j) for rectangular rills. This is an iterative
c         process to solve the uniform flow equation:
c
c         h2(j)=((qr2(j)/alpha/sqrt(slope(j))**(1/(p+1))/width)*
c	    (width+2*h2(j))**(p/(p+1))
C         
          if ( qr2(j)  .le. 0. ) then
           h2(j) = 0.0
          else
           uu = (qr2(j)/alpha/sqrt(slope(j))) ** (1./(p+1.))/width(j)
           h2(j) = 0.2 * (q2(nk) * rsp) ** 0.36
   2       dz = h2(j)
           h2(j) = uu * (width(j) + 2.*dz) ** (p/(p+1.))

           if (abs(dz/h2(j)-1.) .gt. 0.000001) go to 2
          end if

	     ar2(j) = h2(j) * width(j)

        end do
C
C       Probability that overland flow becomes concentrated flow.
C       Equation (14) in Al-Hamdan et al. 2013
C       Risk Assessment of Erosion from Concentrated Flow on Rangelands ..
C
        do j=1,nk
            
           prob(j) = probability (slope(j), bare, h2(j), q2(j))
        
        end do
C       
      if (indx .ne. ilast) then
C                                             store outgoing concentration from
C                                                   end of last fixed time step
        do n = 1, nps
         call store ( outc(n), ilast, cs(n,nk)/(rho(n)*1000.) )
	  end do
         
	  do i = 1, nk
          qso(i) = 0.
          do n = 1, nps
            qso(i) = qso(i) + cs(n,i) * q2(i) * wt / 1000.
          end do
        end do
       
C        write (88, '(F8.2,15(F10.6))') (indx - 1) * delt / 60.,
C     &                                 ( qso(i), i = 1, nk )
      end if

      dt2 = dt / 2.
      
	do n = 1, nps
          csub(n) = 0.
      end do

C                                             lateral inflow is average over dt
        do i = 2, nk
          qil(i) = ql(i)
        end do

      if (nu .gt. 0) then
C                                       compute concentration at upper boundary
C------------------------------------------------------------------------------

        
        if (indx .ne. ilast) then
C                                                        next fixed time step
          do n = 1, nps
            csub1(n) = csub2(n)
          end do
        
	    if (qup(1) .gt. 0.) then
             do n = 1, nps
               call get (upc(1,n), indx, csub2(n))
		   end do
          else

            do n = 1, nps
              csub2(n) = 0.
            end do
          end if

          do n = 1, nps
            dcsub(n) = (csub2(n) - csub1(n)) / delt
          end do

        end if   
          
C                                                                   interpolate
          do n = 1, nps
            csub(n) = csub1(n) + dcsub(n) * dte
		end do

      end if
      
      ilast = indx
C                            BEGIN SEDIMENT ROUTING -- upper boundary condition
C------------------------------------------------------------------------------

      if (a2(1) .lt. 1.E-8) then
C                                                            no upstream inflow
        cst(1) = 0.  
	  do n = 1, nps
           cs(n,1) = 0.
	  end do

      else                     ! positive upper boundary area:
C                                        upstream flow determines concentration
        cst(1) = 0.
	  do n = 1, nps
            cs(n,1) = csub(n) * rho(n) * 1000.
 	      if (cs(n,1). lt. 1.e-15) cs(n,1) = 1.e-15
            cst(1) = cst(1) + cs(n,1)	  
	  end do
         
      end if
C                                                 loop through spatial nodes...
C------------------------------------------------------------------------------

      do i = 2, nk

       im = i - 1                
       distance = float(i-1) * dx

      if (ar2(i) .gt. 1.E-5) then
		do n = 1,nps
		 cp(n) = cl(n,i)
           cpm(n) = cl(n,im)
		end do 
		 cpt = clt (i)
           cpmt = clt (im)
	
C                                           condition at beginning of time step
C------------------------------------------------------------------------------

C                                                      rain splash contribution
	 if (a1(i) .lt. 1.E-8 .and. qil(i) .gt. 0.0) then
		cpt = 0.
		do n=1,nps
    	       cp(n) = (Kibadj * (rfp**1.052)* (qil(i)**0.592))
     &	   /((qil(i)**0.592)+(0.5*vsl(n))) 
		   cpt = cpt + (Kibadj * (rfp**1.052) * (qil(i)**0.592))
     &	   /((qil(i)**0.592)+(0.5*vsl(n)))
  
		  end do
  
	 end if
  
       if (a1(im) .lt. 1.E-8 .and. qil(im) .gt. 0.0) then
		cpmt = 0.
		do n=1,nps
              cpm(n) = (Kibadj * (rfp**1.052) * (qil(im)**0.592))
     &		/((qil(im)**0.592)+(0.5*vsl(n)))
             cpmt = cpmt + (Kibadj * (rfp**1.052) * (qil(im)**0.592))
     &	    /((qil(im)**0.592)+(0.5*vsl(n)))
		  end do
	 end if

C-------------------------------------------------------------------
C     Computation of flow velocity in concentrated flow areas (m/s)
C     V=Q/A (m/s)
C     Q=qr2 (m^3/s)
C     A=ar2 (m^2)
      velocity(i) = qr2(i)/ar2(i)
C---------------------------------------------------------------------------------------------------
C     Unit Sediment Load qs(kg s-1 m-1)
C     Stream Power (Eq 4) paper Nearing et al. 1997
C     w(g/s^3)=density of water(g/cm3) * gravity (cm/s2) * slope * unit discharge (cm2/s), width (cm)
C     stream power was computed using cm and seconds, 

         wStreamPower = 1.*980.7*sin(atan(slope(i)))*(qr2(i)/width(i))*
     &                            10000. ! convert m2 to cm2 in qr2   
        
         wStreamPowerKgpers3=wStreamPower*(1./1000.) ! Convert stream Power from g s-3 to kg s-3
	   vstreampowerKgpers3=streampower_velocity( slope(i), width(i),
     &   h2(i), velocity(i))
	      qs01 = 38.61 * exp( 0.845 + 0.412 * log10(wStreamPower) )   ! wStreamPower is in gm/s^3
	      qs02 =    (1.+ exp( 0.845 + 0.412 * log10(wStreamPower) ) ) ! wStreamPower is in gm/s^3
	      qss = -34.47 + qs01/qs02 
            qs = ( 100./1000. ) * 10.**(qss) ! convert from (g s-1 cm-1) to (kg s-1 m-1)

C ----------------------------------------------------------------------------- 
C                                         compute erosion or deposition term
C--------------------FOR Stream Power -----------------------------------------
	 Dc_v = Dc ( Kcu, kcm, prob(i), adf, sumq2(i),
     &             wStreamPowerKgpers3 )
C---------------------------------------------------------------------------------
      if (qil(i) .gt. 0.) then
          IR_v = Kibadj * (rf**1.052) * (qil(i)**0.592) * rsp
      else
          IR_v = 0.
      end if
C*******************************************************************************
C      rill/interrill detachment, not deposition 
C*******************************************************************************
       if (qr2(i) .gt. 0.) then
           
           eterm = Dc_v * width(i)
           dtrc = eterm + IR_v ! dtrc=kg/(s*m)
           deno = ar2(i)/dt + omega/dx*qr2(i) + Dc_v/qs * qr2(i)
C----------------------------------------------------------------------------------          
	     coef1 = ar1(i)/dt - comega/dx*qr1(i) ! coef1=m^2/s
	     coef2 = omega/dx * qr2(im) ! coef2=m^2/s
	     coef3 = comega/dx * qr1(im) !coef3=m^2/s
           cst(i) = (cpt * coef1 + cst(im) * coef2 + 
     &               cpmt * coef3 + dtrc)/deno ! cst=(kg/m^3)
           
       else
           
           eterm = 0.
           dtrc = eterm + IR_v ! dtrc=kg/(s*m)
           deno = ar2(i)/dt + omega/dx*qr2(i)
           coef1 = ar1(i)/dt - comega/dx*qr1(i) ! coef1=m^2/s
	     coef2 = omega/dx * qr2(im) ! coef2=m^2/s
	     coef3 = comega/dx * qr1(im) !coef3=m^2/s
           cst(i) = (cpt * coef1 + cst(im) * coef2 + 
     &               cpmt * coef3 + dtrc)/deno ! cst=(kg/m^3)
           
       endif

       if (cst(i) .lt. 1.e-15) cst(i) = 1.e-15

C                                                         compute concentration
C------------------------------------------------------------------------------
  
C      Detachment is considered as nonselective process, therefore
C      all detached sediment on every segment is divided between classes
C      according to there fraction in the soil           
	  do n = 1, nps
	 ! incoming sediment is averaged over time
	     cs(n,i) = ( cs(n,im) + cl(n,im) ) / 2. + 
     &		       ( cst(i) - (cst(im) + clt(im)) / 2. ) * pros(n)
	    if (cs(n,i) .lt. 1.e-15) cs(n,i) = 1.e-15  
        end do	 		         
C     eros is in kg/(s*m^2)
       eros(i) = ((cst(i)*qr2(i)+ clt(i)*qr1(i))-
     &           (cst(im)*qr2(im) + clt(im)*qr1(im)))/  
     &           (dx * (width(i) + widthp(i)))       
     
C------------------------------------------------------------------------------
C Use this statement when stream power is used. gs in kg/(s m) ; qr2 (m2/s); rsp (m)  
       if ( ( qs * width(i) ) .lt. ( cst( i ) * qr2( i ) ) ) then
C------------------------------------------------------------------------------
C                                                         there was deposition 
C-----------------------------------------------------------------------------
C                                                   loop thru particle sizes...
C------------------------------------------------------------------------------
        do n = 1, nps        
C                                                       compute deposition term
C------------------------------------------------------------------------------ 		     
C                                                         compute concentration
C------------------------------------------------------------------------------
          qsi(n) = ( qs/nps ) * ( 1.- pros(n) )
		dtrc = 0.5 * vsl(n) * qsi(n) / qr2(i) * width(i)
     &           + Ir_v * pros(n)

C                                                         compute concentration
C------------------------------------------------------------------------------
	    deno = ar2(i)/dt + omega/dx*qr2(i) +
     &		   0.5 * vsl(n) * width(i) 
	    coef1 = ar1(i)/dt - comega/dx*qr1(i)
	    coef2 = omega/dx * qr2(im)
	    coef3 = comega/dx * qr1(im)

          cs(n,i) = (cp(n) * coef1 + cs(n,im) * coef2 +
     &               cpm(n) * coef3 + dtrc)/deno
           if (cs(n,i) .lt. 1.e-15) cs(n,i) = 1.e-15

	  end do
         
C                                             ...end loop thru particle classes
C------------------------------------------------------------------------------
C check if sediment for each particle class after deposition region does not 
C exceed maximum allowable one (sediment in + sediment from interrill erosion)
        do n = 1,nps
	    excd(n) = 1 !1= particle class load doesn't exceed or equals allowable max
	  end do
  50    s1=0.
	  s2=0.
	  s3=0.
	  redist = 0

	  do n =1, nps
	   s1 = s1 + cs(n,i)
	   cmax = (cs(n,im) * qr2(im) + cl(n,im) * qr1(im)) / 2. + 
     &   	       IR_v * pros(n) * dx 
!        
         if (excd(n). eq. 1) then
	     if ((cs(n,i) * qr2(i)) .gt.cmax) then
	  	  cs(n,i) = cmax / qr2(i)
	      excd(n) = 0
	      redist = 1
            else
! Correction for case when cs(n,i)*qr2(i) is equal to cmax
!	      if ((cs(n,i)*qr2(i)). lt. cmax) s3 = s3 + cs(n,i)
	      s3 = s3 + cs(n,i)
	     end if
         end if 
	   s2 = s2 + cs(n,i)
	  end do   
        
	  if (redist .eq. 1) then
         do n = 1,nps
	     if (excd(n) .eq. 1) then
		 cs(n,i) = cs(n,i) + (s1-s2) * cs(n,i) / s3
	     if (cs(n,i) .lt. 1.e-15) cs(n,i) = 1.e-15  
           end if
	   end do
	   go to 50
	  end if 

	  cst(i) = 0.
	  do n = 1,nps
	     cst(i) = cst(i) + cs(n,i)
	  end do
C     eroded sediment per unit rill area eros=kg/(s*m^2)       
	   eros(i) = ((cst(i)*qr2(i) + clt(i)*qr1(i))-
     &        (cst(im)*qr2(im) + clt(im)*qr1(im)))/  
     &        (dx * (width(i)+widthp(i)))         	    
	 end if   
       
C                                                    ... end of deposition part
C------------------------------------------------------------------------------                		 			 		
      else
C                                                       no flow at current node
	 do n=1, nps
	  cs(n,i) = 0.          
       end do
	  eros(i) = 0.
      end if
C                                         erosion / deposition depth accounting
C------------------------------------------------------------------------------
C     depnet will be multiplied by 1000 before printing out;    
C     for WEPP equations sediment source in the right side of the equation
C     is per unit area of the RILL!, not just unit area. depnet in ton/m2 
		 depnet(i) = depnet(i) - 0.5 * ( eros(i) * width(i) + 
     &	  	         erosp(i) * widthp(i)  ) * dt/rsp/1000.           
C------------------------------------------------------------------------------
      end do !                            ...End Spatial Loop Aug2011      
      
!K2_PROFILER
      if (profiler) then
        do i = 2, nk
          do n = 1, nps
            conc2(n,i) = cs(n,i)/(rho(n)*1000.)
            if (ar2(i) .gt. 1.E-5) then
              eros2(n,i) = ((cs(n,i)*qr2(i) + cl(n,i)*qr1(i))-
     &                     (cs(n,i-1)*qr2(i-1) + cl(n,i-1)*qr1(i-1)))/  
     &                     (dx * (width(i) + widthp(i)))
            else
              eros2(n,i) = 0.
            end if
            dnet2(n,i) = dnet2(n,i) - 0.5 * ( eros2(n,i) * width(i) + 
     &    	             eros1(n,i) * widthp(i)  ) * dt / rsp / 1000.
          end do
        end do
      end if
!K2_PROFILER

C-------------------------------------------------------------------------------
C                                             reset variables for next time step
C------------------------------------------------------------------------------
	  rfp = rf                                
	  do i = 1, nk
	     widthp(i) = width(i)
           erosp(i) = eros(i)
		 a1(i) = a2(i) 
           q1(i) = q2(i)
           qr1(i) = qr2(i)
	     ar1(i) = ar2(i)
	     h1(i) = h2(i)
	     clt(i) = cst(i)
          do n = 1, nps
             cl(n,i) = cs(n,i)
          end do
        end do
      
C                    compute contribution to total outflow at current time step
C------------------------------------------------------------------------------

      twso = 0.

      do n = 1, nps
C      for plane concentration 
C      is in kg/m^3, while for other parts it's dimensionless	  
        dwso(n) = qfac * q2(nk) * cs(n,nk)/(rho(n)*1000.) ! dwso (m3/s)
	  wsold = wso(n)
        wso(n) = wso(n) + dt2 * (dwso(n) + dwsol(n)) !wso (m3)
	
 	  dwsol(n) = dwso(n)
        twso = twso + dwso(n) * rho(n) !twso (tons/s)
	  	  
	end do
C                                             check for peak sediment discharge
      if (twso .gt. qsp) then
        qsp = twso !qsp (ton/s)
        tsp = time
      end if
C                                                 compute contribution to total
C                                                   inflow at current time step
      do n = 1, nps
        dwsup(n) =  qfac * q2(1) * csub(n)

        wsin(n) = wsin(n) + dt2 * (dwsup(n) + dwsupl(n))
        dwsupl(n) = dwsup(n)
      end do
    
      return

C------------------------------------------------------------------------------

      entry sedfin_rhem
C                                          finish sediment balance computations
C------------------------------------------------------------------------------
C                                                  store outgoing concentration
      do n = 1, nps
      call store ( outc(n), ilast, cs(n,nk)/(rho(n)*1000.) )     
	end do
      
      wsi = 0.
C                              compute weight of incoming and outgoing sediment
      do n = 1, nps
        wso(n) = wso(n) * rho(n)
        wsi = wsi + wsin(n) * rho(n)
      end do
      
C -----------------------------------------------------------------------------
      wse_neg = 0.
      wse_pos = 0.
      
      do i = 2, nk
        if (depnet(i) .lt. 0.)then
            wse_neg = wse_neg + depnet(i)
        else 
            wse_pos = wse_pos + depnet(i)
        end if
      end do
      wse_neg = dx * qfac * wse_neg
      wse_pos = dx * qfac * wse_pos
C -----------------------------------------------------------------------------
C
C                                 compute weight of eroded / deposited material
      wse = 0.
      do i = 2, nk
        wse = wse + depnet(i)
      end do

  	wse = dx * qfac * wse
C                                          compute weight of suspended material
      wss = 0.
      dwss1 = h1(1) * cst(1)
      do i = 2, nk
        dwss2 = h1(i) * cst(i)
C       For rills on a plane considered small volume of water has a rectangular
C       x-section, plus there are several rills on a plane = plane width/ams           
	  wss = wss + 0.5 * (dwss2 + dwss1) * widthp(i) * dx/1000. 
        dwss1 = dwss2
      end do
      wss = qfac/rsp * wss
C                                                         pass values to writer
C -----------------------------------------------------------------------------
      call swrt ( wso, wsi, wse, wss, qsp, tsp, outc, wse_neg, wse_pos )

      return

C------------------------------------------------------------------------------


      entry sed0_rhem ( msg, inflow, nup, nx, omeg, delx, slp, qf )

      idu = id
C                                        save some arguments as local variables
C------------------------------------------------------------------------------
        nk = nx
        omega = omeg
        comega = 1. - omega
        dx = delx
        xl = dx * float (nk - 1)
        qfac = qf
        ilast = 1    != 0
        nu = nup

       
C                                               obtain upstream inflow block(s)
C------------------------------------------------------------------------------

        do i = 1, nu

          do n = 1, nps

            inv = abs (inflow(i))
            call old (inv, attr(n,1), upc(i,n), nc, ierr)
            if (ierr .gt. 0) call errxit
     &      (msg, 'upstream sediment concentration not found')

          end do

        end do
     
C------------------------------------------------------------------------------

      do n = 1, nps
C                                                     obtain blocks for outflow
C                                                         concentration storage
        call new (idu, attr(n,1), outc(n))
C                                                            set c = 0 at t = 0
        call store (outc(n), 1, 0.)

      end do

      do i = 1, nk
        slope(i) = slp(i)
      end do

      do n = 1, nps
        pros(n) = 0.
      end do
C                                                       get sediment parameters
C------------------------------------------------------------------------------ 
C                                                               get rill spacing
	     call getr4 ('RSP', 0, rsp, ierr)
           if (ierr .gt. 0) call errxit
     &                (msg, 'rill spacing (RSP) not found')
	   
C                                             Splash and sheet erodibility coeff
	     call getr4 ('KSS', 0, Ki, ierr)
	     if (ierr .gt. 0) call errxit
     &	                       (msg, 'KSS (KSS) not found')
          if (Ki .lt. 0.0) call errxit
     &	                       (msg,'Illigal value for KSS')
	      Kibadj = Ki 

C                                 undisturbed concentrated flow erodibility coeff	     
		 call getr4 ('KOM', 0, Kcu, ierr)
		 if (ierr .gt. 0) call errxit 
     &                          (msg, 'Kom (KOM) not found')
           if (Kcu .lt. 0.0) call errxit 
     &	                      (msg,'Illegal value for Kom')      		 

C                                     maximun concentrated flow erodibility coeff
		 call getr4 ('KCM', 0, Kcm, ierr)
		 if (ierr .gt. 0) call errxit 
     &                          (msg, 'Kcm (KCM) not found')
           if (Kcm .lt. 0.0) call errxit 
     &	                      (msg,'Illegal value for Kcm')
           
C                                                              Alfa decay factor
           call getr4 ( 'ADF', 0, Adf, ierr)                   
	     if (ierr .gt. 0) call errxit
     &	                       (msg, 'Illegal value for Adf')

C                                       Bare fraction of bare soil to total area
           call getr4 ( 'BARE', 0, Bare, ierr)                   
	     if (ierr .gt. 0) call errxit
     &	                       (msg, 'Illegal value for bare')
           
C                                                   Manning roughness for rills          
		call getr4 ('RMA', 0, rs, ierr)

          if (ierr .eq. 0) then
	     p = 2./3.
	     alpha = 1. / rs
		else
C                                                     Chezy roughness for rills 
           call getr4 ('RCH', 0, rs, ierr)
	     if (ierr .eq. 0) then
	      p = 0.5
	      alpha = rs
		 else 
      	  call errxit (msg, 'Roughness coefficient is not found')
           end if  
	    end if
C                                      check that class proportions sum to one
C                                                       particle size fractions
	
          do n = 1, nps
            call getr4 ('FR', n, prost, ierr)
            if (ierr .gt. 0) go to 1
	    do nn = 1,nps
	     if (nord(nn). eq. n)  pros(nn) = prost
		end do
c		  pros(nord(n)) = prost                           changed on 06/06 (NB)
		end do
  1       sum = 0.
          do n = 1, nps
            sum = sum + pros(n)
          end do

          if (abs (sum - 1.0) .gt. 0.05)
     &    call errxit (msg, 'particle size fractions do not add to 1')

C                                                          initialize variables
C------------------------------------------------------------------------------
C Before this given event there was no water, sediment
	do n = 1, nps  
        csub1(n) = 0.
        csub2(n) = 0.
        dcsub(n) = 0.
        dwso(n) = 0.
        dwsol(n) = 0.
        wso(n) = 0.
        wsin(n) = 0.
        dwsup(n) = 0.
        dwsupl(n) = 0.
      end do

      wStreamPower = 0.
      qs = 0. 
      qsp = 0.
      tsp = 0.
       
      do i = 1, nk
        depnet(i) = 0.
        depnet_neg(i)=0.
        depnet_pos(i)=0.
        erosp(i) = 0.
	  do n = 1, nps
          cs(n,i) = 0.
          cl(n,i) = 0. 
        end do
        cst(i) = 0.
        clt(i) = 0.
	  a1(i) = 0.
	  ar1(i) = 0.
	  h1(i) = 0.
	  q1(i) = 0.
	  qr1(i) = 0. 
        sumq2(i) = 0.
        widthp(i) = 0.
        end do
C        widthp = 0.
C       write (88, '(F8.2,15(F10.6))') 0., (q1(i), i = 1, nk) ! just to write 0's
      return
      end

C------------------------------------------------------------------------------
      real function Dc ( Ku, Kmax, pro, df, addq2,
     &                   wStreamPowerKgpers3 )

C     Finds detachment capacity
C     tw/h- width/depth of flow at given time at given node, 
      real Ku, kmax, pro, df, addq2, wStreamPowerKgpers3
C------------------------------------------------------------------------------
C     Use stream power to estimate Dc -> Dc=kg/(s*m^2)

         if ( df .le. 0.0 ) then
             
         Dc = ku * wStreamPowerKgpers3    
         
         else
             
	   Dc = ( ( 1. - pro ) * Ku +
     &          ( pro * Kmax * exp( (-1.) * df * addq2 ) ) )
     &          * wStreamPowerKgpers3
         end if
         
	return
C
	end

C------------------------------------------------------------------------------ 

      real function streampower_velocity ( slope,tw, h, velocity )
C
C     Finds effective shear stress acting on soil
C     tw/h- width/depth of flow at given time at given node, 
C     tadj - adjustment factor for shear stress
      
	real slope,tw, h, tadj, R
  	if ( h .gt. 0.0001 ) then
	 R = tw * h / ( tw + 2.*h ) 
C      specific weight of water is taken to be 9807 
C      slope factor is approximated as sin of the slope
	 streampower_velocity = 9807. * sin(atan( slope ) ) * R
     &                                                    * velocity  
      else
	 streampower_velocity = 0.
	end if
	return 
	end
C------------------------------------------------------------------------------

      real function probability ( slope, bare, h, q2_per_unit_width )
C
C     Finds the probability of overland flow to concentrate flow based
C     on a multiple logistic regression
C     Equation (13) in Al-Hamdan et al. (in review)
C     prob = numerator/denominator
C     numerator   = exp(-6.397 + 8.335 * slope + 3.252 * bare + 3440 * q2_per_unit_width)
C     denominator = 1+ numerator
C     slope : slope (m m-1)
C     bare  : fraction of bare soil to total area
C     h     : depth of flow
C     q2_per_unit_width : discharge per unit width overland flow
!      real slope,bare,h,q2_per_unit_width, numerator, denominator
      real slope,bare,h,q2_per_unit_width, arg
      double precision numerator, denominator
	if (h .gt. 0.0001) then
!          numerator = exp( -6.397 + 8.335 * slope + 3.252 * bare
!     &                            + 3440. * q2_per_unit_width )
          arg = -6.397 + 8.335 * slope + 3.252 * bare
     &                 + 3440. * q2_per_unit_width
          numerator = dexp( dble( arg ) )
          denominator = 1.d0 + numerator
          probability = sngl( numerator/denominator )
	else
          probability = 0.0
	end if
	return
	end
C------------------------------------------------------------------------------
