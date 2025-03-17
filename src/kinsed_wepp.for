C   Code type: Fortran subroutine + 2 entry points

C   Compiler: Fortran77

C   Initial code programmed by R. Smith for sediment routine was modified 
C   by N. Bulygina to fit WEPP-like equations for plane

C   Date: 12/2004

C   Description:

C     A four-point explicit finite-difference solution of the equation of
C     supply and conservation of a transported material. The solution is
C     explicit since advantage is taken of the kinematic wave solution for the
C     same time step. This routine is called from PLANE modules.
C     This version treats a sediment supply drawn from as
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
C     ql(20)         real         lateral inflow of water ...
C                                 ... for a plane, the average rate over the time
C                                     interval, m/s or ft/s;
C     qil(20)        real         bed or soil inflow of water (m/s or ft/s),
C     a2(20)         real         area at end of interval divided by 
C                                 rill spacing (ams) or depth of overland flow,
C     q2(20)         real         discharge per unit width or ams at end of time interval,
C     wid(20)        real         width of plane
C     rf             real         current rainfall rate (m/s or ft/s),


C     time           real         current simulation time,

C     nx             int          number of spatial nodes,
C     omeg           real         time weighting factor from finite-diff. eqn.
C                                        of calling routine,
C     delx           real         spatial increment (m or ft),
C     slp(20)        real         bottom slope at each spatial node.
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
C     canopy (CA)        real     canopy cover for Ki adjustment (fraction) 

C     Kib (Ki)           real     base iterrill erodibility parameter (kg*s/m^4)
C     Kr   (Kr)          real     rill erodibility parameter (s/m)   

C     tauc   (Tc)        real     critical shear stress (Pa) 
C     intcov (IC)        real     Interrill cover fraction used for interill erodibility
C                                 adjustments  

C     t_adj (Tadj)       real     factor for shear stress adjustment which equals to the ratio of
C                                 Darcy-Weisbach friction factor for the bare soil to the one 
C                                 for the composite surface  

C     ams (SPA)          real     average microtopographic spacing (m or ft),

C     
C------------------------------------------------------------------------------

C   Subroutines/functions:

C     Intrl    (kinsed_plane.for)
C     shrstress         "  
C     Trc               "
C     Dcap              "  
C     vsetl        (kinsed.for) 
C     switch            "
C     getr4        (reader.for)
C     swrt         (writer.for)
C     errxit       (K2RunW.for)
C------------------------------------------------------------------------------


      subroutine kinsed_wepp ( indx, dt, dte, ql, a2, q2, rf, 
     &                 	qup, time)

      use runpars
      use elpars
      use multip
C                                                                     arguments
C------------------------------------------------------------------------------

      integer  indx, nx
      real dt, dte, rf, time
      real, dimension(20) ::  a2, q2, wid, slp
      real, dimension(20) :: ql
      real, dimension(10) :: qup
      integer, dimension(10) :: inflow

      character msg*(*)

C                                                               local variables
C------------------------------------------------------------------------------
C  qr2- discharge in a rill,
C       while q2-discharge per unit width of the plane (per unit rill space for discharge
C       corresponding to each rill space) 
C  ar2- x-sectional area of flow in a rill,
C       while a2-hydraulic depth of flow in the rill space 
C  h2 - depth of flow at given time and given node    
      real, dimension(20) :: qr2, ar2, h2, eros 
      integer :: ierr, i, im, j, jp, n, na, nup, nlat,tstp, redist
      integer, dimension(5):: excd
      integer,save :: nk, ilast, nu
C
	real Tct, Dc_v, IR_v, coef1,coef2,coef3,int_adj,
     &	 deno1, deno2, eterm1, eterm2, eterm, dtrc, dtrc1, dtrc2, 
     &     cs1, cs2, s1, s2, s3, cmax, width, uu, dz, rs
      real, save::  omega, comega, dx, qsp, tsp, xl, qfac, 
     &              intcov, t_adj, Kib, Kr, tauc, rfp, ams, 
     &              canopy, Kibadj, alpha, p, widthp 
      real,save,dimension(5) ::wso, dwso, dwsol, wsin, dwsup, dwsupl
      real,save,dimension(5) :: csub, pros, csub1, csub2, dcsub
      real,save,dimension(20) :: depnet, slope, a1, q1, qil
      real,save,dimension(5,20) :: cs, cl 
      real,save,dimension(20):: cst, clt, erosp
	real cpt, cpmt
	real,save,dimension(20) ::qr1, ar1, h1
      real, dimension(5):: Tci, cp, cpm
C
      integer, save, dimension(5) :: outc
      integer, save, dimension(10,5) :: upc
      
C
      character(LEN=2),dimension(5) :: attr = reshape((/'S1', 'S2', 
     &'S3', 'S4', 'S5'/),(/5,1/))

      character ka_str*1

C------------------------------------------------------------------------------

11	format(a)

C     width of the rectangular channel- doesn't depend on location, 
C     but changes with time      
      q_max = q2(1)
      do j = 1, nk
        if (q2(j) .gt. q_max) q_max = q2(j)
      end do
	width = 1.13 * (q_max * ams)**0.303
	if (width. gt. ams) width = ams

C      Discharge, x-sectional area and depth of flow per rill space
	  do j=1,nk
	    qr2(j) = q2(j) * ams

c         compute rill flow depth h2(j) for rectangular rills. This is an iterative
c         process to solve the uniform flow equation:
c
c         h2(j)=((qr2(j)/alpha/sqrt(slope(j))**(1/(p+1))/width)*
c                (width+2*h2(j))**(p/(p+1))
c
          if (qr2(j).le.1.E-6) then
           h2(j) = 0.0
          else
           uu = (qr2(j)/alpha/sqrt(slope(j))) ** (1./(p+1.)) / width
           h2(j) = 0.2 * (q2(nk) * ams) ** .36
   2       dz = h2(j)
           h2(j) = uu * (width + 2*dz) ** (p/(p+1.))
           if (abs(dz/h2(j)-1.).gt.0.000001) go to 2
          end if

	   ar2(j) = h2(j) * width 
	  end do
      if (indx .ne. ilast) then
C                                             store outgoing concentration from
C                                                   end of last fixed time step
        do n = 1, nps
         call store (outc(n), ilast, cs(n,nk)/(rho(n)*1000.))
	  end do

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
	 if (a1(i) .lt. 1.E-8 .and. qil(i) .gt. 0) then
		cpt = 0.
		do n=1,nps
     	       cp(n) = Kibadj * rfp * qil(i)/(qil(i)+0.5*vsl(n))
	       cpt = cpt + Kibadj * rfp * qil(i)/(qil(i)+0.5*vsl(n)) 
		  end do
	 end if

       if (a1(im) .lt. 1.E-8 .and. qil(im) .gt. 0) then
		cpmt = 0.
		do n=1,nps
             cpm(n) = Kibadj *rfp *qil(im)/(qil(im)+0.5*vsl(n))
	       cpmt = cpmt + Kibadj * rfp * qil(im)/(qil(im)+0.5*vsl(n)) 
		  end do
	 end if
       
C                                            compute erosion or deposition term
C------------------------------------------------------------------------------
	 call Trc_mult( rho,dpd, slope(i), width, h2(i), 
     &   		        t_adj, nps, pros, Tci, Tct) 
	 Dc_v = Dcap ( Kr, tauc, slope(i), width,
     &             h2(i), t_adj ) 
       
       IR_v = Kibadj * rf * qil(i) * ams
	 if (IR_v. lt. 0.) IR_v = 0.
C*******************************************************************************
C      Suppose there is only rill/interrill detachment, not deposition 
C*******************************************************************************          
	 eterm = 0.  

     	 if (shrstress(slope(i),width, h2(i),t_adj).gt.tauc) 
     &     eterm = Dc_v * width
	    		
	 dtrc = eterm + IR_v 		     
C                                                         compute concentration
C------------------------------------------------------------------------------

	 deno = ar2(i)/dt + omega/dx*qr2(i) 
	 if (shrstress(slope(i),width,h2(i),t_adj).ge.tauc) 
     &     deno = deno + Dc_v/Tct * qr2(i)                    
            	                  
	 coef1 = ar1(i)/dt - comega/dx*qr1(i)
	 coef2 = omega/dx * qr2(im)
	 coef3 = comega/dx * qr1(im)
       cst(i) = (cpt * coef1 + cst(im) * coef2 +
     &           cpmt * coef3 + dtrc)/deno     
   
	 if (cst(i). lt. 1.e-15) cst(i) = 1.e-15  
C      Detachment is considered as nonselective process, therefore
C      all detachted sediment on every segment is divided between classes
C      according to there fraction in the soil           
	 do n = 1,nps
	! incoming sediment is averaged over time
	    cs(n,i) = ( cs(n,im) + cl(n,im) ) / 2. + 
     &		      ( cst(i) - (cst(im) + clt(im)) / 2. ) * pros(n)
	    if (cs(n,i). lt. 1.e-15) cs(n,i) = 1.e-15  
	 end do
	 		         
C     eros is in kg/(s*m^2)
       eros(i) = ((cst(i)*qr2(i)+ clt(i)*qr1(i))-
     &        (cst(im)*qr2(im) + clt(im)*qr1(im)))/  
     &        (dx * (width+widthp))       
     
C      were there really only rill and splash detachments ?              
	 if ((Tct*width). lt. (cst(i) * qr2(i))) then   
C                                                         there was deposition 
C-----------------------------------------------------------------------------
C                                                   loop thru particle sizes...
C------------------------------------------------------------------------------
        do n = 1, nps        
C                                                       compute deposition term
C------------------------------------------------------------------------------
		dtrc = 0.5 * vsl(n) * Tci(n)/ qr2(i) * (width)**2   
     &		 + IR_v * pros(n) 		     
C                                                         compute concentration
C------------------------------------------------------------------------------

C --- "a2" should be ar2 in line 362 Mariano Hernandez 05/10/2009---  
C	    deno = a2(i)/dt + omega/dx*qr2(i) +
C     &		   0.5 * vsl(n) * width
C --- Revised
	    deno = ar2(i)/dt + omega/dx*qr2(i) +
     &           0.5 * vsl(n) * width
C --------------------------------------------------------	 
			   
	    coef1 = ar1(i)/dt - comega/dx*qr1(i)
	    coef2 = omega/dx * qr2(im)
	    coef3 = comega/dx * qr1(im)

          cs(n,i) = (cp(n) * coef1 + cs(n,im) * coef2 +
     &	            cpm(n) * coef3 + dtrc)/deno
           if (cs(n,i). lt. 1.e-15) cs(n,i) = 1.e-15

	  end do
         
C                                             ...end loop thru particle classes
C------------------------------------------------------------------------------
C check if sediment for each particle class after deposition region does not 
C exceed maximum allowable one (sediment in + sediment from interrill erosion)
        do n = 1,nps
	    excd(n) = 1  ! 1= particle class load doesn't exceed or equals allowable max
	  end do
  50    s1=0.
	  s2=0.
	  s3=0.
	  redist = 0

	  do n =1, nps
	   s1 = s1 + cs(n,i)
	   cmax = (cs(n,im) * qr2(im) + cl(n,im) * qr1(im)) / 2. + 
     &   	       IR_v * pros(n) * dx 

         if (excd(n). eq. 1.) then
	     if ((cs(n,i) * qr2(i)) .gt.cmax) then
	  	  cs(n,i) = cmax / qr2(i)
	      excd(n) = 0
	      redist = 1
	     else
	      if ((cs(n,i)*qr2(i)). lt. cmax) s3 = s3 + cs(n,i)
	     end if
         end if 
	   s2 = s2 + cs(n,i)
	  end do   
        
	  if (redist. eq. 1) then
         do n = 1,nps
	     if (excd(n).eq.1) then
		 cs(n,i) = cs(n,i) + (s1-s2) * cs(n,i) / s3
	     if (cs(n,i). lt. 1.e-15) cs(n,i) = 1.e-15  
   	     end if
	   end do
	   go to 50
	  end if 

	  cst(i) = 0.
	  do n = 1,nps
	     cst(i) = cst(i) + cs(n,i) 
	  end do
C     eroded sediment per unit rill area        
	  eros(i) = ((cst(i)*qr2(i) + clt(i)*qr1(i))-
     &        (cst(im)*qr2(im) + clt(im)*qr1(im)))/  
     &        (dx * (width+widthp))         	    
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
C      depnet will be multiplied by 1000 before printing out;    
C     for WEPP equations sediment source in the right side of the equation
C     is per unit area of the RILL!, not just unit area 
		 depnet(i) = depnet(i) - 0.5*(eros(i) * width + 
     &	  	         erosp(i) * widthp)* dt/ams/1000.  
	end do
C                                                           ...end spatial loop
C------------------------------------------------------------------------------

C                                            reset variables for next time step
C------------------------------------------------------------------------------
	  rfp = rf                                
        widthp = width
	  do i = 1, nk
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
        dwso(n) = qfac * q2(nk) * cs(n,nk)/(rho(n)*1000.)
	  wsold = wso(n)
        wso(n) = wso(n) + dt2 * (dwso(n) + dwsol(n))
	
 	  dwsol(n) = dwso(n)
        twso = twso + dwso(n) * rho(n)	  
	end do
C                                             check for peak sediment discharge
      if (twso .gt. qsp) then
        qsp = twso
        tsp = time
      end if
C                                                 compute contribution to total
C                                                   inflow at current time step
      do n = 1, nps
        dwsup(n) =  qfac * q2(1) * csub(n)
c        dwsup(n) = dwsup(n) + qfac * q2(1) * csub(n)

        wsin(n) = wsin(n) + dt2 * (dwsup(n) + dwsupl(n))
        dwsupl(n) = dwsup(n)
      end do
    
      return

C------------------------------------------------------------------------------

      entry sedfin_wepp
C                                          finish sediment balance computations
C------------------------------------------------------------------------------

      do n = 1, nps
C                                                  store outgoing concentration
        call store (outc(n), ilast, cs(n,nk)/(rho(n)*1000.))     
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
	  wss = wss + 0.5 * (dwss2 + dwss1) * widthp * dx/1000. 
	     dwss1 = dwss2
      end do
      wss = qfac/ams * wss

C                                                         pass values to writer
C------------------------------------------------------------------------------
      call swrt ( wso, wsi, wse, wss, qsp, tsp, outc, wse_neg, wse_pos )

      return

C------------------------------------------------------------------------------


      entry sed0_wepp ( msg, inflow, nup, nx, omeg, delx, slp, qf )

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
            call old (inv, attr(n), upc(i,n), nc, ierr)
            if (ierr .gt. 0) call errxit
     &      (msg, 'upstream sediment concentration not found')

          end do

        end do
     
C------------------------------------------------------------------------------

      do n = 1, nps
C                                                     obtain blocks for outflow
C                                                         concentration storage
        call new (idu, attr(n), outc(n))
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
           call getr4 ('CA', 0, canopy, ierr)
C                                                                  Canopy cover
	     if (ierr .gt. 0) call errxit 
     &                          (msg, 'Canopy cover (CA) not found')
	     if (canopy .gt. 1.0 .or. canopy .lt. 0.0) call errxit 
     &	                      (msg,'Illegal value for canopy cover') 

	     call getr4 ('CHT', 0, ht, ierr)
C                                                                 Canopy height
	     if (ierr .gt. 0) call errxit
     &                         (msg, 'Canopy height (CH) not found')

		 call getr4 ('IC', 0, intcov, ierr)
C                                                               Interrill cover
		 if (ierr .gt. 0) call errxit 
     &                   	    (msg, 'interrill cover (IC) not found')
 10        format(6f15.8)	  
		 if (intcov .gt. 1.0 .or. intcov .lt. 0.0) call errxit 
     &	                    (msg,'Illegal value for interrill cover') 

        	 call getr4 ('TC', 0, tauc, ierr)
C                                                          critical shear stress
		 if (ierr .gt. 0) call errxit 
     &                          (msg, 'tc (Tc) not found')
		 if (tauc .lt. 0.0) call errxit 
     &	  (msg,'Illegal value for dry bulk density') 
	
	     
	     call getr4 ('TADJ', 0, t_adj, ierr)
C                                                           shear stress factor
		 if (ierr .gt. 0) call errxit 
     &                     (msg, 'Shear stress factor (TADJ) not found')
		 if (t_adj .lt. 0.0) call errxit 
     &	                 (msg,'Illegal value for shear stress factor') 
      
C                                                   get microtopography spacing        
	     call getr4 ('SPA', 0, ams, ierr)

           if (ierr .gt. 0) call errxit
     &        (msg, 'average microtopographic spacing (SPA) not found')
	   
	     call getr4 ('KI', 0, Kib, ierr)
C                                                        Interrill erodibility
	     if (ierr .gt. 0) call errxit (msg, 'Kib (KI) not found')
	      
	     if (Kib .lt. 0.0) call errxit (msg,'Illigal value for Kib')

C                                                     Adjust Ki? (default = Y)
       call getstr( 'KA', 0, ka_str, l, ierr)

       if(ierr .ne. 0 .or. (ka_str .eq. 'y' .or. ka_str .eq. 'Y')) then

C        Ki adjustment (rangland)

C	    Kibadj = Kib * Exp ( -7.0 * (intcov + canopy))

C        Ki adjustment (crop)

         int_adj = exp(-2.5 * intcov)
         cov_adj = 1.- 2.941 * canopy * (1.- exp(-0.34 * ht)) / ht

	   Kibadj = Kib * int_adj * cov_adj

       else
	     
	   Kibadj = Kib

       end if

		 call getr4 ('KR', 0, Kr, ierr)
C                                                              Rill erodibility
		 if (ierr .gt. 0) call errxit 
     &                          (msg, 'Kr (KR) not found')
          Kr = Kr/100000.      
		 if (Kr .lt. 0.0) call errxit 
     &	                      (msg,'Illegal value for Kr')
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

      do i = 1, nk
        depnet(i) = 0.
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
      end do
      widthp = 0. 
      return
C------------------------------------------------------------------------------
      end

C------------------------------------------------------------------------------ 

      real function shrstress ( slope,tw, h, tadj)

C     Finds effective shear stress acting on soil
C     tw/h- width/depth of flow at given time at given node, 
C     tadj - adjustment factor for shear stress
      
	real slope,tw, h, tadj, R
  	if (h. gt. 0.0001) then
	 R = tw * h / (tw + 2.*h) 
C      specific weight of water is taken to be 9807 
C      slope factor is approximated as sin of the slope
	 shrstress = 9807. * sin(atan(slope)) * R * tadj  
      else
	 shrstress = 0.
	end if
	return 
	end
C------------------------------------------------------------------------------
      real function shield(reyn)

      real reyn
      real y(8), r(8), slope, ycr
      integer i
      save
      data y /0.0772, 0.0579, 0.04, 0.035, 0.034, 0.045, 0.055, 0.057/
      data r /1.0, 2.0, 4.0, 8.0, 12.0, 100.0, 400.0, 1000.0/

      if (reyn.lt.r(1)) then
        i = 2
        slope = (alog(y(i))-alog(y(i-1))) / (alog(r(i))-alog(r(i-1)))
        ycr = alog(y(1)) - slope * (alog(r(1))-alog(reyn))
      else if (reyn.gt.r(8)) then
        i = 8
        slope = (alog(y(i))-alog(y(i-1))) / (alog(r(i))-alog(r(i-1)))
        ycr = y(8) + slope * (alog(reyn)-alog(r(8)))
c
      else
c
        do 100 i = 2, 8
c
          if (reyn.ge.r(i-1).and.reyn.le.r(i)) then
            slope = (alog(y(i))-alog(y(i-1))) / (alog(r(i))-
     1          alog(r(i-1)))
            ycr = alog(y(i-1)) + slope * (alog(reyn)-alog(r(i-1)))
            go to 200
          end if
c
  100   continue
c
      end if
c
  200 shield = exp(ycr)
c
      return
      end


C------------------------------------------------------------------------------
      subroutine Trc ( ss, d, slope, tw, h, 
     &    	         tadj, Tci, dlt)

C     Calculates transport capasity using Yalin's equation
C	ss - specific gravity, d- diameter,
C     tw/h- width/depth of flow at given time at given node, 
C     tadj - adjustment factor for shear stress, knu- kinematic viscosity of water
C     Tci - transport capacity for the ith p. class
C     dlt - weight for the future total transport capacity calculations

      real ss, d, slope, tw, h, tadj,
     &	 ts, delta, reyn, Ycr, Y, beta, s, dlt, Tci
      
 10   format(6f15.8)
 11   format (a)	

	ts = shrstress (slope, tw, h, tadj)

	if (ts. gt. 0.) then
	 reyn = sqrt(ts/1000.) * d / 0.000001
C      Ycr - dimensionless critical shear stress from Shield's diagram
	 Ycr = shield(reyn)
	 Y = (ts / 1000.) / ((ss-1.) * 9.81 * d)

	 if ( Y .gt. Ycr) then

	 delta = Y / Ycr -1

	 dlt = delta

	 beta = 2.45 * ss **(-0.4) * Ycr ** 0.5 * delta

	 Tci = 0.635 * ss * d * sqrt (1000.) * sqrt (ts) * delta * 
     &  	  (1.- 1. / beta * alog (1. + beta)) 
      
       else 
	  Tci = 0.
        dlt = 0. 
	 end if
      else
       Tci = 0.
       dlt = 0. 
	end if
	return
	end 
C------------------------------------------------------------------------------

      subroutine Trc_mult ( ss, d, slope, tw, h, tadj,
     &	                     nps, pros, Tci, Tct)

C     Calculates total and distributed between particle classes 
C     transport capasity using Yalin's equation
C	ss - specific gravity, d- diameter,
C     tw/h- width/depth of flow at given time at given node, 
C     tadj - adjustment factor for shear stress, knu- kinematic viscosity of water
C     nps - number of particle classes
C     pros - fraction of each particle class in the soil
C     Tci - transport capacity for ith p. class
C     Tct - total transport capacity  
      real slope, tw, h, tadj,
     &	 ts, Ycr, Y, beta, s, grav, Tct, sumdlt
	real, dimension(5):: Tci, pros, delta, ss, d
	integer nps, i

 10   format(6f15.8)
 11   format (a)	
	sumdlt = 0.
	do i=1,nps 
         call Trc (ss(i),d(i),slope,tw,h,tadj, Tci(i),delta(i))
	   sumdlt = sumdlt + delta(i)
	end do
	
	Tct = 0.
	if (sumdlt. gt. 0.) then
	do i=1,nps
	   Tci(i) = Tci(i) * delta(i) / sumdlt * pros(i) * real(nps)
	   Tct = Tct + Tci(i)
	end do
	end if
      return
	end 
C------------------------------------------------------------------------------
      real function Dcap ( K, tauc, slope, tw, h, tadj )
C     Finds detachment capacity
C     tw/h- width/depth of flow at given time at given node, 
C     tadj - adjustment factor for shear stress
      real K, tauc, grav, slope, tw, h, tadj,ts ,tc
	
	tc=tauc
      
	ts = shrstress (slope, tw, h, tadj )
	Dcap = K * ( ts - tauc )
	 
	return

	end
C------------------------------------------------------------------------------
  