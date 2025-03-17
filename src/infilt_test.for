C This file contains infiltration related subroutines and functions
C
C   Code type: Fortran subroutine with additional entry points

C   Compiler: Lahey Fortran 90  This version with modules.

C   Programmed by: R.E. Smith/C. Unkrich

C   Origin Date: 10/93

C   Last modification: 2/12 by RES
C!!  modifications in constant f=Ks case, see below line 230+ , 02/12/02
C!!  modified 12/02 to add ht and qu as separate parameters
C!!  modified 5/05 to handle th < thr in rectr and deal with small th in flams
C!!  modified 2/12 to allow selection of one of 3 infiltration functions using alf

C   Description:

C     Computes an average infiltration rate for each spatial node
C     using a three-parameter model developed by Parlange et al
C     (1982). Modifed for treating general rain patterns, surface
C     head effect, and wet initial conditions, 6/92. Modified for a
C     two-layer soil profile, 10/93. A value of zero thickness for
C     the upper layer indicates a single layer. A value of zero for
C     the lower capillary parameter signals a constant rate through
C     the lower layer when water reaches the interface. Revised again
C     to handle parallel processing for compound channel sections.
C     Revised 10/01 to correct interactions between effects of CV, 
C     surface relief, and two layer systems under complex rainfalls.
C     NOTE: 99 is assumed the number of the diagnostic output file, if
C     logical DIAG is .T.
C------------------------------------------------------------------------------

C   Entry points:

C     infil0       obtains parameter values from input block and initializes
C                  variables,
C------------------------------------------------------------------------------

C   Arguments:

C     afac(20)     real        in: relative area covered by flowing water
C                              out: relative area to which negative vav is to
C                              be applied.  the rest of the area infiltrates at
C                              rate rfi
C
C     rfi          real        rainfall rate
C
C     dtt          real        time step over which computations apply
C
C     h1(20,k)     real        mean flow depths (m or ft),

C     qu(20,k)     real        flow supply rate (m/s, ft/s),
C                              which includes lateral and upstream in channels
C     vav(20,k)    real        dt-average infiltration rate (m/s or ft/s),
C
C     ichan        int         1 for planes and main channels, 2 for call
C                              from overbank area

C------------------------------------------------------------------------------

C   Input Block (input file labels in parenthesis):

C     sk1 (KS)           real      Ks, upper layer, mm/hr or in/hr,

C     sk2                          Ks, lower layer,

C     g1 (G)             real      capillary drive, upper layer, mm or in,

C     g2                           capillary drive, lower layer,

C     al1 (DI)           real      pore size distribution index, upper layer,

C     al2                          pore size distribution index, lower layer,

C     por1 (PO)          real      porosity, upper layer,

C     por2                         porosity, lower layer,

C     rock1 (RO)         real      volumetric rock fraction, upper layer,

C     rock2                        volumetric rock fraction, lower layer,

C     sint: thin (SA)    real      fraction of saturated pore space, initially

C     cv (CV)            real      coefficient of variation of Ks, {0,1},
C                                  [less than 0.1 treated as 0]
C     zfl (THI)          real      thickness of upper soil layer, m or ft,
C
C GLOBAL:
C     units              int       1 = metric,
C                                  2 = english,

C     nk                 int       number of spatial nodes,

C     msg              char*(*)    identifies the current element in case an error

C                              is encountered,

C     diagn              log       .true. indicates diagnostic output is desired,
C------------------------------------------------------------------------------

C   Subroutines/Functions:

C     errvf        (infilt.for)
C     rectr             "
C     errly             "
C     setg              "
C     th2oth1           "
C     th1oth2           "
C     rkohy             "
C     fpond             "
C     funk              "
C     flams             "
C     fgex              called by fcapv
C!!   fAvdI             "  called in infilt
C!!   fIsots            called by fAvdI
C     vcurf             "
C     vfm               "
C     fcapv             "  gets f(I), fixed or random
C     iter         (miscel.for)
C     getr8        (reader.for)
C     qprw         (writer.   )  send information on soil sat. to writer
C     errxit       (K2RunW.for)
C - - - - - - - - - - - - - - - - - - - - - - - -- -- - - - - - - - --- - - -
C
C     Working Variables
C
C     thiu    effective current initial water content at pulse front
C     fi()    infiltrated water above fmin counted as prewet depth
C     fs()    overall infiltrated depth
C     fv()    infiltrated depth of pulse during r > K
C     fr()    infiltrated water subject to redistribution during hiatus
C     f2l     infiltrated depth at beginning of time step
C     f2n     infiltrated depth after time step
C     f2b     infiltrated depth of main pulse, end of dtt
C     f2p     infiltrated depth of rewet pulse
C     ffl     infiltrated depth above soil interface
C     cumcm   max. capacity of upper layer:  cumcm >= ffl
C     gu      effective G depending on location of wetting front.
C     fkp     net infiltration due to flux of initial water
C     qki     relative flux of initial pulse during rewet
C     ref     effective surface supply rate including surface depth
C     erfi    rainrate equivalent to maintain initial soil water
C     rhc     coverage factor, based on % area covered by flowing water
C               -related to afac
C     tho()    water content at surface
C     zp()     depth of rewet pulse
C     zb()     depth of main wetting pulse
C     vmin     aerial effective final f for infilt at surface
C     vmu      value of vmin used for base or rewet calculation
C     twolay   logical .T. when there are two soil layers
C     subcon   logical .T. when lower layer has Ks2 .le. 2.*Ks1
C     lowet    logical .T. when wetting reaches layer interface
C     topfil   logical .T. when upper layer is saturated.

C------------------------------------------------------------------------------

      subroutine infilt ( afac, rfi, h1, qu, dtt, vav, ichan )

      use multip          ! this module for LaheyF version only
      use itpars
      use elpars
      use runpars

      implicit double precision (a-h, o-z)
      real rfi, dtt, sat
C                                                                     arguments
C------------------------------------------------------------------------------

      logical :: diagn, op

      integer ::  nx, ichan   !, k

      real,dimension(20) :: afac
      real,dimension(20,nchan) :: h1, vav, qu
C                                                               local variables
C------------------------------------------------------------------------------
C
      character msg*(*)
      character*6 trace, source, fmode, set*2
      character(LEN=20) :: calldat
      character*14 notice

      real  :: sumfs, sumv(20)

      logical :: twolay, ifl, topfil(20,2), subcon(2), lowet(20,2)
     &, notify(2), dbg  !

      integer :: ierr, nk, i, j, niter, idu(2)

      dimension fs(20,2), fv(20,2), thiu(20,2), tho(20,2), fr(20,2),
     &          fi(20,2), fkp(20,2), qki(20,2), fcr(20,2), zp(20,2),
     &          zb(20,2), gu(2), cumcm(2), thin(2), ziu(20), t(2),
     &          def1(2), def2(2), por1(2), por2(2)

      common /infil/ twolay(2), zfl(2), f2l, vmu, sk1(2), sk2(2), g1(2),
     &         g2(2), ths1(2), ths2(2), al1(2), al2(2), apu, vmin,
     &      thr1(2), thr2(2), thin2(2), rock1(2), rock2(2), alf, pb1(2),
     &      pb2(2), eps1(2), eps2(2), erfi(2), dtu, trp, i, dbg, intyp

      common /layr/ g2u(2), thm1(2), thm2(2), thi2, rc(2), vmin2(2),
     &              ref, wcf(2), ffl(20,2), sku(20), cv(2), fil, idl, j

C------------------------------------------------------------------------------

      external errvf, errly

C-----------------------------------------------------------------------------
      
      save  fs, fv, thiu, tho, fr, fi, fkp, qki, fcr, zp, zb, 
     &      gu, cumcm, thin, ziu, def1, def2, por1, por2, t, 
     &      nk, jd, notice, sumrf, sumif, topfil, subcon, lowet, 
     &      notify, niter, po, idu, defl, rxv
C------------------------------------------------------------------------------

      i = ichan
C   
      t(i) = t(i) + dtt
      trp = t(i)/60.
C      
      trace = 'infilt'
      intyp = infm    ! transfer from elpars module to common /infil/
C
        do j = 2, nk
C                                                loop thru 2 -> nk locations...
C------------------------------------------------------------------------------
          rfj = rfi   !(j,i)
          rfdt = rfj*dtt
          ht = h1(j,i)
          if(j .eq. jd) sumrf = sumrf + rfdt
          dtu = dtt

          if (sk1(i) .le. 0.) then
C                                           impermeable soil option
            vav(j,i) = 0.
            afac(j) = 1.0
            fmode = 'imperv'
            rxv = 0.

          else if (g1(i) .eq. 0.) then    ! fixed loss option
C   note this option does not now recognise two-layer systems (could be added)
C** following line partially modified 12/02 for steady validation case
            refs = rfj + qu(j,i)      ! constant rate, but
            rst = rfj / sk1(i)              ! varies with rst when CV > 0.1
            vfmf = 1.0
            if(cv(i) .gt. 0.05) vfmf = vfm (rst, cv(i))
            vmin = sk1(i) * vfmf
            if(vmin .lt. rfj) afac(j) = 1.0
            rhc = afac(j)
C!! the following 2 lines are modifications.  
            vav(j,i) = min(vmin,rfj)
            if(refs .gt. 1.e-8 .or. ht .gt. 1.e-6) 
     &         vav(j,i) = min(vmin, refs/afac(j))  
            fmode = 'fix_fc'
            fs(j,i) = fs(j,i) + vav(j,i)*dtu
C                                       need some changes for microtopog
C  normal option
          else
            sku(j) = sk1(i)
            defs = def1(i)
C
            if(lowet(j,i)) then
!              if(subcon(i) .and. fup .eq. 0.) then
              if(subcon(i)) then
C    do flux match when control shifts to lower layer
                sku(j) = sk2(i)
                if(sk2(i) .gt. 0.) then
                  fst = vav(j,i) / sk2(i)
                  if(fst .le. 1) then
                    fup = 0.           !?
                  else
                    uapu = g2u(i) * def2(i)
                    vIeq = fIofs ( intyp, fst, uapu, alf ) !
                    fup = ffl(j,i) * (1. - (vIeq/ffl(j,i)))
                  end if
C                else
                end if
              else
                fup = 0.
              end if
C 
            else
              fup = 0.
            end if
C
            if(subcon(i)) then  ! subsurface control condition:
              if(lowet(j,i)) then
                sku(j) = sk2(i)
                if((ffl(j,i) .ge. cumcm(i)) 
     &              .and. .not. topfil(j,i)) then
                  topfil(j,i) = .true.
                  tho(j,i) = ths1(i)
C  
                  if(j .eq. nk .and. .not. notify(i)) then
                    dum = 1.0
                    izr = 0
                    vfill = cumcm(i)
C                    call qwrt (idl, izr, trace, vfill, dum, dum,
C     &       dum, idu, 0, i, t(i))  ! pass cumcm and t to writer
C                    if(diag .and. j .eq. jd) write(99,*) t(i)
                    notify(i) = .true.
                  end if  
                end if
              else 
                tho(j,i) = ths1(i)
                sku(j) = sk1(i)
              end if
            else
C                     surface layer control conditions
              if (fs(j,i) .ge. cumcm(i)) then   ! was fv and fcr  12/99
C  this assumes twolay, since otherwise cumcm is too large 
                         
                sku(j) = vmin2(i)
                defs = def2(i)
              end if
            end if
C
           if(sku(j) .le. 0.) then  ! impervious lower layer reached
              vav(j,i) = 0.
              rxv = 0.
              fmode = 'imperv'
              go to 100
            end if
            rst = rfj / sku(j)
C
            ref = rfj
C                                         capillary drive includes surface head
            apu = defs * (ht + gu(i))
            if( apu .le. 0.) then
              iunit = 77              ! output file should always be opened
C              inquire (99, opened = op)
C              if(op) iunit = 99
              write(iunit,"(' APU negative in INFILT: '/,i3,(7g12.4))") 
     &           j, g1(i), zfl(i), trp, ht, gu(i), thiu(j,i), defs
              stop    ' negative value of cap. drive in INFILT for this
     &element '
            end if
C  
            rhc = 1.0
C                        factor applicable to spatial distrib. of Ksat:
CCC            if(topfil(j,i)) afac(j) = 1.0   ! new 6/02
            rhv = 0.1 + 0.9*afac(j)   !assumes afac always represents relative coverage
            vfmf = rhv + (1.-rhv) * vfm(rst,cv(i))
C 
            vmin = sku(j) * vfmf      ! new method. sku is found from setg
C                        when surface is covered with water, Keff = Kmean

            viu = fs(j,i) - fup 
            fapt =  favdI( Infm, apu, vmin, viu, alf, dtu)
            fap = fapt/dtu
C!!
            if((qu(j,i)+ht) .gt. 0.) then
C  the following is new criterea [10/01]
              famx = 10.       ! default upper lim.
              if(fs(j,i) .le. 1.e-5) then
                famx = vmin*2.*(1.+dsqrt(0.5*apu/dtu/vmin))  ! new limit
              end if

              if(zp(j,i) .gt. 0.) viu = fv(j,i) - fup
C!!             
C
              if(lowet(j,i)) fap = fap + (cumcm(i) - ffl(j,i))/dtu    ! 10/02
              if(rfj .le. fap ) then !  .and. rfj .le. vav(j,i)
C!! treat runon and lateral inflows on fractional area afac:
                ref = rfj + ht/dtu
C                if(qu(j,i) .le. 0. .and. ht .gt. 0.) ref = ref + ht/dtu ! commented out 5/26/15 by CU
                ref = min(ref,1.01*fap,famx)
                rhc = afac(j) ! reduction when infil from runon
C                                                flow area or during recession
              else
               afac(j) = 1.0
              end if
C              
            else
C 
              fap = -0.
              afac(j) = 1.0  ! this is returned for use with ql in plane
            end if
C!!
            if(diag .and. j .eq. jd) write(99,877) i, trp, ref*po,
     &  fap*po, vmin*po, ziu(j), zb(j,i), fr(j,i), rhc !, rc(i), zb(j,i)
 877  format(/i1,'  time        ref         fpot       vmin        ziu    
     &       zb          fr          rhc'/9g12.4)
C
C   + + + + + + + + + +  +  +  +  +  +  +   +   +   +   +    +    +    +
            if (ref .lt. vmin ) then
C                                                  ponding NOT possible
              vav(j,i) = ref
              refdt = ref*dtu
              if (fv(j,i) .gt. 0.) then
C                                                 begining of redistribution --
C                                                  low rainfall and positive fv
                fr(j,i) = fv(j,i) + fi(j,i) + fr(j,i) 
                fup = 0.
                zp(j,i) = 0.
                fkp(j,i) = 0.
                fv(j,i) = 0.
                fi(j,i) = 0.
                topfil(j,i) = .false.
              end if

              fadd = rhc*refdt + (1.-rhc)*rfdt
              reffa = fadd/dtu
              if (fr(j,i) .gt. 0.) then
C                                                           redistribute during
C                                                                this time step
                fmode = ' redis'
                thor = tho(j,i)
C!!  new trace variable
                write(calldat,'("r",i2,"/",i2,i14)')j,nk,itm
C!!     new trace added to parameter list
                call rectr (fr(j,i), tho(j,i), thin(i), reffa, tho2, 
     &             calldat)
C
                fs(j,i) = fs(j,i) + fadd
                thiu(j,i) = tho(j,i)
C!!  add check  12/02
                deltho = tho(j,i) - thin(i)
                if(deltho .gt. 1.e-6) zb(j,i) =  fs(j,i)/deltho
C!! end add
                ziu(j) = zb(j,i)
C  depth to which rewet pulse must go to find initial conditions:
                if (twolay(i)) then
                  fcr(j,i) = zfl(i)*(ths1(i) - thiu(j,i))
                  if(lowet(j,i)) then
C
                    if(.not. topfil(j,i)) 
     &                   ffl(j,i) = zfl(i)*(tho(j,i) - thin(i))
                    zb(j,i) = zfl(i) + (fs(j,i) - 
     &                       ffl(j,i))/(tho2 - thin2(i))
                  else  
                    ffl(j,i) = ffl(j,i) + fadd
                     if(zb(j,i) .ge. zfl(i)) lowet(j,i) = .true.
                  end if
                end if   ! end twolayer case
                srat = (thiu(j,i) - thr1(i)) / (ths1(i) - thr1(i))

                if (srat .le. 0.) then    ! upper layer relative redistr. flux

                  qki(j,i) = 0.

                else

                  qki(j,i) = srat**eps1(i)

                end if
C                  if(diag .and. j .eq. jd) then
C                    write(99,'(" redist ",8g12.4)') 
C     &    rfj*po, for, fr(j,i), thor, tho(j,i), tho2, qki(j,i)
C                  end if
C
              else    !  no redistributible volume:
C                        prewet -- low rf, no fr or fv:  accumulate fi()
C
                fmode = ' prewt'
                if(fi(j,i) .gt. 0. .and. tho(j,i) .gt. thin(i)) then
                  zz = fi(j,i)/(tho(j,i) - thin(i))
                else
                  tho(j,i) = thin(i) + .001  ! need a kickoff value for RECTR
                  zz = 0.001
                end if
                set = 'pr'
                thiv = thiu(j,i)
C!! twolay test new 12/02
                if(twolay(i) .and. zz .ge. zfl(i)) then
C                  ffl(j,i) = (tho(j,i) - thin(i))*zfl(i)
                  thiv = th2oth1(thiu(j,i))
                  if(subcon(i)) then
                    lowet(j,i) = .true. ! control in lower layer due to init. wetting
                    zb(j,i) = zz
                  end if
                end if
C here the value of sku and apu can change when zz reaches zfl
                if(sku(j) .gt. 0.)
     &          call setg(thiv, zz, fi(j,i), 0.D0, gu(i), i, sku(j),set)
C 
                srat = (thiu(j,i) - thr1(i)) / (ths1(i) - thr1(i))
C
                if (srat .le. 0.) then

                  qki(j,i) = 0.

                else

                  qki(j,i) = srat**eps1(i)   ! this is dimensionless

                end if

                thur = tho(j,i)
C                thiur = thur  ! temporary
                f2t = fi(j,i)  
C!!
                write(calldat,'("p",i2"/",I2,i14)')j,nk,itm
C!! new trace parameter in call list
                call rectr(f2t, thur, thin(i), reffa, tho2, calldat)  ! growing initial pulse
                if(diag .and. thur .le. 0.) then
                  write(99,'(" rerr:",a6,i2,4g12.4)') fmode, j, f2t 
                end if
C
                tho(j,i)= dmin1 (thur, 0.99*ths1(i))
C
                if( .not. twolay(i)) then
                  thiu(j,i) = 0.6 * tho(j,i) + 0.4 * thin(i)
                else
                  thiu(j,i) = tho(j,i)
                end if
                ziu(j) = zz
                dth = thiu(j,i) - thin(i)
                if(f2t .gt. 1.e-8 .and. dth .gt. 0.) ziu(j) = 
     &                  f2t/(thiu(j,i) - thin(i))
C
                if (twolay(i)) fcr(j,i) = zfl(i)*(ths1(i) - thiu(j,i))
                fi(j,i) = fi(j,i) + fadd
                fs(j,i) = fs(j,i) + fadd
                if (subcon(i)) then
                  if( topfil(j,i) ) then
                    ffl(j,i) = zfl(i)*(tho(j,i) - thin(i))
                  else
                    ffl(j,i) = fs(j,i)
                  end if
                end if
C
C        if(diag .and. j .eq. jd) write(99,'(" prew ",i2,7g13.4)')
C     &  jd, f2t, fi(j,i), thiur, thur, zz, ziu(j), reffa*po

              end if
C   low rf case redistributed.  Now cycle:
              go to 100

            end if
C--------------------------------------------------------------------------
C   higher rf case ************************ potential for ponding
            if (ref .gt. vmin .and. vav(j,i) .le. vmin) vav(j,i) = ref
C    
            fbl = fs(j,i) - fup
            f2l = max(fbl, 0.)
            if(sk2(i) .le. 0.) then
              if(zb(j,i) .ge. zfl(i) .or. fr(j,i) .ge. cumcm(i)) then
                lowet(j,i) = .true.
                topfil(j,i) = .true.
                f2n = f2l
                vav(j,i) = 0.
                go to 100
              end if
            end if
C!!  change: compare with small value, not zero
            if (fi(j,i) .gt. 0.001) then
C!! new 12/02: value of storage until prewet pulse is passed:
              capr = (ths1(i) - tho(J,i))*max(zb(j,i),ziu(j))
              if(subcon(i) .and. lowet(j,i)) then
               fi(j,i) = 0.
               fv(j,i) = fs(j,i)
               qki(j,i) = 0.
               if(sk2(i) .le. 0.) fr(j,i) = 0.
C!!
              else ! if(ref*dtu .le. capr) then
               zb(j,i) = ziu(j)
               ziu(j) = 0.
               fr(j,i) = fr(j,i) + fi(j,i)  ! put initial input water into redist.
               fi(j,i) = 0.
               bav = ref
!               f2b = fbl + ref*dtu
               fmode = 'rewet '
              end if
              f2b = fbl + ref*dtu
              go to 75                     ! go to redistr. section
            end if
C                                                        end preponding section
            thiv = thin(i)
            hv = ht                                 ! new 6/02
C!!  add twolay to test  12/02 RES
            if(twolay(i) .and. zb(j,i) .ge. zfl(i)) then
C!! end change
              thiv = thin2(i)
              thi2 = thiv
              if(subcon(i)) then
                lowet(j,i) = .true.     ! RES 6/15
                sku(j) = sk2(i)          ! RES 6/15
              end if
CCC              if(topfil(j,i)) hv = ht + 0.5*zfl(i)  ! proposed improvement 6/02   
            end if
C
            set = 'cd'
            fmode = 'base  '
            call setg (thiv, zb(j,i), fs(j,i), hv, gu(i), i,sku(j), set)
C  apu is also set
            vmb = sku(j)
C                                     set appropriate capillary drive parameter
C!!  add twolay test  12/02 RES
            if(twolay(i) .and. zb(j,i) .gt. zfl(i)) then
C!! end change                                            preponding for base case
              
              if (g2(i) .le. 0.0001) then
C                                                         fixed lower loss rate
                f2b = f2l + vmin2(i) * dtu

                go to 75

              end if
C
            end if
C
C
            If(.not. lowet(j,i)) then  ! single soil or else unfilled:
              if (ref .gt. vmin .and. tho(j,i) .lt. ths1(i)) then
                f2t = fs(j,i) !- fkp(j,i) + (rfj - vmin*qki(j,i))*dtu
                if (f2t .lt. 0.) f2t = 0.
C
                thur = max(tho(j,i),(thin(i)+.01))
                f2o = f2t
C!! new trace parameter for rectr call  18/12/02
                write(calldat,'("b",i2,"/",i2,i14)')j,nk,itm
C                                                               estimate rising
C                                                            surface saturation
                call rectr(f2t, thur, thin(i), ref, tho2, calldat)
C!!
                tho(j,i) = dmin1 (thur, ths1(i))
C
              end if
            else
              f2r = ffl(j,i)/cumcm(i)
              thotry = f2r*(ths1(i) - thin(i)) + thin(i)
              if(thotry .ge. tho(j,i)) tho(j,i) = thotry
            End If
C

            ifl = .false.

            Vol = f2l
            vst = vol/apu  ! dimensionless depth
            vmu = vmin
C      
            call fcapv (infm, Vol, fo, ifl, dum)  ! uses vmu for Ksat
C                    now calculate base infiltration, original
C                                    initial conditions
            fil = fo   !  
C
            If(vst .gt. 0.05 .or. cv(i) .ge. 0.1) then
C!!  new explicit estimator for change in infil depth  18/12/02
              vx =  favdI( Infm, apu, vmin, vol, alf, dtu)
              f2b = f2l + vx
              if( f2l .le. 0.001) then
                vl = f2l + 1.2*vx       !! max value for iteration
              else
                vl = f2l + dtu * vmin * ( apu + f2l) / f2l ! G-A ext for max
              end if
C!!  end change
            else
C
              tsl = vst + dexp(-vst) - 1.  ! parlange/ dimensionless
              tsn = tsl + dtu*vmu/apu
              vsn = tsn + dsqrt(2.*tsn)    ! normalized philip 
              f2b = vsn*apu               ! rescaled infil depth
              vl = 1.5*f2b
            end if
C                                                              vl = upper limit
            f2m = f2l + vmin * dtu
C                                                             f2m = lower limit
C
      if((f2m .lt. 0. .or. f2l .lt. 0. .or. vl .lt. 0.) .and. diag) then
C         write(99,'(" b",i3,4g13.4)') j, fs(j,i), fup, ffl(j,i), fv(j,i)
         stop ' neg stor values in infilt.  see diagno.out' 
       end if
            if (f2b .lt. 1.e-5) f2b = 1.e-5
C    errvf uses vmu for effective Ksat:
C   iterate: try f2b with limits vl(max) and f2m(min): 
            trial = f2b
            source = 'infl  '
            write(source(5:6),'(i2)') j

CLU 2/2013 don't iterate if f2m >= vl

            if (f2m .ge. vl) then

              f2b = f2m

            else

              call iter8 (errvf, f2b, f2m, vl, ierr, source)

              if (ierr .gt. 0) then
                ftav = (f2b - f2l)/dtu
                write(77,788)ierr, j, t(i)/60., f2l, trial, f2m, vl,  
     &                       f2b, vmu*po, ftav*po, vst
 788            format(" infilt err:",2I3,f7.1,8g12.4)
                stop ' error - no convergence (errvf) base'
              else if(j .eq. jd) then
                niter = niter - ierr
              end if

            end if
C
            if(diag .and. j .eq. jd) then
              ftav = (f2b - f2l)/dtu
             write(99,789) j, trial, f2l, f2m, vl, f2b, ftav*po, fil*po,
     &          apu
 789  format(" t, etc.:",I5,8g12.4) !,f7.1
            end if
C  compare capacity with available:
            flm = f2l + ref * dtu
C
            if (f2b .gt. flm) f2b = flm
            bav = (f2b - f2l) / dtu
C
            divz = tho(j,i) - thin(i)
!            if(divz .le. 0.) then
C             write(99,789) j, trial, f2l, f2m, vl, f2b, bav*po, fil*po
C              write(99,'(2g12.5)') tho(j,i), vmu*po
!              stop ' infilt: divz zero '
!            end if
C this is the depth assuming wetting to saturation: (not true after lowet)
!            zb(j,i) = f2b  / divz ! - fkp(j,i)

            IF( divz .gt. 0. ) zb(j,i) = f2b  / divz ! CU 9/11/2013
C
            If(twolay(i)) then
C
              if(subcon(i)) then
C
                if(lowet(j,i)) then
                  if(sk2(i) .gt. 0.) then
                    zb(j,i) = zfl(i) + (f2b )/(thm2(i)-thin2(i)) !- fkp(j,i)
                  else if (zb(j,i) .gt. zfl(i)) then
                    vav(j,i) = 0.
                    f2n = f2l
                    go to 100
                  end if
                else if(zb(j,i) .ge. zfl(i)) then
                  lowet(j,i) = .true.
                  if(sk2(i) .le. 0.) then
                    vav(j,i) = 0.
                    f2n = f2l
                    go to 100
                  end if
                end if
C
              else if (f2b .gt. cumcm(i)) then
C                          ! pulse beyond crust into lowerlayer
                if (g2(i) .le. 0.) then

                  trial = (f2b - cumcm(i)) / (f2b - f2l)*sk2(i)*dtu
                  if (trial .le. cumcm(i)) f2b = cumcm(i)

                else

                  zb(j,i) = zfl(i) + (f2b - cumcm(i)) /
     &                    (thm2(i) - thin2(i))
                end if
C
              end if
C     
            End If
C
  75        continue
C
CCC          if(lowet(i) .and. ref .gt. sk1(i))  ...to be developed
            if (fr(j,i) .gt. 0.) then
C  *************************** rewet or secondary pulse ***********
              zs = zp(j,i) / zb(j,i)
              thiv = thiu(j,i) 
              fmode(5:6) = 'rw'
C
              gb = gu(i)
              vqu = vmin  ! vqu is used for rewet pulse
              if(zs .lt. 1. .and. zp(j,i) .lt. zfl(i)) then
                fupr = 0. 
                vqu = sk1(i)*vfmf
              else if( subcon(i)) then  
                if(zp(j,i) .gt. zfl(i)) then
                  fv(j,i) = f2l - fkp(j,i)
                  fr(j,i) = 0.  ! to bypass rewet infil calcs.
C 
                end if
              end if
              vqki = vqu* qki(j,i)
              f2p = f2l + ref*dtu
C              
              if(bav .lt. vqki .or. vqu .gt. ref ) then !.or. fr(j,i) .le. 0. 
                zp(j,i) = zb(j,i) + 1.e-8
C
              else  if(vqu .lt. ref) then
C                                              create transitional drive term
                if (twolay(i) .and. zp(j,i) .ge. zfl(i)) then    
C                                        ! wetting into second layer
                  if (g2(i) .eq. 0.) then
C                                                       fixed lower loss rate
                    f2n = f2l + vmin2(i) * dtu
                    if(rfj .gt. vmin2(i)) afac(j) = 1.0
                    go to 100

                  end if

                  if(sk2(i) .gt. 0.) then
                    thi2 = th2oth1 (thiu(j,i))
                    thiv = thi2
                  end if
C           
                end if
C
                set = 'rw'
                write(set,'(I2)') j
                call setg (thiv, zp(j,i), fv(j,i), ht, gu(i), i, vmb,
     &                set)
C 
                vqu = vmb * vfmf      ! ? new
                if (abs (vmin - vqu) .le. 0.001 * vqu) then ! vmb is sku found from setg in base calcs

                  fth = (ths1(i) - thiu(j,i)) / (ths1(i) - thin(i))

                  if (fth .lt. 0.15) then
                    gama = 1. - 0.15 * zs
                    if (gama .le. 0.5) gama = 0.5
                    zs = zs * gama
                  end if

                  exfun = 0.

                  if (zs .lt. 0.999 .and. zs .gt. 0.) exfun = zs**1.5
C    effective net drive value
                  ge = gu(i) + (gb - gu(i)) * exfun
                  apu = (ge / gu(i)) * apu
                  if( apu .le. 0.) then
                   iunit = 77   ! output file should always be opened
C                   inquire (99, opened = op)
C                   if(op) iunit = 99
            write(iunit,"(' APU negative in infilt-rw: ',2i3,(6g13.5))") 
     &           j, i, trp, zb(j,i),ht, gb, gu(i), thiv 
                 stop ' negative value of cap. drive in INFILT/rw for th
     &is element '
                  end if

                end if


C Begin edit 08/29/2019 by CU

C                f2l = fv(j,i) - fkp(j,i)  ! - fupr

                if(fv(j,i) .gt. fkp(j,i)) then
                  f2l = fv(j,i) - fkp(j,i)  ! - fupr
                else
                  fkp(j,i) = fv(j,i)
                  f2l = 0.
                end if

C End edit 08/29/2019 by CU

                vmu = vqu * (1. - qki(j,i))
C                                                             
                Vol = f2l
                vls = Vol/apu
                ifl = .false.

                call fcapv (infm, Vol, fo, ifl, dum)  ! uses vmu set just above.
C                               Now make explicit estimate of f2p (vmin > vmu)
                fil = fo
C                 
                if(vls .ge. 0.05 .or. cv(i) .gt. 0.1) then
C!!                    ! replaced with the next line, new estimator  12/02
                  Vx =  favdI(Infm, apu, vmu, vol, alf, dtu)
                  f2p = f2l + Vx         ! predicted
                  if( f2l .le. 0.001) then
                    vl = f2l + 1.2*vx       !! max value for iteration
                  else
                    vl = f2l + dtu * vmin * ( apu + f2l) / f2l ! G-A ext for max
                  end if
C!!
                else
                  tsl = vls + dexp(-vls) - 1.  ! parlange/ dimensionless
                  tsn = tsl + dtu*vmu/apu
                  vns = tsn + dsqrt(2.*tsn)    ! normalized philip
                  f2p = vns*apu               ! rescaled estimated new infil depth
                  vl = 1.5*f2p
                end if
C
                f2m = f2l + vmu * dtu
C iterate: try f2p within limits vl(hi) and f2m(lo): 
C                                                                   rewet pulse
      if((f2m .lt. 0. .or. f2l .lt. 0. .or. vl .lt. 0.) .and. diag) then
        write(99,'(" p",i3,6g13.4)') j, fs(j,i), fupr, ffl(j,i), fv(j,i)
     & , f2l
        end if
C                                                          infiltration depth
                source = 'infr  '
                write(source(5:6),'(i2)') j
                trial = f2p
C
CLU 6/2016 don't iterate if f2m >= vl

                if (f2m .ge. vl) then

                  f2p = f2m

                else

                  call iter8 (errvf, f2p, f2m, vl, ierr, source)
C                                        ervf uses vmu for the effective Ksat
                  if (ierr .gt. 0) then
                    write(77,969) j,ref*po, f2l,f2p,f2m,vl,trial,vns,tsn
                    stop ' error - no convergence (errvf) rewet'
                  end if
                  if (f2p .lt. f2l) then
                  write(77,969) j,ref*po,f2l,f2p,f2m,vmu,vl,bav*po,
     &                          vqki*po
969                 format(i3,8g13.4)
                    stop ' error - infilt/iter: Fnew less than Flast'
                  end if

                end if

                dft = ref * dtu
                flm = f2l + dft - vqki * dtu
                if (flm .lt. f2p) f2p = flm
                rwav = (f2p - f2l) / dtu + vqki
C                                                       rewet flux should never
C                                                        be less than base flux

                if (rwav .lt. bav) f2p = f2l + (bav - vqki) * dtu

                fkp(j,i) = fkp(j,i) + vqki * dtu  !sum of vol moving thru init. pulse
                zp(j,i) = f2p / (thm1(i) - thiu(j,i)) !assumes rewet pulse is saturated
C!! add twolay(i) test 12/02
                if (twolay(i) .and. (zp(j,i) .gt. zfl(i))) then

                  if (g2(i) .eq. 0.) then

                  trial = (f2p - fcr(j,i)) / (f2p - f2l) * sk2(i) * dtu
                    if (trial .le. fcr(j,i)) f2p = fcr(j,i)
                    zp(j,i) = zfl(i) + .01

                  else

                    zp(j,i) = zfl(i) + (f2p - fcr(j,i)) /
     &                      (thm2(i) - thi2)
                  end if
                end if
C 
                if(diag .and. j .eq. jd) then    
        write(99,790) j, vl,f2p,fv(j,i), ht, ref*po,
     &  rwav*po, apu, famx*po
 790  format(" t, rew:",I5,8g12.4)
                end if
              end if    

              if(zp(j,i) .ge. zb(j,i)) then
C                                                           rewet pulse catches
C                                                                original pulse
                fr(j,i) = 0.
                fkp(j,i) = 0.

C Begin edit 08/29/2019 by CU

C                f2l = fs(j,i) - fup

                if(fs(j,i) .gt. fup) then
                  f2l = fs(j,i) - fup
                else
                  fup = fs(j,i)
                  f2l = 0.
                end if

C End edit 08/29/2019 by CU

                thiu(j,i) = thin(i)
                f2n = f2b
                qki(j,i) = 0.
                zp(j,i) = 0.
                fcr(j,i) = cumcm(i)

              else

                f2n = f2p

              end if
C end rewet pulse section
            else 
C                            no rewet pulse
              zp(j,i) = 0.
              f2n = f2b

            end if
C                                  For two-layer system, !  use "excess" to
C                                  fill upper layer until saturated:
            rxv = 0.
            vavu = (f2n - f2l) / dtu + vmin * qki(j,i)
            if(vavu .lt. 0.) then
              write(77,*)' vavu neg ',j,qki(j,i),vmin,f2n,f2l,dtu,f2b,f2p
              stop ' vavu neg at 932 in infil'
            end if 
C   
            vav(j,i) = vavu
            if(subcon(i)) then
              if( lowet(j,i) .and. .not. topfil(j,i)) then
C      fill upper layer storage from infil into upper layer
                defl = cumcm(i) - ffl(j,i)
                rxv = vavu * dtu
                if(rxv .le. defl) then
                  ffl(j,i) = ffl(j,i) + rxv
                else   ! now filled upper layer
                  ffl(j,i) = cumcm(i) + 1.e-8
                  topfil(j,i) = .true.
                  tho(j,i) = ths1(i)
                end if
              else if(.not. lowet(j,i)) then  ! 
                rfin = min(rfdt,(vavu*dtu))
                ffl(j,i) = ffl(j,i) + rhc*(f2n-f2l) + (1.-rhc)*rfin
              end if
C                     ! single layer system or else top already full
            end if
C                     weight increase of storage based on rhc
       rmb = 0.   ! temporary: what is this??
      fv(j,i) = rxv + (fup +f2l) + rhc*(f2n -f2l) + fkp(j,i) + 
     & (1.-rhc)*(rfdt - rmb)
      fs(j,i) = rxv + (fup+fbl) + rhc*(f2n -f2l) + (1.-rhc)*(rfdt-rmb) 
     & + vmin*qki(j,i)*dtu

          IF( fs(j,i) < 0. ) fs(j,i) = 0.  ! CU 9/4/2019

          end if
C          
100       continue

          if(sk1(i) .gt. 0.) then
            if (diag .and. j .eq. jd) then   ! put any diagnostics here:
C              sumif = sumif + vav(j,1)*dtu
C              rmu = ref*po
C              if(twolay(i)) then
C                write(99,198) 
C     & lowet(j,i),topfil(j,i),  qki(j,i), tho(j,i), thiu(j,i), defl, rxv
C              else if(g1(i) .gt. 0.) then
C                write(99,'(" tho/thiu: ",2g12.4)') tho(j,i), thiu(j,i)
C              end if
C              write(99,199)
C              write(99,'(f6.1,2i4,a6,8g12.4)') trp,jd,niter,fmode,
C     &  rfj*po,vav(j,1)*po,vmin*po,fs(j,i),fv(j,i),ffl(j,i),afac(2) ! , sumrf, sumif
        if(diag .and. sk1(i) .gt. 0.) 
     &        write(99,'(" Tot fs:"/,(8g12.4))') (fs(k,1),k=2,nk)
C
            end if
          end if

        end do
C------------------------------------------------------------------------------
C                                                        ...end 1 -> nk loop
        if(diag) then
          sumfs = 0.
          nkm = nk - 1
          do j = 2,nk
            sumv(j) = sumv(j) + vav(j,i) * dtu  ! not used for now
            sumfs = sumfs +  fs(j,i) 
          end do
          sumfs = sumfs/real(nkm,4)  !*qbal(1)   vol = area * av.depth
          write(99,'(" net ave. inf:",g12.4," m ")') sumfs
        end if
  198 format(" lowet topfil    qki      tho         thiu         defl      
     &    rxv"/3x,L2,4x,L2,5g12.4)
  199  format(" infl:t  j iter  mode    rf         f          Ks(eff)      
     &   Fs         Fv          F(1)        alph")
C                                                        
      if(itm .ge. itlim ) then
        sumtab(ltab)%flowr = 0.
        if(diag .and. sk1(i) .gt. 0.) 
     &         write(99,'(" Tot fs:"/,(8g12.4))') (fs(j,1),j=2,nk)
        nkm = nk-1
        if(subcon(i) .and. sk2(i) .gt. 0.) then
          finlot = 0
C       
          do j = 2, nk 
            finlo = 0.
            if(lowet(j,i)) finlo = max((fs(j,i) - ffl(j,i)),0.)
            finlot = finlot + finlo
          end do

          finlom = finlot/real(nkm,4)
C
          if(finlom .gt. 1.e-6) then
            jmns = -1
C           write(99,'(" finlo ",3g13.4)') finlot, fs(2,1), ffl(2,1)
C            call qwrt (idl, jmns, trace, finlom, dum, dum,
C     &       dum, idu, 0, i, t(i))  ! pass cumcm and t to writer
          end if
C
          ni = ltab + i - 1
          sumtab(ni)%flowr = finlom
          if(diag)   write(99,'(" fl:",8g12.4)') (ffl(j,1),j=2,nk)
        end if
C here one could recalc and get thend = (fs - finlo)/zfl
      end if
C
C------------------------------------------------------------------------------
      return

C------------------------------------------------------------------------------


      entry infil0 (iv, ki, nx, msg, diagn, sat)
C                                                    get / initialize variables
      trace = 'infl0 '
      t = (/0.,0./)
      idl = iv
      idu(ki) = iv
      niter = 0
      rsoil = .true.   ! designates real soil: i.e. not constant or 0 Ks
      notice = 'element       '
      write(notice(9:14),'(I6)') idl
      notify = .false.
      sumrf = 0.
      sumif = 0.
C
      diag = diagn
      dbg = diagn

      if (ki .eq. 1) then              ! plane or main channel.  2 for overbank

        nk = nx

        if (units .eq. 1) then

 !          conv = 1000.
          con2 = 1.

        else

 !          conv = 12.
          con2 = 3.28

        end if
C                                                                values assumed
C        smax = 0.95
        po = conv*3600. ! printout coefficient to get mm/h(in./hr) from m/sec (ft/sec)
C        alf = 0.8        

      end if
C 
      jd = 2    ! node for which diagnostics will be written
C
      if(jd .gt. nk) jd = nk
      i = ki    ! this is nchan for initialization
      zfl(i) = 0.
      nvi = ltab + i - 1
C                                                             Ks of upper layer
C                                                              (mm/hr or in/hr)
      call getr8 ('KE', 1, sk1(i), ierr)

      if (ierr .gt. 0) then
        call getr8 ('KS', 1, sk1(i), ierr)
        if (ierr .gt. 0) sk1(i) = 0.
      end if
C                                                           impervious surface?
      if (sk1(i) .eq. 0.) then
        rsoil = .false.
        sumtab(nvi)%twola = .false.
        sumtab(nvi)%thst = 0.
        return
      end if

      sk1(i) = rmks * sk1(i) / (3600. * conv)          ! meters/sec
C                                                          capillary suction of
C                                                        upper layer (mm or in)
      call getr8 ('G', 1, g1(i), ierr)

      if (ierr .gt. 0) g1(i) = 0.
      if(rmg .le. 0.) stop ' zero multiplier'

      g1(i) = rmg * g1(i) / conv                       ! convert to m or ft

      if (g1(i) .gt. 0.) then
C                                                         pore size dist. index
C                                                                of upper layer
        call getr8 ('DI', 1, al1(i), ierr)

        if (ierr .gt. 0) call errxit
     &  (notice, 'pore size dist. index (DI, upper layer) not found')
        if(al1(i) .gt. 1.5) call errxit(notice,' pore size dist. index c
     &annot be > 1.5')

        call getr8 ('PO', 1, por1(i), ierr)
C                                                       porosity of upper layer
        if (ierr .gt. 0) call errxit
     &  (notice, 'error - porosity (PO, upper layer) not found')

        call getr8 ('RO', 1, rock1(i), ierr)
C                                                             vol. rock content
C                                                                of upper layer
        if (ierr .gt. 0) rock1(1) = 0.

C                                                            maximum saturation
        call getr8('SM', 1, smax, ierr)
C
        if(ierr .gt. 0) then
          smax = 0.95
        end if

        ths1(i) = por1(i) * smax
C
        call getr8('SR', 1, ser1, ierr)
C
        if(ierr .gt. 0) then
          ser1 = -1.
        end if

        if(sat .lt. 0.) then       ! init. sat not found in rain file
C
          call getr8 ('SA', 0, sint, ierr)
C                                                      initial saturation
          if (ierr .gt. 0) then
C        
            call errxit
     &    (notice, 'initial soil saturation (SA) not specified')
          end if
          
        else
          sint = sat
        end if

        sint = sint * rmsat
          
        if (sint .ge. smax) g1(i) = 0.

        call getr8 ('CV', 0, cv(i), ierr)
C                                                     coeff. of variation of Ks
        if (ierr .gt. 0) cv(i) = 0.

        cv(i) = cv(i) * rmcv    ! CV multiplier, default 1.0; from module 'multip'
C new section to allow variable alf, within range 0.05 to 0.95
C
        call getr8( 'ALF', 0, alfi, ierr)
        if( ierr .gt. 0) then
          alf = 0.8
          infm = 3
        else
          if( alfi .ge. 0.96) then
            alf = 1.0
            infm = 2
          else if( alfi .le. 0.04) then
            alf = 0.
            infm = 1
          else
            alf = alfi
            infm = 3
          end if
        end if
        alfinf = alf
C
C end new section        
      else          ! no soil - impervious
        rsoil = .false.
        ths1(i) = .4                               ! for g1 = 0.
        sint = 0.2
        thin(i) = .16
        cv(i) = 0.
C   can return here?
      end if


      call getr8 ('THI', 0, zfl(i), ierr)
C                                                            thickness of upper
C                                                              layer (mm or in)
      if (ierr .gt. 0) zfl(i) = 0.
C!!  next two moved to here 3/05  default values     
      g2(i) = 0.
      g2u(i) = 0.
 !
      if (zfl(i) .le. 0.) then
C                                                                  single layer
C!! for possible diagn print
        thin2(i) = 0.
        ths2(i) = 0.
        twolay(i) = .false.
        subcon(i) = .false.
        sumtab(nvi)%twola = .false.
      else
C                                                                    two layers
        zfl(i) = zfl(i) / conv ! now in meters/ft
        if((units .eq. 1 .and. zfl(i) .lt. 1.e-2) .or.  ! 1 cm limit
     &       units .ne. 1 .and. zfl(i) .lt. .05) 
     &     call errxit(notice,' upper soil layer too thin for use: check 
     & units ')
C
        twolay(i) = .true.
        depmax = max(zfl(i),depmax)    ! find watershed max val.ue of zfl
        sumtab(nvi)%twola = .true.

        call getr8 ('KE', 2, sk2(i), ierr)

        if (ierr .gt. 0) then
          call getr8 ('KS', 2, sk2(i), ierr)
          if (ierr .gt. 0) sk2(i) = 0.
        end if

C!!    changed 3/05:                                                         Ks of lower layer
C!!        if (sk2(i) .eq. 0.) then
C                                                        impervious lower layer
C!!          g2(i) = 0.01   ! meters
        if ( sk2(i) .gt. 0.00001) then
C!!        else
C                                                                      pervious
          sk2(i) = rmks * sk2(i) / (3600. * conv)

          call getr8 ('G', 2, g2(i), ierr)
C                                                capillary drive of lower layer
          if (ierr .gt. 0) g2(i) = 0.
          g2(i) = rmg * g2(i) / conv


          if (g2(i) .gt. 0.) then
C                                                         pore size dist. index
C                                                                of lower layer
            call getr8 ('DI', 2, al2(i), ierr)

            if (ierr .gt. 0) call errxit
     &     (notice, 'pore size dist. index (DI, lower layer) not found')
            if(al2(i) .gt. 1.5) call errxit
     &      (notice,' pore size dist. index cannot be >1.5')

            call getr8 ('PO', 2, por2(i), ierr)
C                                                       porosity of lower layer
            if (ierr .gt. 0) call errxit
     &      (notice, 'porosity (PO), lower layer) not found')
C                                                       residual sat of lower layer
            call getr8 ('SR', 2, ser2, ierr)
            if(ierr .gt. 0) ser2 = -1.

            call getr8 ('RO', 2, rock2(i), ierr)
C                                                      vol. rock of lower layer
              if (ierr .gt. 0) rock2(i) = 0.
              ths2(i) = smax * por2(i)

          else if (sk2(i) .ge. sk1(i)) then

            call errxit
     &      (notice, 'fixed lower seepage rate > Ksat of upper soil')

          end if
        end if
      end if
C                                                       initial calculations...
C------------------------------------------------------------------------------

C                                                                values assumed
      if (g1(i) .gt. 0.) then
C!                  estimating function, g in m or ft
        if(ser1 .lt. 0.) ser1 = 0.1 + 0.3/(1.+(con2/g1(i))**4.)**.25  
        thr1(i) = ser1*por1(i)             ! could be *ths1(i)
        dth1 = ths1(i) - thr1(i)
        thin(i) = sint * dth1 + thr1(i)  ! this is more correct, but:
C        thin(i) = sint * por1(i)   ! Carl's version  9-3-2009
        if(thin(i) .lt. thr1(i)) thin(i) = thr1(i) + 0.01 ! revised 9-3-2009 by CLU
!C
        sumtab(nvi)%thst = thin(i)
C
        def1(i) = (ths1(i) - thin(i)) * (1. - rock1(i))
        if(al1(i) .le. 0.) al1(i) = 0.2  ! insert default value for lambda if necessary
        if(al1(i) .gt. 1.5) al1(i) = 1.5
C  estimate a value for Pb:
        pb1(i) = g1(i) * (2. + 5. * al1(i)) / (3.4 + 3. * al1(i))  ! rev 4/02
        eps1(i) = (2. + 3. * al1(i)) / al1(i)
        erfi(i) = sk1(i)*sint**eps1(i)
C        write(77,'(" rfi ",g12.4)') erfi(i)
        thwilt = thr1(i) + dth1*(pb1(i)/153.)**al1(i)  ! 15 bar, wilting point
        thfc = thr1(i) + dth1*(pb1(i)/3.33)**al1(i)    ! est. "field capcy" at 1/3 bar
C report values
        pbv1 = pb1(i)
        thru1 = thr1(i)
        thsv1 = ths1(i)
        vks1 = sk1(i)*3600.*conv
        zlay = zfl(i)
C
        if(twolay(i)) then
C                     this is max water held in upper layer at field capacy:
          sumtab(nvi)%flcap = zfl(i)*(thfc - thin(i)) * (1.-rock1(i))
C
        end if
        thm1(i) = ths1(i)

      end if

      if(twolay(i) .and. g2(i) .gt. 0.) then

        pb2(i) = g2(i) * (2. + 5. * al2(i)) / (3.4 + 3. * al2(i))  ! rev 4/02
C             ! estimating function, g in m
        if(ser2 .lt. 0.) ser2 = 0.1 + 0.3/(1.+(con2/g2(i))**4.)**.25  
        thr2(i) = por2(i)*ser2  
        thin2(i) = th2oth1 (thin(i))
        def2(i) = (ths2(i) - thin2(i)) * (1. - rock2(i))
        if(al2(i) .le. 0.) al2(i) = 0.2  ! insert default value for lambda if necessary
        if(al2(i) .gt. 1.5) al2(i) = 1.5
        eps2(i) = (2. + 3. * al2(i)) / al2(i)
C   output write values
        pbv2 = pb2(i)
        thru2 = thr2(i)
        thsv2 = ths2(i)
        vks2 = sk2(i)*3600.*conv
        thinu2 = thin2(i)
C!!  new lines 3/05
      else if (g2(i) .le. 0.) then  ! imperv or fixed Ks2
        thr2(i) = 0.
        def2(i) = 0.
        pb2(i) = 0.

      end if
C
      gu(i) = g1(i)
      cumcm(i) = 1000.
      thm1(i) = ths1(i)

      if (twolay(i)) then
C                                                        analyze flow potential
        cumcm(i) = zfl(i) * (ths1(i) - thin(i))
C
        sumtab(nvi)%pored = zfl(i)*def1(i)*conv   !zfl(i)*por1(i
C
        vmin2(i) = sk2(i)

        if(g2(i) .gt. 0.) then

          hc = zfl(i) * (1. - sk2(i) / sk1(i))   ! can be large, negative
          thm2(i) = ths2(i)

          if (sk1(i) .lt. sk2(i)) then
C !!   add 1/03
            vmu = sk1(i)
            defl = 0.
            rxv = 0.   ! also named in save statement
C!!                                                            get vmin and thm's
C                                                              for crusted layered case
            rc(i) = sk1(i) * (1. + alf / (dexp (alf * def1(i)*
     &              zfl(i) / g1(i)) - 1.))  !  critical flux
C!! abs() added
            hc = abs(hc)
            if (hc .gt. pb2(i))  then  ! 

             source = 'errlay'
             call iter8 (errly, hc, pb2(i), 100. * pb2(i), ierr, source)
C    this iteration finds effective K for crusted case:  vmu
              if (ierr .gt. 0) stop ' error - no convergence (errly)'

              vmin2(i) = vmu   ! effective composite K for crust
              thm2(i) = thr2(i) + (ths2(i) - thr2(i)) *
     &                  (pb2(i) / hc)**al2(i)

              if (hc .gt. pb1(i)) thm1(i) = thr1(i) + (ths1(i) -
     &                            thr1(i)) * (pb1(i) / hc)**al1(i)

              g2u(i)= pb2(i)*flams(al2(i),thm2(i),ths2(i),thr2(i),trace)

            else
C
              g2u(i) = g2(i) - hc           

            end if

            rat = sk2(i) / sk1(i)
            bf = (dlog (vmin2(i) / sk1(i)) + 1.) / (dlog (rat) + 1.)

            if (bf. lt. 1.0) then

              wcf(i) = bf / (1. - bf)

              if (wcf(i) .gt. 2. * rat) wcf(i) = 2. * rat

            else

              wcf(i) = 2. * rat

            end if

          else  ! sk1 > sk2
C                                                       restrictive lower layer
            g2u(i) = g2(i) !  - hc
            vmu = sk2(i)
            if(sk1(i) .gt. 2.*sk2(i)) subcon(i) = .true.
C
          end if
C  fixed Ks lower bound
        else
C!!  new line
          g2u(i) = g2(i)
          if(sk1(i) .gt. 2.*sk2(i)) subcon(i) = .true.
        end if

      else
C                                                               no second layer
        vmu = sk1(i)
        thm1(i) = ths1(i)
        zfl(i) = 1000.
        cumcm(i) = 1000.
        sk2(i) = sk1(i)

      end if
C
      thinu = thin(i)
      fsz = 0.
      satin = sint
      
      If(g1(i) .gt. 0.) then
        satfc = (thfc-thr1(i))/dth1
        satwl = (thwilt-thr1(i))/dth1
      end if

      do j = 1, nk
        sku(j) = sk1(i)
        fs(j,i) = fsz
        thiu(j,i) = thinu
        tho(j,i) = thinu
        fr(j,i) = 0.
        fi(j,i) = 0.
        qki(j,i) = 0.
        fkp(j,i) = 0.
        ffl(j,i) = fsz
        fcr(j,i) = cumcm(i)
        zp(j,i) = 0.
        zb(j,i) = 0.
        ziu(j) = 0.
        fv(j,i) = 0.
        topfil(j,i) = .false.
        lowet(j,i) = .false.

      end do
C
      if(diagn) then
        write(99,999) vmu*po, cumcm(i), g1(i), g2(i), g2u(i),gu(i),cv(i)
        if(rsoil) write(99,998) thin(i),thin2(i),ths1(i),ths2(i)
     &   ,vmin2(i), rc(i)
      end if
 999  format(2x,' vmu          cumc         G1(m)       G2         G2U 
     &        GU         CV(K)'/8g12.4)
 998  format(2x,' thi          thi(2)      ths1        ths2'/8g12.4)
      return

      end
C------------------------------------------------------------------------------


      double precision function th2oth1 (th)

C   Computes saturation in SUBso7il to match head in SURFACE soil

      implicit double precision (a-h, o-z)

      common /infil/ twolay(2), zfl(2), f2l, vmu, sk1(2), sk2(2), g1(2),
     &         g2(2), ths1(2), ths2(2), al1(2), al2(2), apu, vmin,
     &      thr1(2), thr2(2), thin2(2), rock1(2), rock2(2), alf, pb1(2),
     &      pb2(2), eps1(2), eps2(2), erfi(2), dtu, trp, i, dbg, intyp

      logical :: twolay, dbg

      integer i

      thrto = (th - thr1(i)) / (ths1(i) - thr1(i))
      if (thrto .le. 0.) then
        thrto = 0.
        h1 = 1.e5
      else
        h1 = pb1(i)/thrto**(1./al1(i))
      end if  
      if(h1 .gt. pb2(i)) then
        th2oth1 = thr2(i) + (ths2(i)-thr2(i))*(pb2(i)/h1)**al2(i)
      else
       th2oth1 = 0.99*ths2(i)
      end if
      return

      end
C------------------------------------------------------------------------------


      double precision function th1oth2 (th)

C   Computes saturation in SURFACE soil to match head in SUBsoil

      implicit double precision (a-h, o-z)

      common /infil/ twolay(2), zfl(2), f2l, vmu, sk1(2), sk2(2), g1(2),
     &         g2(2), ths1(2), ths2(2), al1(2), al2(2), apu, vmin,
     &      thr1(2), thr2(2), thin2(2), rock1(2), rock2(2), alf, pb1(2),
     &      pb2(2), eps1(2), eps2(2), erfi(2), dtu, trp, i, dbg, intyp

      logical twolay, dbg

      integer i

      thrto = (th - thr2(i)) / (ths2(i) - thr2(i))
      if(thrto .le. 0.) then
        thrto = 0.
        h2 = 1.e5
      else
        h2 = pb2(i)/thrto**(1./al2(i))
      end if
      if (h2 .gt. pb1(i)) then
        th1oth2 = thr1(i) + (ths1(i)-thr1(i))*(pb1(i)/h2)**al1(i)
      else
        th1oth2 = 0.99*ths1(i)
      end if
      return

      end
C------------------------------------------------------------------------------

      subroutine errvf ( F2n, ferr, dferr)

      implicit double precision (a-h, o-z)

C   Computes the residual and derivative at F2 of the error function form of
C     the 3-parameter infil equation for infiltrated volume as a f(I). 
C     This is an external function for ITER()
C   Cum. Infil. and infiltrability fil at beginning of step, F2l, are
C      assumed known.
C   fv2 is infiltrability at end of step
C   fil is infiltrability at beginning
C
      common /infil/ twolay(2), zfl(2), f2l, vmu, sk1(2), sk2(2), g1(2),
     &         g2(2), ths1(2), ths2(2), al1(2), al2(2), apu, vmin,
     &      thr1(2), thr2(2), thin2(2), rock1(2), rock2(2), alf, pb1(2),
     &      pb2(2), eps1(2), eps2(2), erfi(2), dtu, trp, i, dbg, intyp

      common /layr/ g2u(2), thm1(2), thm2(2), thi2, rc(2), vmin2(2),
     &              ref, wcf(2), ffl(20,2), sku(20), cv(2), fil, idl, j

      logical :: twolay, ifn, dbg

      integer :: i, idl, j, intyp

      data ifn/.true./
C
      if(cv(i) .ge. 0.1) then
C                        ! gets fv2 with spatial capacity function:
        call fcapv (Intyp, F2n, fv2, ifn, dfdi) 
        ferr = F2n - F2l - 0.5 * dtu * (fv2 + fil)
        dferr = 1. - 0.5 * dtu * dfdi
C
      else  ! simpler and more precise for no spatial variance case:
C
        select case (intyp)
          case (1)  ! G-A
            ferr = F2n - F2L - vmin*dtu - 
     &          apu*dlog((F2n + apu)/(F2L + apu))
            dferr = 1.0 - apu/(F2n + apu)
          case (2)  ! S-P
            xex2 = 1. - dexp(-F2n/Apu)
            xex1 = 1. - dexp(-F2L/Apu)
            ferr = F2n - F2L - vmin*dtu - apu*(xex2 - xex1)
            dferr = xex2
          case (3)
C!! new checks for overflo:  12/02
           xtermL = alf*F2L/apu
           xtermN = alf*F2N/apu
           if(xtermL .gt. 80. .or. xtermN .gt. 80.) then
             ferr = F2n - F2L - vmin*dtu
             dferr = 1.0
             if(dbg) write(77,99) idl, j, trp, F2L, apu, F2n
           else
             FG2 = dexp(xtermN)          ! substitute xterms
             ratl = (FG2 - 1. + alf)/(dexp(xtermL) - 1. + alf)
             if(ratl .le. 0. .or. F2l .lt. 0.) then
               write(77,99) idl, j, ratl, F2l, apu, F2n
  99  format(' neg. F2L or logf invalid in errvf: ',2i4,4g14.4)
               stop ' invalid argument for logf '
             end if
             ferr = F2n - F2l - apu*dlog(ratl) - vmin*dtu*(1.- alf)
             dferr = 1.0 - alf*FG2/(FG2 - 1. + alf)
             if( dbg .and. j .eq. 2) then
               write(99,98) trp, ratl, f2n, xtermn, ferr, dferr
  98  format( 8g13.5)
             end if
           end if
C
        end select
      end if
      return
C
      end
C------------------------------------------------------------------------------
C!! new dummy variable (origin) added 12/02  RES
      subroutine rectr (fv, th, thi, rf, th2, origin)

      implicit double precision (a-h, o-z)

C   Computes the redistribution (reduction in saturation) of wetted soil pulse,
C   which is assumed to elongate with time, according to Smith, et al,
C   WRR (1992), modified for a two-layer profile.

      common /infil/ twolay(2), zfl(2), f2l, vmu, sk1(2), sk2(2), g1(2),
     &         g2(2), ths1(2), ths2(2), al1(2), al2(2), apu, vmin,
     &      thr1(2), thr2(2), thin2(2), rock1(2), rock2(2), alf, pb1(2),
     &      pb2(2), eps1(2), eps2(2), erfi(2), dtu, trp, i, dbg, intyp

C!! diagn added 12/02
      logical :: twolay, op, dbg

      integer           :: i, nloop, ncut             !, md
      character(LEN=4)  :: note, look
      character(LEN=6)  :: trace
      character(LEN=20) :: origin
      character(LEN=14) :: notes = 'element       '

      trace = 'rect  '
      note = 'norm'
      look = 'norm'
      write(notes(9:14),'(I6)') i
      il = 1
      fl = fv
      thu = th
      if (th .le. thi .and. rf .lt. erfi(i)) then
        this = thi
        go to 90

      end if

      zn = fv / (th - thi + .001)

      if (twolay(i) .and. zn .ge. zfl(i)) then                ! two-layer case

        if (g2(i) .eq. 0.) then
C                                               fixed flux - lower soil
          dfl = sk2(i) * dtu
          fl = fv - dfl
          this = th - dfl / zfl(i)

          if (this .lt. thi) then

            dfl = 0.
            this = thi

          end if

          go to 90

        end if

        il = 2
        th1 = th
C        phi = ths2(i) - thr2(i)
        alm = al2(i)
        ths = ths2(i)
        thr = thr2(i)
        pb = pb2(i)
        thu = th2oth1 (th)
        if(thu .lt. thr) stop ' thu<thr in 2lay '
        thit = thin2(i)

      else
C                                            redistribute in upper soil
        il = 1
        alm = al1(i)
        pb = pb1(i)
        ths = ths1(i)
        thr = thr1(i)
        if( thu .lt. thr) stop 'thu<thr in 1lay'
C        phi = ths1(i) - thr1(i)
        thit = thi

      end if

      fac = 2.
      if (rf .gt. 1.E-4) fac = 1.5
      if(rf .lt. -1.E-4) fac = 3.
C
      thwt = 0.99 * thu    ! corrected 11/99
C
      trace(5:6) = 'i '
      write(trace(6:6),'(i1)') il

      fog = flams (alm, thwt, ths, thr, trace)    

      gup = pb * fog
      if(thwt .gt. ths .or. thwt .lt. thr) then
        iunit = 77   ! output file should always be opened
C        inquire (99, opened = op)
C        if(op) iunit = 99
        write(iunit,'(" rect bounds ",i3,5g13.4)') 
     &    il, trp, thr, thu,ths, gup 
      end if

 
 47   continue
      this = thu  
      t = 0.
      odtt = dtu
      nloop = 0
      ncut = 0
      if(dbg) note = 'diag'
C get starting time step:
C  G for water content this:
      avgz = flams (alm, this, ths, thr, trace) * pb            !flams
      fk0 = funk (fl, rf, this, avgz, fac, thit, il, note)      ! funk
      note = 'norm'

      IF( fk0 .NE. 0. )THEN ! CU 9/11/2013
        dtz = 0.002/abs(fk0)
      ELSE                  ! CU 9/11/2013
        dtz = dtu           ! CU 9/11/2013
      END IF                ! CU 9/11/2013

      th2 = this -  dtz * fk0
      flz = fl + rf*dtz
      fkn = funk (flz, rf, th2, avgz, fac, thit, il, note)      ! funk
C** new
      dfav = 0.5*(abs(fkn) + abs(fk0))
C                          Runga-Kutta move thru adjustment of saturation
C!!  new self-adjusting time step version of loop
      do while (t .lt. dtu)
        if(this .lt. thr) then
          st = vmin*((thit - thr)/(ths - thr))**eps1(i)
          write(77,277) origin,nloop,  fv,fl,thr,thit,thu,thls,rf/st,
     &    fk1, fk2, fk3, fk4, th2, th3, th4,   dfav  !,  avg0, avg1, avg3
 277  format(1x,a20,i3,7g12.4/7g12.4/7g12.4)
C          this = 0.5*(thit+thr)
          exit
        end if
C
        trace(5:6) = '01'
        thls = this
        avg0 = flams (alm, this, ths, thr, trace) * pb            !flams
        note = 'diag'
        fk1 = funk (fl, rf, this, avg0, fac, thit, il, note)      ! funk
C        if( fk1 .le. 0. .and. rf .le. vmin) then
C          write(99,*) origin, dtu, fl,rf,this,avg0,fac,thit
C          call errxit(notes,' fk1 is negative from funk in rectr')
C        end if
C!! set time step

        IF( fk1 .NE. 0. )THEN ! CU 9/11/2013
          dtt = 0.002/abs(fk1)
          if(dtt .gt. 1.5*odtt) dtt = 1.5*odtt  ! don't grow too fast
        ELSE                  ! CU 9/11/2013
          dtt = 1.5*odtt      ! CU 9/11/2013
        END IF                ! CU 9/11/2013

        if((t+dtt) .gt. dtu) dtt = dtu - t + 1.e-4
C
  50    continue
        odtt = dtt
        th2 = this - 0.5 * dtt * fk1
        th2 = max(th2,thit)

        trace(5:6) = '02'
        avg1 = flams (alm, th2, ths, thr, trace) * pb           !flams

        hft = fl + 0.5 * rf * dtt

        fk2 = funk (hft, rf, th2, avg1, fac, thit, il, note)           ! funk

        th3 = this - 0.5 * dtt * fk2
        th3 = max(th3,thit)

        fk3 = funk ( hft, rf, th3, avg1, fac, thit, il, note)          ! funk

        th4 = this  - dtt * fk3
        th4 = max(th4, thit)

        trace(5:6) = '04'
        avg3 = flams (alm, th4, ths, thr, trace) * pb           ! flams

        ft = fl + rf * dtt

        fk4 = funk (ft, rf, th4, avg3, fac, thit, il, note)            ! funk

        dfac = (fk1 + 2.*(fk2 + fk3) + fk4)/6.
C reset theta and storage depth fl
        this = thls - dtt * dfac
        nloop = nloop + 1
C!! new 4/05
        if(this .lt. thr) then
          if(ncut .lt. 10) then
            dtt = dtt*0.5
            ncut = ncut + 1
            go to 50
          else
            stop ' unable to resolv small th in RECTR with 10 time step 
     &cuts '
          end if
        end if
C!!
        fl = ft
        t = t + dtt
        if(dbg .and. nloop .gt. 100) then         !origin 
          if(rf .gt. vmin .and. this .ge. 0.99*ths) go to 90
          write(99,'(a20,i3,4f9.4/5g13.4/10x,4g12.4)')  origin, nloop,
     &     dtt,t,ft,rf, this, fk1, fk2, fk3, fk4,  th2, th3, th4, dfac 
          if(nloop .gt. 200) stop ' too many loops in rectr'
        end if
C        if(dbg) then
C          write(99,91) nloop, dtt, thls,this, rf, fk1, fk2, avg1, avg3
C  91  format(" in rectr:",i4,9g12.4)
C        end if

      end do
      if(this .lt. thit .or. this .gt. ths) then
        if(dbg) write (99,'(a20," rect out of bounds ",/i4,7g12.4)') 
     &   origin,  nloop,dtu,this,thit,thr,rf,fv,ft
        if(this .lt. thr) stop ' debug rect stop'
      end if

90    continue

      th = this

      if (il .ge. 2) then !th = th1oth2 (this)                        !th1oth2
        th2 = th
        th = th1 + this - thu   ! needed for vol. balance
      end if
      fv = fl

      return

      end
C------------------------------------------------------------------------------

      double precision function flams (a, thw, ths, thr, trace)

      implicit double precision (a-h, o-z)

C   Integral of kr(psi) up to a limit thw, scaled on the value of PB
C   Revised version good up to h = 0.  value of thw first interpreted
C   into equivalent h, and parameterized values of integral used.

C            a   is lambda
C            ths and thr are sat and resid thetas, respectively
C
      character(len=6) :: trace
      logical :: op

      if(thw .gt. ths) thw = ths
      s = (thw - thr) / (ths - thr)
      dp = 1.0
      if(a .lt. 0.2) dp = 0.96 + 0.2*a   ! empirical second order improvement
      GP = GofP(a)

      if (s .gt. 1.e-10) then
C
        hs = hsofs( s, a )                                           ! hsofs
        be = BParm( a )                                              ! BParm
        ae = Aparm( a )                                              ! AParm
        ce = 2.7                                         ! curve coefficient
        if(abs(hs -.05) .le. 1.e-10) then
          flams = GP
        else

          term = (hs-.05)/ae                    ! CU 9/11/2013
          IF( term .GT. 1.E10 )THEN             ! CU 9/11/2013
            flams = GP * term**(-be)            ! CU 9/11/2013
          ELSE                                  ! CU 9/11/2013
            flams = GP*(1.+ term**ce)**(-be/ce) ! CU 9/11/2013
          END IF                                ! CU 9/11/2013

!          flams = GP*(1. +((hs-.05)/ae)**ce)**(-be/ce)       ! GofP; 
        end if
C
      else
C
        iunit = 77   ! output file should always be opened
        inquire (99, opened = op)
        if(op) then  ! only write if diagnostic activated
          iunit = 99
          Write(iunit,"(' Se too small in flams from ',a6,/(5g13.4))") 
     &        trace,  thw, ths, s, a
C        stop ' neg flams '
        end if
C!! new default
        flams = GP*(1. - dp)

      end if

      return

      end

C ----------------------------------------------------------------------

      double precision function hsofs( Se, alam )

      implicit double precision (a-h, o-z)

C  Transitional Brooks-Corey function to obtain abs(scaled h),
C  assuming C = 5 and Sh/Pe = .05:
C
      If(Se .le. 0.  .or. Se .gt. 1.0) then
        stop ' Call to hsofs with Se out of range'
      Else
        C = 5.
        pw = C/alam
        pwi = 1./C
        denom = Se**pw
        if(denom .gt. 0.) then
          term = 1./denom
          hsofs = .05 + (term - 1.0)**pwi
        else
          hsofs = .05 + 1./Se**(1./alam)
        end if
      End If
      Return
      end
C
C ---------------------------------------------------
C
      double precision function GofP( alam )

      implicit double precision (a-h, o-z)
C
      if(alam .gt. 0.) then
        Gofp = (3.4 + 3.*alam)/(2. + 5.*alam) ! Carl has (1.9 + 4.5*alam)
      else
        stop ' call to Gofp with alam out of range'
      end if
      return
      end
C ---------------------------------------------------
C
      double precision function Aparm( alam)
C  fitted curve intercept

      implicit double precision (a-h, o-z)

      if(alam .gt. 0.) then
        Aparm = 0.55 + 0.16* alam
      else
        stop ' call to Aparm with alam out of range'
      end if
      return
      end
C ---------------------------------------------------
C
      double precision function Bparm (alam)
C  fitted curve slope value

      implicit double precision (a-h, o-z)

      if(alam .gt. 0.) then
        Bparm = 1.02 + 3.0*alam
      else
        stop ' call to Bparm with alam out of range'
      end if
      return
      end

C------------------------------------------------------------------------------

      double precision function funk (f, r, thw, gu, fac, thi, il, str)

C   Decay of water content as block of infiltrated water redistributes
C   with surface influx r: d(thet)/dt

      implicit double precision (a-h, o-z)

      common /infil/ twolay(2), zfl(2), f2l, vmu, sk1(2), sk2(2), g1(2),
     &         g2(2), ths1(2), ths2(2), al1(2), al2(2), apu, vmin,
     &      thr1(2), thr2(2), thin2(2), rock1(2), rock2(2), alf, pb1(2),
     &      pb2(2), eps1(2), eps2(2), erfi(2), dtu, trp, i, dbg, intyp

      character(len=4) :: str
      logical twolay, dbg

      integer i, il

      thwt = thw

      if (il .eq. 1) then

        eps = eps1(i)
        ths = ths1(i)
        thr = thr1(i)

      else

        eps = eps2(i)
        ths = ths2(i)
        thr = thr2(i)

      end if

      if (thwt .lt. thr) thwt = thr
      thwu = thwt - thi
C
      thwr = thwt - thr
      if (thwr .lt. 0.) thwr = 0.
C      
      phi = ths - thr
      bf = 0.85
C
      fu = f 
      tm1 = (thwr/phi)**eps
C
      if(fu .gt. 0. .and. thwu .gt. 1.e-4) then
       zuinv = (thwu + .0001)/(fu + .0001)
      else
       zuinv = 5.  ! this is d(theta) divided by infil depth for tiny depths
      endif
C
      tm2 = bf * fac * gu * zuinv
      funk = zuinv * (vmin * (tm1 + tm2) - r)
C
C      if(str .eq. 'diag') then
C        write(99,'(" funk",5g12.4/5g12.4)')funk, thwu, thwt, fu, zuinv, 
C     &   tm1, tm2, gu, vmin, r
C      end if
C 
      return

      end
C------------------------------------------------------------------------------


      double precision function fpond (rr)

      implicit double precision (a-h, o-z)

C   Current ponding F as a function of current rainrate RR, excluding constant
C   total due to KI, if any.  Use 3=parameter infiltrability function.

      common /infil/ twolay(2), zfl(2), f2l, vmu, sk1(2), sk2(2), g1(2),
     &         g2(2), ths1(2), ths2(2), al1(2), al2(2), apu, vmin,
     &      thr1(2), thr2(2), thin2(2), rock1(2), rock2(2), alf, pb1(2),
     &      pb2(2), eps1(2), eps2(2), erfi(2), dtu, trp, i, dbg, intyp

      logical :: twolay, dbg

      integer i


      if (rr .le. vmin) then

        fpond = 1.E6

      else

        fpond = apu / alf * dlog ((rr - vmin
     &          * (1. - alf)) / (rr - vmin))

      end if

      return

      end
C------------------------------------------------------------------------------

      double precision function rkohy (p, al, h, ifd, dkdh)

C   Finds simple relative K as f(H) and dK/dH, using Brooks and Corey.
C   Assumes head and entry head P are both positive.

      implicit double precision (a-h, o-z)

      if (h .lt. p) then

        rkohy = 1.
        if (ifd .gt. 0) dkdh = 0.

      else

        eta = 2. + 3. * al
        rh = h / p
        rkohy = rh**(-eta)
        if (ifd .gt. 0) dkdh = -eta / p * rh**(-eta - 1.)

      end if

      return

      end
C------------------------------------------------------------------------------

      subroutine setg (thi, zb, fv, ht, gef, k, sk, set)

C   Sets appropriate values of G and AP depending on location
C   of wetting front, for 1 or 2 soil layers

      implicit double precision (a-h, o-z)

      common /infil/ twolay(2), zfl(2), f2l, vmu, sk1(2), sk2(2), g1(2),
     &         g2(2), ths1(2), ths2(2), al1(2), al2(2), apu, vmin,
     &      thr1(2), thr2(2), thin2(2), rock1(2), rock2(2), alf, pb1(2),
     &      pb2(2), eps1(2), eps2(2), erfi(2), dtu, trp, i, dbg, intyp

      common /layr/ g2u(2), thm1(2), thm2(2), thi2, rc(2), vmin2(2),
     &              ref, wcf(2), ffl(20,2), sku(20), cv(2), fil, idl, j

      logical :: twolay, op, dbg
      character*6 trace, set*2

      integer i, k, idl

      trace = 'sg   0'
      trace(3:4) = set
      if (.not. twolay(i) .or. zb .lt. zfl(i)) then
C                                                              upper layer only
        sk = sk1(i)

        gef = g1(i) - pb1(i) * flams (al1(i),thi,ths1(i),thr1(i),trace)

        apu = (gef + ht) * (ths1(i) - thi) * (1. - rock1(i))
C 
        if(gef .le. 0.) then
           iunit = 77   ! output file should always be opened
C           inquire (99, opened = op)
C           if(op) iunit = 99
          write(iunit,"(' Gu negative in setg: '/a6,i3,(5g13.4))") 
     &           trace, i, gef, pb1(i),thi, g1(i), ht 
          stop ' negative value of cap. drive in setg for this element '
        end if
      else

        if (sk1(i) .lt. sk2(i)) then
C                                                        upper layer control --
C                                                          matching value for g
          cumi = zfl(i) * (ths1(i) - thi) * (1. - rock1(i))
          defth2 = (thm2(i) - thi2) * (1. - rock2(i))
          trace(5:6) = ' a'
          g2e = g2u(i) - pb2(i)*flams(al2(i),thi2,ths2(i),thr2(i),trace)
C    the following expression corrected 3/2013:
          gfm = alf * cumi / dlog (1. + alf *
     &          vmin2(i) / (rc(i) - vmin2(i))) / defth2
          sk = vmin2(i)
          rx = fv / cumi - 1.
          if (rx .lt. 0.) rx = 0.
          omwt = dexp (-wcf(k) * rx / 16.7)
          gef = omwt * gfm + (1.- omwt) * g2e
C          thwf = thm2(i)
          apu = (gef + ht) * defth2
          if(gef .le. 0.) then
           iunit = 77   ! output file should always be opened
C           inquire (99, opened = op)
C           if(op) iunit = 99
           write(iunit,"(' Gu negative in setg: '/a6,i3,(8g13.4))") 
     &           trace, i, gef, thi, gfm, g2e, omwt, zb, zfl(i)
          stop ' negative value of cap. drive in setg for this element '
          end if

        else
C                                             possible lower layer control
          sk = sk2(i)

          if (g2(i) .ge. 0.0001) then

            trace(5:6) = ' b'
            gef = g2u(i) -pb2(i)*flams(al2(i),thi,ths2(i),thr2(i),trace)

            apu = (gef + ht) * (thm2(i) - thi2) * (1. - rock2(i))  ! ht new, 6/02

            if( apu .le. 0.) then
              iunit = 77   ! output file should always be opened
C              inquire (99, opened = op)
C              if(op) iunit = 99
          write(iunit,"(' APU negative in INFILT/s: '/a6,i3,(5g13.4))") 
     &           trace, i, gef, thi, thi2, thm2(i), gef 
              stop ' negative value of cap. drive in setg for this
     &element '
            end if
          else
            gef = 0.
            apu = 0.
          end if
        end if
      end if

      return

      end
C------------------------------------------------------------------------------


      subroutine errly (hc, ferr, dferr)

C   External function for iteratively finding subcrust head and final f
C   for flow thru Brooks-Corey crust.  Adapted for use by KINEROS2.

      implicit double precision (a-h, o-z)

      common /infil/ twolay(2), zfl(2), f2l, vmu, sk1(2), sk2(2), g1(2),
     &         g2(2), ths1(2), ths2(2), al1(2), al2(2), apu, vmin,
     &      thr1(2), thr2(2), thin2(2), rock1(2), rock2(2), alf, pb1(2),
     &      pb2(2), eps1(2), eps2(2), erfi(2), dtu, trp, i, dbg, intyp

      logical :: twolay, dbg

      integer i
C  values of H and PB are presumed both positive transformed

      hct = hc
      itrr = 0
      dh = abs(hct / 50.)
      vsum = 0.
      ssum = 0.
      hn = hct

      rk2 = rkohy (pb2(i), al2(i), hct, 1, drkdh)

      dkdh = drkdh * sk2(i)
      vkh2 = rk2 * sk2(i)
      vmu = vkh2

      rk1c = rkohy (pb1(i), al1(i), hct, 0, dum)

      vkh1c = rk1c * sk1(i)

10    continue

      ho = hn
      hn = hn - dh

      if (hn .lt. 0.) hn = 0.

      hb = 0.5 * (hn + ho)

      rk1 = rkohy (pb1(i), al1(i), hb, 0, dum)

      vkh1 = sk1(i) * rk1
      term = dh * vkh1 / (vkh2 - vkh1)
      sterm = term / (vkh2 - vkh1)

      if (term .le. 0.) then

        ferr = 10. * hc
        dferr = 50.

        return

      end if

      vsum = vsum + term
      ssum = ssum + sterm
      itrr = itrr + 1

      if (hn .gt. 0.01) go to 10

      ferr = vsum - zfl(i)
      dferr =  vkh1c / (vkh2 - vkh1c) - dkdh * ssum
      vmu = vkh2

      return

      end
C------------------------------------------------------------------------------
C
      double precision function fgex (Ifunc, I, g, a)
C  dimensionless infiltrability as fn. of cum I, 3-par model, no CV.
C  This is the upper envelope of all CV-affected infiltrability curves.
C  a is the soil type parameter alpha, g is the capillary drive

      implicit double precision (a-h, o-z)

      double precision ::   I, Is

      integer:: Ifunc
C
      Is = I/g
      If( Ifunc .ge. 2) then
        vex = a * Is
        if( Is .le. 1.e-8) then
          vpart = 1.e6
        else if(Is .le. 1.e-5) then
          vpart = 1./Is
        else If(vex .gt. 80.) then
          vpart = 0.
        else
          vpart = a / (dexp (a * Is) - 1.)
        end if
      Else   ! G-A
        If (Is .le. 1.e-8 ) then
          vpart = 1.e6
        Else
          vpart = 1.0 / Is
        End If
      End if
      
      fgex = 1. + vpart

      return

      end

C!!-------------------------------------------------------------new 12/02
      double precision function favdI( Infim, Apu, AK, Iv, gam, dt)
C finds increment of potential infil depth in interval dt

      implicit double precision (a-h, o-z)

      double precision :: Iv, Is1, Is2, Isz  !, gam = 0.8

      integer :: Infim
C Iv could be intent IN or returned as new value of cumI      
C
C  get scaling parameters
      Ino = Apu
      Tn = Apu/AK
C  scale Iv
      If( Apu .le. 1.1e-6) then
        write(77,'(3g13.3)') Iv,apu
C        stop ' invalid arg. in favDI '
        FavdI = AK*dt
        return
      end if
      Is1 = Iv/Apu
      if(Is1 .gt. 100.) then
        FavdI = AK*dt
        return
      end if
C
      select case(Infim)
       case( 1 )  ! G-A
        ts1 = Is1 - dlog(Is1 + 1.0)
       case( 2 )  ! S-P
        ts1 = Is1 -1.0 + dexp(-is1)
       case( 3 )
        ts1 = (Is1 -dlog((dexp(gam*Is1) -1. +gam)/gam))/(1. - gam) ! 3-par funct
       end select
C
      Isz = fIsots( ts1, gam)
C  step scaled time
      ts2 = ts1 + dt/Tn
C  approximate decay function
      Is2 = fIsots(ts2, gam)
C  rescale deltaI:
      favdI = (Is2 - Isz)*Apu
C
      return
      end function favdI
C * --------------------------------------------------------------------
      double precision function fIofs(Infim, fs, Apu, gam )
C finds infil depth corresponding to a flux for givenparameters

      implicit double precision (a-h, o-z)

        double precision :: Apu, fs, gam, va, ra

        integer :: Infim
C
        select case (Infim)
        case ( 1 )  ! G-A
         va = 1./(fs - 1.)
        case ( 2 )  ! S-P
         va = dlog( fs/(fs - 1.))
        case ( 3 )  ! 3 - par S-P
         ra = (fs - 1. + gam)/(fs - 1.)
         if(ra .gt. 0.) then
           va = dlog(ra)/gam
         else
           write(99,'(" fs, gam: ",2g12.4)') fs, gam
           stop '  bad value in fsofi '
         end if
        end select
        fIofs = va * Apu
C 
      end function fIofs
C * --------------------------------------------------------------------
      double precision function fIsots(ts, gam)
C
C   explicit CumIs (ts) function of Smith (integral of Eq. 6.45),
C      with omega interpreted from gam. 
C   ts is scaled time, infiltrability from ponding 
C   gam is three-par. weighting function: 1 = S-P, 0 = G-A
C   should work for all gam
C
      implicit double precision (a-h, o-z)

      om = 1./(8.- 6.*gam*gam)
C
      if(ts .le. 0.) then
        fIsots = 0.
      else
        t8 = ts/8.
        tb = 2.*ts
        term1 = ts*(0.75 - om)
        term2 = dsqrt(ts*ts + 8.*ts)/4.
        term3 = dsqrt(t8) + sqrt(t8+1.)
        term4 = dsqrt(om*om*ts*ts + ts/2.)
        term5 = om*dsqrt(tb) + dsqrt(om*om*tb + 1.)
        fIsots = term1 + term2 - 2.*dlog(term3) + term4 + 
     &           0.5*dlog(term5)/om
      end if
      return
      end
C------------------------------------------------------------------------------

      double precision function vcurf (rs, cov)

C   Finds infilt curve coeff, cc, as a function of Cov(Ksat) and R*

      implicit double precision (a-h, o-z)

      if (cov .lt. 0.1) then

        vcurf = 21.

        return

      else if (rs .ge. 0.1) then

        rsf = 1. - dexp (-0.6 * (rs - .1))
        cvf = 0.75 / cov**1.3

      else

        rsf = 0.
        cvf = 1.

      end if

      vcurf = 1.+ rsf * cvf

      if (vcurf .lt. 1.) stop ' error - infilt: curve coeff. < 1'

      return

      end
C------------------------------------------------------------------------------


      double precision function vfm (rs, cov)

C   Finds dimensionless reduction factor to apply to ks to get areal fmin

      implicit double precision (a-h, o-z)

      if (rs .lt. 0.2) then

         vfm = 1.

      else if (cov .lt. 0.1) then
C
          vfm = 1.0

      else

        p = 1.8 / cov**.85
        vfm = (1. + (1. / rs)**p)**(-1. / p)

      end if

      return

      end
C------------------------------------------------------------------------------

      subroutine fcapv (Icode, Ic, fu, ifl, dfdi)

C   Finds infiltrability fu(Ic) for randomly heterogeneous area, and for
C   flag IFL true, finds derivative w/r cumulative infilt Ic

C   Parameters:
C  in common /layr/:
C     ref     rainrate [L/T]
C     vmu     effective final f 
C     cc      curvature coefficient
C
C     Icode: 1=G-A, 2=S-P, 3= 3 param function
C     Ic     cumul. infil depth, [L]
C     fu     infil capacity at infiltrated depth I.....................................................................................................c [L/T]
C     ifl    flag signaling the requirement for a derivative value
C     dfdi   derivative of fv w/r ic

      implicit double precision (a-h, o-z)

      common /infil/ twolay(2), zfl(2), f2l, vmu, sk1(2), sk2(2), g1(2),
     &         g2(2), ths1(2), ths2(2), al1(2), al2(2), apu, vmin,
     &      thr1(2), thr2(2), thin2(2), rock1(2), rock2(2), alf, pb1(2),
     &      pb2(2), eps1(2), eps2(2), erfi(2), dtu, trp, i, dbg, intyp

      common /layr/ g2u(2), thm1(2), thm2(2), thi2, rc(2), vmin2(2),
     &              ref, wcf(2), ffl(20,2), sku(20), cv(2), fil, idl, j

      logical :: twolay, ifl, dbg

      integer :: i, idl, icode

      double precision :: ic, fu

C                                                               define function
C      fexi (al, fu, ag) = exp (al * fu / ag)

      res = ref / vmu

      if (Ic .le. 1.E-6 .and. cv(i) .lt. 0.1) then
C  ! randomness not effective at first:
        fu = vmu * fgex( Icode, Ic, apu, alf)
C        fe = fexi(alf, Ic, apu)
        if(ifl) then
          dfdi = vmu / apu * dfstar(Icode, Ic, apu, alf)
        end if

      else

        if (res .le. 1.) then     ! use .le. rather than .lt.

          fu = ref
          dfdi = 0.

        else

          rat = Ic / apu
          cc = vcurf(res, cv(i))

          if (rat * alf .ge. 38.) then

            fs = 1.
            fe = -10.
            d1 = 0.

          else

C            fe = fexi (alf, Ic, apu)

            rsm = res - 1.
C
            core = fcore ( icode, IC, apu, alf )
C            if (ifl)  d1 = vmu * rsm**2 / apu
            if( ifl ) d1 = vmu / apu * dfstar( icode, IC, apu, alf )
C
C            core = (fe - 1.)/ alf

            if (core .le. 0.) then      ! new test

              fs = 1. + rsm

            else

              vu = rsm * core 
              test = cc * dlog (vu)

              if (test .lt. 39. .and. cc .lt. 20.) then

                if (test .gt. -40) then

                  fs = 1. + rsm * (1.+ vu**cc)**(-1. / cc)

                  if (ifl) d1 = d1 * vu**(cc - 1.) /
     &                          (1. + vu**cc)**(1. + 1. / cc)

                else

                  fs = 1. + rsm
                  d1 = 0.

                end if

              else  !  insignificant CV:

C                fs = 1. + alf / (fe - 1.)
                fs = fstar( Icode, Ic, apu, alf)
                if (ifl) d1 = vmu/apu * dfstar( Icode, Ic, apu, alf )

              end if
            end if
          end if

          fu = fs * vmu
          if(ifl) dfdi = d1 

        end if
      end if

      return

      end subroutine fcapv
C---------------------------------------------
      double precision function fexi (al, fu, ag) 

      implicit double precision (a-h, o-z)

        double precision :: al, fu, ag
C
        fexi  = dexp (al * fu / ag)

        if (fexi .eq. 1.) fexi = 1.00001
      
      end function fexi
C---------------------------------------------
      double precision function fstar (Infc, Ic, apu, alf )

      implicit double precision (a-h, o-z)

        integer :: Infc

        double precision    :: apu, alf, fe, Ic, Istar
C
        select case (Infc)
C
          case(1)
            Istar = Ic / apu
            fstar = (1. + Istar ) / Istar
          case( 2,3)
            fe = fexi( alf, Ic, apu)
            fstar = alf / (fe - 1.) + 1.
        end select
C      
      end function fstar
C----------------------------------------------
      double precision function dfstar ( Infc, Ic, Apu, alf )

      implicit double precision (a-h, o-z)

        integer :: Infc

        double precision    :: apu, alf, fe, Ic
C
        select case (Infc)
          case(2,3)
            fe = fexi(alf, Ic, apu)  ! exp function
            dfstar = - alf*alf*fe/( fe-1.) ** 2
          case(1)
            dfstar = -(Ic / apu) ** 2
        end select
      end function dfstar
C
C----------------------------------------------
C
      double precision function fcore (Infc, Ic, Apu, alf )

      implicit double precision (a-h, o-z)

        integer :: Infc

        double precision    :: apu, alf, fe, Ic
C
        select case (Infc)
          case(2,3)
            fe = fexi(alf, Ic, apu)  ! exp function
            fcore = (fe - 1.) / alf
          case(1)
            fcore = 1.0
        end select
C
      end function fcore
C
C
