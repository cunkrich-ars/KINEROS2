C   modules for K2
      
      module multip
        real :: rmks, rmn, rmcv, rmg, rmin, rmcoh, rmspl, rmwco, rmlen
        real :: rmn_chan, rmsat
      end module multip
C
      module runpars           ! generally global variables
        logical :: silent

        real :: imperv ! total impervious area
        real :: bed    ! total channel bed area
        real :: bedvol ! total channel sand volume

C TEMP SNOW
        CHARACTER(LEN=5) :: mmdd_str
        CHARACTER(LEN=4) :: yr_str
        REAL :: ice
        LOGICAL :: snow
C TEMP SNOW

        integer elapsed_time

        type sumval
          integer :: idel, itype 
C  idel   is the identifying id number
C  itype  is the saved value of ltyp
          logical :: twola
C  twola  is true when the element has a two layer soil (planes, urban, main chan)
          real :: are, abot, cumare, volin, volrn, volro, ropeak,
     &            ftot, flowr, flcap, vbot, pored, thst, sedout 
C                this is a structure to hold the summary data for each element:
C  are    is the element area in m**2 or ft**2 which receives rainfall
C  abot   is the channel bottom area
C  cumare is the cumulative area at the outlet of the element
C  volin  is the volume of water entering the element exclusive of rain; m3 or ft3
C  volrn  is the volume of rain falling on the element in m**3 or ft**3
C  volro  is the volume of runoff, m3 or ft3
C  ropeak is the peak runoff rate in mm/h or in/hr
C  ftot   is the total infiltration in mm or in.
C  flowr  is the infiltration into the lower layer, if any, in mm or in.
C  flcap  is the field cap storage in the upper layer, if twola
C  vbot   is the channel infil into the channel bottom, m^3
C  pored  is the init. pore storage of the upper layer, if 2 layers. mm or in.
C  thst   is the starting surface water content 
C  sedout is the sediment output for each element in the output table
C
        end type sumval
C
        integer, parameter  :: files(4) = (/55, 66, 77, 88/) ! parameter, rain, output, print=3
        logical :: sed, cour, chrain, tabl, vapi, pond_rain, gage_rain
        logical :: sed_debug
!   sed    .true. = route sediment,
!   cour   .true.   adjust time step for Courant condition,
!          .false.  do not adjust (but report),
!  chrain  .true. = include rain on all channels
!   tabl   .true. = print table of hydrologic statistics for each element
C   vapi   .true. = use API value and setup initial surface pulse using SAT
        integer :: nele, limit, units, nps, nelt, ltab, p1, p2, p3
!   nele      total number of elements 
!   limit     number of time steps,
!   p1, etc.  quarter points of limit
!   units     1 = metric, 2 = english
!   nps       number of particle size classes:  if SED is simulated
!   nelt      cumulative count of number of elements simulated
        real, parameter :: dxlmt = 30.
C** new dx max here for reality check
        real :: clen, tfin, delt, conv, wt, arconv, wconv, xnu, depmax,
     &          thfc, thwilt, thinu, satfc, satwl, satin, thru1,
     &          thru2, thinu2, thsv1, thsv2, pbv1, pbv2, vks1, vks2
!   clen     characteristic length (m or ft),
!   tfin     elapsed time of event simulation, min
!   delt     user-specified time step (sec)(input in min.),
C      delt changed to seconds, tfin in min.
C Units conversions, depending on units of run:
!   conv     conversion for length from in to ft or mm to m
!     wt     water density in kg/m^3 or lb./ft^3
! arconv     area conversion:  m^2/ha or ft^2/ac
!  wconv     mass conversion:  kg/TN or lb./ton
!  xnu       kinematic viscosity m^2/s or ft^2/s
! depmax     max. upper layer depth over watershed, m
!  thfc      upper soil water content at -0.3 bars pressure head
!  thwilt    upper soil water content at -15 bars
!  thinu     upper soil water content used at start of runoff-producing rain
 ! satfc     rel. (scaled) value of thfc
 !  (etc.)
        real, dimension(4) :: dtm
C dtm(4)     smallest dt to satisfy the Courant condition
C            during each quarter of the simulation
        real, dimension(5) :: rho, dpd, vsl, pwt
        integer, dimension(5) :: nord
        real :: grav, temp, bdep   
!  rho   is solid density of particles by class
!  dpd   is particle class diameter in mm or in; used internally in m or ft.
!  vsl   is settling velocity
!  nord  is order of particle classes by erodability
        character(LEN=70) :: title, fname(3), filrun
C  title    is 60 character run identifying title, in run file or from 
C           options dialog.
C  fname(1) is input parameter file name
C       (2) is input rain file
C       (3) is output
C  filrun is the runfile name
        character(LEN=16),dimension(0:8) :: typname = (/
     &        'Plane Element   ','Channel Elem.   ','Pipe Element    ',
     &        'Pond Element    ','Compound Chan   ','Inject Elem.    ',
     &        'Urban Element   ','Overbank Area   ','Adder Element   '/)
C  typname is Name associated with ltyp, as defined below
        character(LEN=2) :: dlab, dilb
 !  dlab   abbrev. major length units [m or ft]
 !  dilb           minor length units [ mm or in]
        type (sumval), dimension(1000000) :: sumtab
 !
!  summary statistics kept for tabular reporting
        integer, parameter :: xfile = 1955

!TEMP 21AUG2013
!        INTEGER :: IDs(1000)
!        REAL    :: Hav(1000,1000)
!        INTEGER :: iplane
!TEMP 21AUG2013
      end module runpars
C
      module elpars       ! generally element particular parameters
        logical ::  diag, rsoil 
!   diag     .true. = write diagnostic info to unit 99,
!   rsoil           = real soil (not fixed f or impervious)
        integer :: nprint, ltyp, id, nchan , ndv, infm, iexceed
C   nprint is code for printout option: 0 = none, 1 = minimum, 2 = hydrograph,
C                                       3 = spreadsheet hydrograph output(unit 119)
C                                           to FILE name given in parameter file
C   ltyp is element type:  0 = plane
C                          1 = simple channel
C                          2 = pipe
C                          3 = pond
C                          4 = compound channel
C                          5 = inject element
C                          6 = urban element
C                          7 = overbank area
C                          8 = adder element
C   nchs  is storage dimension, 2 for channels which may have nchan=2
C   nchan is 1 in general, only 2 for operations on the overbank of compound channel
!   ndv   is no. of surface divisions, used in writer
!   iexceed is time step pond rating was exceeded
        real :: barea, avai, apa, alfinf !dx, xlen
C   barea is channel bottom area
C   avai, for channels, is the weighted mean infiltrating area
C   apa is, for api option, the actual depth of api pulse used
C
        real,dimension(15) :: qbal
C    qbal(1)  = area                  (sq.m/sq.ft),
C    qbal(2)  = rainfall              (cu.m/cu.ft),
C    qbal(3)  = baseflow              (cu.m/cu.ft),
C    qbal(4)  = injection             (cu.m/cu.ft),
C    qbal(5)  = initial pond storage  (cu.m/cu.ft),
C    qbal(6)  = plane infiltration    (cu.m/cu.ft), (or overbank)
C    qbal(7)  = channel infiltration  (cu.m/cu.ft),
C    qbal(8)  = interception          (cu.m/cu.ft),
C    qbal(9)  = storage               (cu.m/cu.ft),
C    qbal(10) = outflow               (cu.m/cu.ft),
C    qbal(11) = chan. bottom infil    (cu.m/cu.ft),
C    
        character(LEN=20):: sumro, sumsed    ! set in writer & forwarded to graf

!K2_PROFILER
        logical :: profiler
        real :: conc2(5,20), dnet1(5,20), dnet2(5,20)
        real :: eros1(5,20), eros2(5,20)
!K2_PROFILER

        integer                 :: xpond_n, xchan_n
        real,dimension(1000)    :: xpond_pct, xchan_pct
        integer,dimension(1000) :: xpond_id, xchan_id

      end module elpars
C
      module itpars
        integer :: itm, itlim
C  iteration status parameters
      end module itpars
C
      module grafdat
        integer ::  ndat 
        logical :: ifad = .false.    ! flag indicating input of actual data
        character(len=12) :: qlabl=' Flow, mm/hr', tlabl=' time, min. '
        character(LEN=9),dimension(2) :: slabl = (/'Sed. Flux',
     &           '  kg/s   '/)
        real, dimension(500) :: tdat, qdat, sdat   !  actual plot data, if any
        real, dimension(:),allocatable :: tgr, rgr, qd, sd
      end module grafdat
C  the following are modules used in the windows grafiks shell.
      module messagew
        character (len = 70), dimension (3) :: msgstring
        CHARACTER (LEN =210) :: OUTSTR
        Character (LEN = 20) :: titlbar = 'Warning Message     '
        character (len = 70) :: blankline = 
     &'                                                                 
     &     '
        equivalence(outstr,msgstring)
        real     :: xvalu
        integer  :: ibutt, iconx, nlm
        logical  :: errxt
      end module messagew
C 
      module putgets
        character(len=200)  :: stringtx
        character(len=50), dimension (4) :: partstrng
        character(len=12)   :: charvalu
        integer             :: ihand1
        logical             :: errout
        equivalence (partstrng,stringtx)
      end module putgets
C
