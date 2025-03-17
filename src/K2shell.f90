! Run K2 from a command shell
!
! Usage:
!
! k2 -b file
!
!   k2 is the K2shell executable file
!   -b specifies batch mode
!   -s specifies batch mode and suppresses screen output
!   file is the run file (default = kin.fil)
!
! Batch mode
!
! The run file is a series of lines, one line for each simulation, with the following comma-separated values:
!
! parameter-file, rainfall-file, output-file, title, tfin, delt, cour?, sed?, multiplier-file?, table?
!
! The title must be enclosed in double quotes.
!
! Entries after delt are optional, will default to 'N'

! Optional file: mlist.fil - list subset of elements on which to apply multipliers (1 column)

! Citation:
! Goodrich, D.C., I.S. Burns, C.L. Unkrich, D.J. Semmens, D.P. Guertin, M. Hernandez, S. Yatheendradas, J.R. Kennedy,
! and L.R. Levick. 2012. K2/AGWA: model use, calibration, and validation. Transactions of the ASABE 55(4):1561-1574.
! 
! Disclaimer:
! No warranties, expressed or implied, are made that this program is free from errors or will meet the requirements
! for any particular application.  The U.S. Department of Agriculture disclaims all liability for direct or
! consequential damages resulting from the use of this program.

program K2shell

  use runpars
  use elpars
  use multip
  use kinsed_def


  integer, parameter  :: PROMPT = 1
  integer, parameter  :: READ_FROM_FILE = 2
  integer, parameter  :: NO_MULTIPLIERS = 3
  integer, parameter  :: INTERACTIVE = 4
  integer, parameter  :: BATCH = 5

  real, parameter     :: RMULT = 1.0

  integer             :: arg_len ,rfile, file_error, action, ierr, ngage, nd, mlist(100000), nlist
  integer             :: tvals(8), t0, te
  logical             :: file_exists, char2log, call_writer, unif_sat, cmults, extra_mults
  real                :: t(5000), z(5000), sat, sumbal(10)
  real                :: rmks_, rmn_, rmcv_, rmg_, rmin_, rmcoh_, rmspl_
  real                :: rmks_alt, rmn_alt, rmcv_alt, rmg_alt, rmin_alt, rmcoh_alt, rmspl_alt
  double precision    :: x, y
  character(len=1)    :: answer, cour_char, sed_char, tabl_char, c, bell
  character(len=2)    :: flag
  character(len=10)   :: block_name
  character(len=20)   :: tfin_char, delt_char, rmks_char, rmn_char, rmcv_char
  character(len=20)   :: rmg_char, rmin_char, rmcoh_char, rmspl_char
  character(len=20)   :: rmks_chan_char, rmg_chan_char, rmn_chan_char, rmwco_char, rmlen_char, rmsat_char
  character(len=11)   :: version
  character(len=150)  :: parfile, rainfile, outfile, multfile, auxfile, runfile, fields(20)
  character(len=1000) :: line

real :: v_urban, v_plane, f_plane, f_pond, f_urban, f_chan, s_pond

! TEMP SNOW
  INTEGER, PARAMETER :: DAYS(12) = (/0,31,60,91,121,152,182,213,244,274,305,335/)
!  snow = .TRUE.
  snow = .FALSE.
!  WRITE( 99, '("year, month, day, sat, rainfall volume (mm), runoff volume (mm), sediment yield (kg/ha)")' )
! TEMP SNOW

bisects = 0.
ncalls_bisects = 0.
ncalls  = 0.

v_urban = 0.
v_plane = 0.
f_plane = 0.
f_pond  = 0.
f_urban = 0.
f_chan  = 0.
s_pond  = 0.


10 format(a) ! for general character i/o

  version = '25-Oct-2019'
  bell    = char(7)
  silent  = .FALSE.

  extra_mults = .FALSE.

! Default values

  cour = .false.
  sed  = .false.
  tabl = .false.

  rmks_chan = 1.0
  rmg_chan  = 1.0
  rmn_chan  = 1.0
  rmwco     = 1.0
  rmlen     = 1.0
  rmsat     = 1.0
     
! Command line & run file

  if(COMMAND_ARGUMENT_COUNT() == 0) then
    mode = INTERACTIVE
    inquire(file = 'kin.fil', exist = file_exists)
    if(file_exists) then
      write(*, "('$ Repeat previous run? ')")
      read(*, 10) answer
      if(char2log(answer)) then
        open(4, file = 'kin.fil', status = 'old', iostat = file_error)
        if(file_error .ne. 0) call errxit('Kineros2', "Can't open kin.fil")
        action = READ_FROM_FILE
      else
        action = PROMPT
      end if
    else
      action = PROMPT
    end if
  else
    call GET_COMMAND_ARGUMENT( 1, flag, arg_len )
    if(flag == '-b' .or. flag == '-s') then
      mode = BATCH
      if(COMMAND_ARGUMENT_COUNT() == 2) then
        call GET_COMMAND_ARGUMENT( 2, runfile, arg_len )
      else
        runfile = 'kin.fil'
        arg_len = 7
      end if
      open(4, file = runfile(1:arg_len), status = 'old', iostat = file_error)
      if(file_error .ne. 0) call errxit('Kineros2', "Can't open "//runfile(1:arg_len))
      if(flag == '-s') silent = .TRUE.
    end if
  end if

! File with list of elements on which to apply multipliers (optional)

  nlist = 0
  OPEN( 1, FILE = 'mlist.fil', STATUS = 'old', IOSTAT = file_error )
  IF( file_error == 0 )THEN
    READ( 1, * ) rmks_alt, rmn_alt, rmcv_alt, rmg_alt, rmin_alt, rmcoh_alt, rmspl_alt
    DO
      READ( 1, *, IOSTAT = file_error ) i
      IF( file_error == 0 )THEN
        nlist = nlist + 1
        mlist(nlist) = i
      ELSE
        EXIT
      END IF
    END DO
  END IF
  CLOSE( 1 )

! File with channel-only multipliers (optional)

!  cmults = .FALSE.
!  rmwco  = 1.
!  rmn_chan = 1.
!  OPEN( 1, FILE = 'cmult.fil', STATUS = 'old', IOSTAT = file_error )
!  IF( file_error == 0 )THEN
!    cmults = .TRUE.
!    READ( 1, * ) rmks_chan
!    READ( 1, * ) rmg_chan
!    READ( 1, * ) rmwco
!  END IF
!  CLOSE( 1 )

! Begin simulation loop

  do

    CALL DATE_AND_TIME( values = tvals )

    t0 = tvals(5) * 3600 + tvals(6) * 60 + tvals(7)

    if(mode == INTERACTIVE) then
    
      write(*, 10)
      write(*, 10) '  KINEROS2'
      write(*, 10) '  Kinematic Runoff and Erosion Model'
      write(*, 10) '  Version '//version
      write(*, 10) '  U. S. Department of Agriculture'
      write(*, 10) '  Agricultural Research Service'
      write(*, 10)
      write(*, 10) '  Citation: Goodrich, D.C., I.S. Burns, C.L. Unkrich, D.J. Semmens,'
      write(*, 10) '  D.P. Guertin, M. Hernandez, S. Yatheendradas, J.R. Kennedy, and'
      write(*, 10) '  L.R. Levick. 2012. K2/AGWA: model use, calibration, and validation.'
      write(*, 10) '  Transactions of the ASABE 55(4):1561-1574.'
      write(*, 10)
      write(*, 10) '  Disclaimer: No warranties, expressed or implied, are made that this'
      write(*, 10) '  program is free from errors or will meet the requirements for any'
      write(*, 10) '  particular application.  The U.S. Department of Agriculture'
      write(*, 10) '  disclaims all liability for direct or consequential damages'
      write(*, 10) '  resulting from the use of this program.'
      write(*, 10)

      select case(action)
      case(READ_FROM_FILE)
        read(4, 10) parfile
        read(4, 10) rainfile
        read(4, 10) outfile
        read(4, 10) title
        read(4, 10) tfin_char
        read(4, 10) delt_char
        read(4, 10) cour_char
        read(4, 10) sed_char
        read(4, 10) multfile
        read(4, 10, iostat = file_error) tabl_char
        if(file_error .ne. 0) tabl_char = 'n'
        close(4)
      case(PROMPT)
        write(*, "('$            Parameter file: ')")
        read(*, 10) parfile
        write(*, "('$             Rainfall file: ')")
        read(*, 10) rainfile
        write(*, "('$               Output file: ')")
        read(*, 10) outfile
        write(*, "('$               Description: ')")
        read(*, 10) title
        write(*, "('$            Duration (min): ')")
        read(*, 10) tfin_char
        write(*, "('$           Time step (min): ')")
        read(*, 10) delt_char
        write(*, "('$ Courant Adjustment? (y/n): ')")
        read(*, 10) cour_char
        write(*, "('$           Sediment? (y/n): ')")
        read(*, 10) sed_char
        write(*, "('$   Multipliers? (y/n/file): ')")
        read(*, 10) multfile
        write(*, "('$    Tabular Summary? (y/n): ')")
        read(*, 10) tabl_char
      end select
    
!     Open parameter, rainfall and output files - the rainfall file is treated
!     as optional (the file unit is set to -1 as a flag for subroutine rain to
!     create a single raingage entry with zero depths)
     
      open(files(1), file = parfile, status = 'old', iostat = file_error)
      if(file_error .ne. 0) call errxit('K2shell', "Can't open "//parfile)
     
      if(len_trim(rainfile) .eq. 0) then ! no rainfall file specified     
        rfile = -1
        rainfile = 'None'
      else
        open(files(2), file = rainfile, status = 'old', iostat = file_error)
        if(file_error .ne. 0) call errxit('K2shell', "Can't open "//rainfile)
        rfile = files(2)
      end if
     
      open(files(3), file = outfile, status = 'unknown', iostat = file_error)
      if(file_error .ne. 0) call errxit('K2shell', "Can't open "//outfile)
     
!     Decode character input into floating and logical values
     
      read(tfin_char , *) tfin
      read(delt_char , *) delt
      cour = char2log(cour_char)
      sed = char2log(sed_char)
      tabl = char2log(tabl_char)
    
!     Three possibilities for parameter multipliers:
    
!     1) Prompt
!     2) Read from file specified by user or from mult.fil
!     3) No multipliers 
    
      if(len_trim(multfile) .gt. 1) then ! user specified a file
        open(4, file = multfile, status = 'old', iostat = file_error)
        if(file_error .ne. 0) call errxit('K2shell', "Can't open "//multfile)
        action = READ_FROM_FILE
      else
        if(char2log(multfile)) then ! user wants multipliers
          inquire(file = 'mult.fil', exist = file_exists)
          if(file_exists) then
            write(*, "('$ Use previous multipliers? ')")
            read(*, 10) answer
            if(char2log(answer)) then ! read multipliers from mult.fil
              open(4, file = 'mult.fil', status = 'old', iostat = file_error)
              if(file_error .ne. 0) call errxit('K2shell', "Can't open mult.fil")
              action = READ_FROM_FILE
            else
              action = PROMPT
            end if
          else
            action = PROMPT
          end if
        else ! user doesn't want multipliers
          action = NO_MULTIPLIERS
        end if
      end if
    
      select case(action)
      case(READ_FROM_FILE)
        read(4, 10) rmks_char
        read(4, 10) rmn_char
        read(4, 10) rmcv_char
        read(4, 10) rmg_char
        read(4, 10) rmin_char
        read(4, 10) rmcoh_char
        read(4, 10) rmspl_char
        close(4)
      case(PROMPT)
        write(*, "('$            Ks multiplier (0-1): ')")
        read(*, 10) rmks_char
        write(*, "('$ Manning/Chezy multiplier (0-1): ')")
        read(*, 10) rmn_char
        write(*, "('$            CV multiplier (0-1): ')")
        read(*, 10) rmcv_char
        write(*, "('$             G multiplier (0-1): ')")
        read(*, 10) rmg_char
        write(*, "('$  Interception multiplier (0-1): ')")
        read(*, 10) rmin_char
        write(*, "('$ Soil cohesion multiplier (0-1): ')")
        read(*, 10) rmcoh_char
        write(*, "('$   Rain splash multiplier (0-1): ')")
        read(*, 10) rmspl_char
      case(NO_MULTIPLIERS)
        rmks_char  = '1.0'
        rmn_char   = '1.0'
        rmcv_char  = '1.0'
        rmg_char   = '1.0'
        rmin_char  = '1.0'
        rmcoh_char = '1.0'
        rmspl_char = '1.0'
      end select
     
!     Decode multipliers into floating values (module multip)
     
      read(rmks_char,  *) rmks_
      read(rmn_char,   *) rmn_
      read(rmcv_char,  *) rmcv_
      read(rmg_char,   *) rmg_
      read(rmin_char,  *) rmin_
      read(rmcoh_char, *) rmcoh_
      read(rmspl_char, *) rmspl_
     
!     Save event info to kin.fil
     
      lpar  = len_trim(parfile)
      lrain = len_trim(rainfile)
      lout  = len_trim(outfile)
      lmult = len_trim(multfile)
     
      open(4, file = 'kin.fil', status = 'unknown')
     
      write(4, 10) parfile(1:lpar)
      write(4, 10) rainfile(1:lrain)
      write(4, 10) outfile(1:lout)
      write(4, 10) title
      write(4, 10) tfin_char
      write(4, 10) delt_char
      write(4, 10) cour_char
      write(4, 10) sed_char
      write(4, 10) multfile(1:lmult)
      write(4, 10) tabl_char
      close(4)
     
!     Save multipliers in mult.fil
     
      open(4, file = 'mult.fil', status = 'unknown')
     
      write(4, 10) rmks_char
      write(4, 10) rmn_char
      write(4, 10) rmcv_char
      write(4, 10) rmg_char
      write(4, 10) rmin_char
      write(4, 10) rmcoh_char
      write(4, 10) rmspl_char
      close(4)
     
!     Also write info to output file
     
      write(files(3), 20) version,           &
                          title,             &
                          parfile(1:lpar),   &
                          rainfile(1:lrain), &
                          tfin_char,         &
                          delt_char,         &
                          cour_char,         &
                          sed_char,          &
                          multfile(1:lmult), &
                          tabl_char
     
20    format(' KINEROS2 ',a//                   &
             ' Title: ',a//                     &
             ' Parameter File Used....... ',a/  &
             ' Rainfall File Used........ ',a// &
             ' Length of Run, minutes.... ',a/  &
             ' Time Step, minutes........ ',a/  &
             ' Use Courant criteria?..... ',a/  &
             ' Simulate Sed. Transport?.. ',a/  &
             ' Multiplier file (if any).. ',a/  &
             ' Tabular Summary?.......... ',a)
      
      write(files(3), 30) rmks_char, rmn_char, rmcv_char, rmg_char, &
                          rmin_char, rmcoh_char, rmspl_char
    
30    format(/' Multipliers Used:'//            &
              ' Saturated Conductivity... ',a/  &
              ' Manning n................ ',a/  &
              ' CV of Ksat............... ',a/  &
              ' Capillary Drive Coeff.... ',a/  &
              ' Intercepted Depth........ ',a/  &
              ' Sediment Cohesion Coeff.. ',a/  &
              ' Sediment Splash Coeff.... ',a/)
    
    else ! BATCH mode
    
      read(4, 10, iostat = file_error) line
   
      if(file_error .ne. 0) exit
   
!     Replace double quotes around title field with spaces
   
      i = index(line, '"')
      line(i:i) = ' '
      j = index(line, '"')
      line(j:j) = ' '
   
!     Temporarily replace any commas in title field with bell characters
     
      k = index(line(i:j), ',')
      do while( k > 0)
        k = i + k - 1
        line(k:k) = bell
        k = index(line(i:j), ',')
      end do
   
!     Identify comma-separated fields
     
      fields = ' '
      i = 1
      j = 1
      k = index(line, ',')
      do while( k > 0)
        line(k:k) = ' '
        fields(i) = line(j:k-1)
        i = i + 1
        j = k + 1
        k = index(line, ',')
      end do
      fields(i) = line(j:)
   
!     Restore commas to title field
     
      k = index(fields(4), bell)
      do while( k > 0)
        fields(4)(k:k) = ','
        k = index(fields(4), bell)
      end do
   
!     Open parameter, rainfall and output files - the rainfall file is treated
!     as optional (the file unit is set to -1 as a flag for subroutine rain to
!     create a single raingage entry with zero depths)
   
      open(files(1), file = trim(fields(1)), status = 'old', iostat = file_error)
   
      if(file_error .ne. 0) call errxit('K2shell', "Can't open "//trim(fields(1)))
   
      if(len_trim(fields(2)) == 0) then ! no rainfall file specified     
        rfile  = -1
        fields(2) = 'None'
      else
        open(files(2), file = fields(2), status = 'old', iostat = file_error)
        if(file_error .ne. 0) call errxit('K2shell', "Can't open "//fields(2))
        rfile  = files(2)
      end if

! TEMP SNOW
      IF( snow )THEN
        yr_str  = fields(2)(8:11)
        READ( fields(2)(13:15), * ) k ! doy
        DO j = 12, 1, -1
          IF( DAYS(j) < k ) EXIT
        END DO
        i = k - DAYS(j) ! day of month
        WRITE( mmdd_str, '(I2.2,",",I2.2)' ) j, i
      END IF
! TEMP SNOW

      open(files(3), file = fields(3), status = 'unknown', iostat = file_error)
      if(file_error .ne. 0) call errxit('K2shell', "Can't open "//fields(3))
   
!     Decode character input into floating and logical values (module runpars)
   
      read(fields(5) , *) tfin
      read(fields(6) , *) delt
     
!     If no value for a logical variable is specified, it remains at its default value
     
      if(len_trim(fields(7))  > 0) cour = char2log(fields(7))
      if(len_trim(fields(8))  > 0) sed  = char2log(fields(8))
      if(len_trim(fields(10)) > 0) tabl = char2log(fields(10))
   
!     Two options for parameter multipliers:
   
!     1) Use file specified in kin.fil
!     2) No multipliers 
   
      if(len_trim(fields(9))      == 0 .or. &
           verify(fields(9),'N ') == 0 .or. &
           verify(fields(9),'n ') == 0) then ! no multipliers
   
        rmks_char  = '1.0'
        rmn_char   = '1.0'
        rmcv_char  = '1.0'
        rmg_char   = '1.0'
        rmin_char  = '1.0'
        rmcoh_char = '1.0'
        rmspl_char = '1.0'
    
        rmks_chan_char = '1.0'
        rmg_chan_char  = '1.0'
        rmn_chan_char  = '1.0'
        rmwco_char     = '1.0'
        rmlen_char     = '1.0'
        rmsat_char     = '1.0'

        fields(9) = 'N'
    
      else
    
        open(7, file = fields(9), status = 'old', iostat = file_error)
        if(file_error .ne. 0) call errxit('K2shell', "Can't open "//fields(9))
    
        read(7, 10) rmks_char
        read(7, 10) rmn_char
        read(7, 10) rmcv_char
        read(7, 10) rmg_char
        read(7, 10) rmin_char
        read(7, 10) rmcoh_char
        read(7, 10) rmspl_char

        read(7, 10, iostat = file_error) rmks_chan_char
        if(file_error .ne. 0)THEN
          rmks_chan_char = '1.0'
          extra_mults = .FALSE.
        else
!          cmults = .TRUE.
          extra_mults = .TRUE.
        end if

        read(7, 10, iostat = file_error) rmg_chan_char
        if(file_error .ne. 0) rmg_chan_char = '1.0'

        read(7, 10, iostat = file_error) rmn_chan_char
        if(file_error .ne. 0) rmn_chan_char = '1.0'

        read(7, 10, iostat = file_error) rmwco_char
        if(file_error .ne. 0) rmwco_char = '1.0'

        read(7, 10, iostat = file_error) rmlen_char
        if(file_error .ne. 0) rmlen_char = '1.0'

        read(7, 10, iostat = file_error) rmsat_char
        if(file_error .ne. 0) rmsat_char = '1.0'

        close(7)
    
      end if
    
!     Decode multipliers into floating values (module multip)
     
      read(rmks_char,  *) rmks_
      read(rmn_char,   *) rmn_
      read(rmcv_char,  *) rmcv_
      read(rmg_char,   *) rmg_
      read(rmin_char,  *) rmin_
      read(rmcoh_char, *) rmcoh_
      read(rmspl_char, *) rmspl_

      read(rmks_chan_char, *) rmks_chan
      read(rmg_chan_char,  *) rmg_chan
      read(rmn_chan_char,  *) rmn_chan
      read(rmwco_char,     *) rmwco
      read(rmlen_char,     *) rmlen
      read(rmsat_char,     *) rmsat
     
!     Write info to output file
     
      write(files(3), 20) version, trim(fields(4)), trim(fields(1)), &
                          trim(fields(2)), (trim(fields(i)), i = 5, 10)
   
      write(files(3), 30) rmks_char, rmn_char, rmcv_char, rmg_char, &
                          rmin_char, rmcoh_char, rmspl_char


      if(extra_mults) then

        write(files(3), 40) rmks_chan_char, rmg_chan_char, rmn_chan_char, &
                            rmwco_char, rmlen_char, rmsat_char

40      format(/' Additional Multipliers:'//            &
                ' Channel Ks............... ',a/  &
                ' Channel G................ ',a/  &
                ' Channel n................ ',a/  &
                ' Woolhiser Coeff.......... ',a/  &
                ' Channel Length........... ',a/  &
                ' Initial Saturation....... ',a/)
      end if

    end if

!   Initialize global volumes
    
    do j = 1, 10
      sumbal(j) = 0.
    end do
    
!   Compute number of time steps
    
    limit = int (tfin / delt) + 1
    
!   Convert time increment to seconds
    
    delt = delt * 60.
    
!   Read global parameter block
    
    call reader(files(1), block_name, ierr)
    if(ierr .gt. 0) call errxit('K2shell', "Invalid global block")
    
    call getr4('C', 0, clen, ierr)
    if(ierr .ne. 0) call errxit('K2shell', "char. length not found")
    
!   System of units
    
    call getstr ('U', 0, c, l, ierr)
    if(ierr .ne. 0) call errxit('K2shell', "units not specified")
    
    if(c .eq. 'm' .or. c .eq. 'M') then ! metric
      units  = 1
      dlab   = 'm.'
      dilb   = 'mm'
      conv   = 1000.
      wt     = 1000.  ! kg/m3 water density
      arconv = 10000. ! m^2 per ha.
      wconv  = 1000.  ! kg/tonne
      grav   = 9.81
      bdep   = 656.
    else if(c .eq. 'e' .or. c .eq. 'E') then ! english
      units  = 2
      dlab   = 'ft'
      dilb   = 'in'
      conv   = 12.
      wt     = 62.4    
      arconv = 43560.
      wconv  = 2000.
      grav   = 32.2
      bdep   = 200.
    else
      call errxit('K2shell', "units not recognized")
    end if
    
    write (*, '(/)')
    
!   Diagnostic output to unit 99?
    
    call getstr( 'DIAG', 0, c, l, ierr)
    if(ierr .eq. 0 .and. (c .eq. 'y' .or. c .eq. 'Y')) then
      diag = .true.
      open(99, file = 'diagnostics.txt', status = 'unknown')
    else
      diag = .false. ! default
    end if
    
!   Rain on channels?
    
    call getstr( 'CHR', 0, c, l, ierr)
    if(ierr .eq. 0 .and. (c .eq. 'y' .or. c .eq. 'Y')) then
      chrain = .true.
    else
      chrain = .false. ! default
    end if
    
!   Rain on ponds?
    
    call getstr( 'PR', 0, c, l, ierr)
    if(ierr .eq. 0 .and. (c .eq. 'y' .or. c .eq. 'Y')) then
      pond_rain = .true.
    else
      pond_rain = .false. ! default
    end if

!   Assign rain gages sequentially to elements?
    
    call getstr( 'GR', 0, c, l, ierr)
    if(ierr .eq. 0 .and. (c .eq. 'y' .or. c .eq. 'Y')) then
      gage_rain = .true.
    else
      gage_rain = .false. ! default
    end if

!   Uniform saturation for all elements (ignore rain gage values)

    call getr4('SA', 0, usat, ierr)
    if (ierr .eq. 0) then
      unif_sat  = .true.
    else
      unif_sat = .false. ! default
    end if

!   Monitor overtopping of channels and ponds

    xpond_n = 0
    xchan_n = 0

!   Need temperature (C) to compute kinematic viscosity for laminar flow option in plane
    
    call getr4('T', 0, temp, ierr)
    if(ierr .gt. 0) temp = 0.
    if(units .eq. 2) temp = (temp - 32.) * 5. / 9.
    xnu = .0000017756 / (1.+ 0.03368 * temp + 0.000221 * temp*temp) ! m**2/s
    if(units .eq. 2) xnu = xnu * 10.76496 ! convert to ft**2/s
    
    if(sed) call sed00()
    
    call wrt00(files(3))
    
    call clerk(limit)
    
    call rain(rfile)
    
    do j = 1, 4
      dtm(j) = delt
    end do

    nelt = 0
    ltab = 0
    nrain = 0
    imperv = 0.
    bed = 0.
    bedvol = 0.
    
!TEMP 21AUG2013
!    NPLANES = 53
!    iplane = 0
!    Hav = 0.
!TEMP 21AUG2013

!   Element processing loop...

    do
    
      call reader(files(1), block_name, ierr)
    
      if(ierr .gt. 0) then
        if(ierr .eq. 1) exit ! end of file
        call errxit('K2shell', "error reading parameter file")                                                                   
      end if
     
      nprint = 0
      nchan  = 1
      nelt   = nelt + 1 ! these will not agree when there is a
      ltab   = ltab + 1 ! compound channel with overbank infiltration
    
!     Print = 3 option
    
      call geti4 ('PRI', 0, nprint, ierr)
    
      if(nprint .eq. 3) then
        call getstr ('FI', 0, auxfile, j, ierr)
        if (ierr .gt. 0) then
          call errxit('K2shell', "PRINT=3, but no file specified")                                                                   
        else
          open (files(4), file = auxfile, status = 'unknown', iostat = file_error)
          if (file_error .ne. 0) call errxit('K2shell', "can't open "//auxfile)                                                                   
        end if
      end if
    
!     Profiler
    
      profiler = .FALSE.                                                                   
      if (block_name(1:5) == 'PLANE' .or. block_name(1:7) == 'CHANNEL') then
        call getstr ('PRO', 0, auxfile, j, ierr)
        if (ierr .eq. 0) then
          open (789, file = auxfile, status = 'unknown', iostat = file_error)
          if (file_error .eq. 0) profiler = .TRUE.                                                                   
        end if
      end if
    
!     Report progress

      call geti4 ('ID', 0, id, ierr)

      if(.not. silent) write (*, "('+ Processing ', a8, i8)") block_name, id

!     Rainfall & initial soil moisture

      call getr4 ('SA', 0, sat, ierr)

      if (ierr .gt. 0) sat = -1.

      if(block_name(1:5) == 'PLANE'                 .or. &
         block_name(1:7) == 'CHANNEL'               .or. &
        (block_name(1:4) == 'POND' .and. pond_rain) .or. &
         block_name(1:5) == 'URBAN') then ! Rainfall
    
        call geti4 ('RA', 0, ngage, ierr)

        if(ierr .gt. 0) then
          if(gage_rain) then
            if(block_name(1:7) == 'CHANNEL')then
              if(chrain) then
                nrain = nrain + 1
                ngage = nrain
              else
                ngage = 0
              end if
            else
              nrain = nrain + 1
              ngage = nrain
            end if
          else
            ngage = 0
            call getr8 ('X', 0, x, ierr)
            if(ierr .gt. 0) x = 0.
            call getr8 ('Y', 0, y, ierr)
            if(ierr .gt. 0) y = 0.
          end if
        end if

        if(block_name(1:7) == 'CHANNEL')then
!          if(chrain .AND..NOT. gage_rain) then
          if(.NOT. gage_rain) then
            call interp (x, y, t, z, nd, sat, ngage)
          end if
        else
          call interp (x, y, t, z, nd, sat, ngage)
        end if

        do i = 1, nd
          z(i) = z(i) * RMULT
        end do

! TEMP SNOW
        IF( snow )THEN
          limit = INT( t(nd) * 2. / delt ) + 1 ! t, delt in seconds
          CALL clerk(limit)
          CALL getr4('SM', 1, sm, ierr)
          CALL getr4('G', 1, g, ierr)
          CALL getr4('PO', 1, p, ierr)
          vr = (0.1 + 0.3/(1.+(1./g)**4.)**.25) * p
          vm = sm * p
          vi = sat * vm
          sat = (vi - vr) / (vm - vr)
          IF( sat < 0. ) sat = 0.
!          sat = -1.
        END IF
! TEMP SNOW

      end if

!     Initial saturation

      if(unif_sat) sat = usat
    
!     Multipliers

      rmks  = rmks_ 
      rmn   = rmn_  
      rmcv  = rmcv_ 
      rmg   = rmg_  
      rmin  = rmin_ 
      rmcoh = rmcoh_
      rmspl = rmspl_

      IF( nlist > 0 )THEN
        DO i = 1, nlist
          IF( mlist(i) == id )THEN
            rmks  = rmks_alt
            rmn   = rmn_alt
            rmcv  = rmcv_alt
            rmg   = rmg_alt
            rmin  = rmin_alt
            rmcoh = rmcoh_alt
            rmspl = rmspl_alt
            EXIT
          END IF
        END DO
      END IF

      IF( extra_mults .AND. block_name(1:7) == 'CHANNEL' )THEN
        rmks = rmks_chan
        rmg  = rmg_chan
        rmn  = rmn_chan
      END IF

!     Elements

      call_writer = .TRUE.

      if(block_name(1:5) .eq. 'PLANE') then
!TEMP 21AUG2013
!        iplane = iplane + 1
!TEMP 21AUG2013
        call plane(t, z, nd, sat)
v_plane = v_plane + qbal(10)
f_plane = f_plane + qbal(6)
      else if(block_name(1:7) .eq. 'CHANNEL') then
        call channl(t, z, nd, sat)
f_chan = f_chan + qbal(7)
      else if(block_name(1:4) .eq. 'PIPE') then
        call pipe()                      
      else if(block_name(1:6) .eq. 'INJECT') then
        call inject()
      else if(block_name(1:4) .eq. 'POND') then
        call pond(t, z, nd)
f_pond = f_pond + qbal(7)
s_pond = s_pond + qbal(9)
      else if(block_name(1:5) .eq. 'URBAN') then                       
        call urban(t, z, nd, sat)
v_urban = v_urban + qbal(10)
f_urban = f_urban + qbal(6)
      else if(block_name(1:5) .eq. 'ADDER') then                       
        call adder()
!        call_writer = .FALSE.
      else if(block_name(1:8) .eq. 'DIVERTER') then                       
        call diverter()
        call_writer = .FALSE.
      else
        call errxit(block_name, 'invalid element')
      end if
      
      if(call_writer) then

!       Add element contribution to global volume balance components
    
        do i = 1, 9
          sumbal(i) = sumbal(i) + qbal(i)
        end do
    
        call writer(files(4))
    
      end if

      if (profiler) close(789)

!   End of element processing loop
    
    end do

!TEMP 21AUG2013
!    WRITE( 666, '(2I6)' ) NPLANES, limit
!    WRITE( 666, '(1000I8)' ) (IDs(j), j = 1, NPLANES)
!    DO i = 1, limit
!      WRITE( 666, '(1000F8.2)' ) (Hav(j,i)*1000., j = 1, NPLANES) ! mm
!    END DO
!TEMP 21AUG2013

!    WRITE( 44, '(A3,F16.4)' ) doy, sumbal(2)
    
!   Outflow from final element
    
    sumbal(10) = qbal(10)
    
!   Write  event summary
    
    call event(sumbal)

    do i = 1, 4
      close(files(i))
    end do

    CALL DATE_AND_TIME( values = tvals )

    te = tvals(5) * 3600 + tvals(6) * 60 + tvals(7) - t0

    write (*, "(' Elapsed time (sec) ', i3.3)") te

!   End of simulation loop

    if(mode == INTERACTIVE) exit

! write(99, '(4F16.2)') v_urban, v_plane, bed, bedvol
! write(99, '(7F16.0)') v_urban, f_urban, v_plane, f_plane, s_pond, f_pond, f_chan

! TEMP SNOW
1976 CONTINUE
! TEMP SNOW

  end do
  
end program K2shell


logical function char2log(char)

! Converts character values Y|N or y|n into true|false
! Returns char as uppercase

  character(LEN=1) :: char

  if(char .eq. 'Y' .or. char .eq. 'y') then
    char2log = .true.
    char = 'Y'
  else
    char2log = .false.
    char = 'N'
  end if
end

subroutine errxit (id, message)

! Write error message to stdout & output file

   use runpars

   character(len=*) id, message

   write(files(3), 10) id, trim(message)
   write(*, 10)        id, trim(message)
10 format(//, ' error - ', a, /, 1x, a)
   stop ' '
end
