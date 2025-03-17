! Description:
!
!   Diverts water and sediment to as many as 10 downstream elements
!
!-------------------------------------------------------------------------------
!
! Global parameters:
!
!   sed                    log    .true. = route sediment (module runpars)
!   nps                    int     number of particle classes (module runpars)
!   limit                  int     number of time steps (module runpars)
!
!-------------------------------------------------------------------------------
!
! Input Block (input file labels in parenthesis):
!
!   id    (ID)             int*4  identifier (up to 6 digits)
!   idup  (UP)             int*4  identifier of upstream contributing element
!   idiv  (DIV)            int*4  up to 10 identifiers for diverted flow
!
!   qin   (IN)             real   tabulated discharge values of inflow, and
!   qdiv  (DIV1,DIV2,...)  real   corresponding values of flows diverted
!
!   nrows (N)              int*4  number rows in table
!
!-------------------------------------------------------------------------------
!
! Subroutines/functions:
!
!   getr4         (reader.for)
!   geti4              "
!   new           (miscel.for/clerk)
!   old                "
!   get                "
!   store              "
!   geta               "
!   stora              "
!   errxit        (errxit.for)
!-------------------------------------------------------------------------------

  subroutine diverter()

    use runpars
    use elpars
  
    implicit none
  
    integer i, j, k, n 
    integer ierr, ndivs, nrows 
    integer idup, idiv(10) 
    integer upr, outr, divr(10) 
    integer upq, outq, divq(10) 
    integer upc(5), outc(5), divc(10,5) 
  
    real qin(1000), qdiv(1000,10), conc(5)
    real coarea, qr, qu, qd, qsum
  
    character(len=15) msg
    character(len=5)  div(10)
    character(len=2) attr(5)
  
    data div/'DIV01','DIV02','DIV03','DIV04','DIV05','DIV06','DIV07','DIV08','DIV09','DIV10'/
    data attr/'S1','S2','S3','S4','S5'/
  
    qbal = 0.
  
    call geti4 ('ID', 0, id, ierr)
    if (ierr .gt. 0) call errxit ('DIVERTER', 'missing/invalid identifier (ID)')
  
    write (msg, '("DIVERTER ",i6)') id
  
    call geti4 ('UP', 0, idup, ierr)
    if (ierr .gt. 0) call errxit (msg, 'missing/invalid upstream identifier (ID)')
  
    ndivs = 0
    do k = 1, 10
      call geti4 ('DIV', k, idiv(k), ierr)
      if (ierr .eq. 3) call errxit(msg, 'invalid diversion identifier (DIV)')
      if (ierr .eq. 0) then
        ndivs = ndivs + 1
      else
        exit
      end if
    end do
  
    call geti4 ('N', 0, nrows, ierr)
    if (ierr .gt. 0) call errxit(msg, 'number of table rows (N) missing or invalid')
    if (nrows .gt. 999) call errxit(msg, 'too many table rows')
  
    do  j = 1, nrows
      call getr4 ('IN', j, qin(j), ierr)
      if (ierr .gt. 0) call errxit(msg, 'inflow (IN) missing or invalid')
      do k = 1, ndivs
        call getr4 (div(k), j, qdiv(j,k), ierr)
        if (ierr .gt. 0) call errxit(msg, 'diversion ('//div(k)//' missing or invalid')
      end do
    end do
  
    upr = 0
    call old (idup, 'RA', upr, n, ierr)
  
    call old (idup, 'RO', upq, n, ierr)
    if (ierr .gt. 0) call errxit(msg, 'upstream inflow not found')
  
    call geta (idup, coarea)
    call stora (id, coarea)
    do k = 1, ndivs
      call stora (idiv(k), 0.) ! zero area for diversions
    end do
  
    if (sed) then
      do k = 1, nps
        call old (idup, attr(k), upc(k), n, ierr)
        if (ierr .gt. 0) call errxit(msg, 'upstream sediment concentration not found')
      end do
    end if
  
    call new (id, 'RA', outr)
    call new (id, 'RO', outq)
  
    if (sed) then
      do n = 1, nps
        call new (id, attr(n), outc(n))
      end do
    end if
  
    do k = 1, ndivs
      call new (idiv(k), 'RA', divr(k))
      call new (idiv(k), 'RO', divq(k))
      if (sed) then
        do n = 1, nps
          call new (idiv(k), attr(n), divc(k,n))
        end do
      end if
    end do
  
    qr = 0.
  
    do i = 1, limit - 1
  
      if (upr .gt. 0)  call get (upr, i, qr)
  
      call get (upq, i, qu)
      if (sed) then
        do n = 1, nps
          call get (upc(n), i, conc(n))
        end do
      end if
  
      if (qu .ge. qin(nrows)) then
        j = nrows
      else
        do j = nrows-1, 1, -1
          if (qin(j) < qu) exit
        end do
      end if
  
      qsum = 0.
      do k = 1, ndivs
        if (j .eq. nrows) then
          qd = qdiv(nrows,k)
        else if (j .eq. 0) then
          qd = 0.
        else
          qd = qdiv(j,k) + ((qu - qin(j))/(qin(j+1) - qin(j)))*(qdiv(j+1,k) - qdiv(j,k))
        end if
        call store (divq(k), i, qd)
        qsum = qsum + qd
      end do
  
      call store (outq, i, qu-qsum)
  
      call store (outr, i, qr)
      do k = 1, ndivs
        call store (divr(k), i, 0.) ! zero rainfall for diversions
      end do

      if (sed) then
        do n = 1, nps
          call store (outc(n), i, conc(n))
          do k = 1, ndivs
            call store (divc(k,n), i, conc(n))
          end do
        end do
      end if

    end do

  end
