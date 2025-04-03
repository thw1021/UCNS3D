module communications
  ! @brief
  ! this module includes all the subroutines related to the mpi communications for various procedures
  ! including exchange of halo cells from stencils, direct-side neighbours, and the boundary extrapolated
  ! values for the conserved variables and their gradients at each surface/edge gaussian quadrature point
  use mpi_info
  use declaration
  use transform
  implicit none

contains
  subroutine renumber_neighbours(n, ielem, xmpie, xmpielrank, iexchanger, iexchanges)
    ! @brief
    ! this subroutine renumbers the neighbours indexing for cross referencing between different cpus
    ! it is a process that is performed once the beginning of each run
    implicit none
    type(element_number), allocatable, dimension(:, :), intent(inout)::ielem
    integer, intent(in)::n
    integer, allocatable, dimension(:), intent(in)::xmpie
    integer, allocatable, dimension(:), intent(in)::xmpielrank
    type(exchange), allocatable, dimension(:), intent(inout)::iexchanger, iexchanges
    integer::i,j,k,l,m,e,kmaxe,ineedt,tneedt,ix1,ix2,cnbt,ivt,itax,inum_points,inn,jjj,kk,itogg,ifdn,ifdn2,ifdn3,iii,iouf,jj1
    integer::c_n1, c_n2, c_n3, c_n4, d_n1, d_n2, d_n3, d_n4, kvf, itor4
    ineedt = iexchanger(1)%tot
    tneedt = iexchanges(1)%tot
    kmaxe = xmpielrank(n)
    do i = 1, kmaxe
      if (ielem(n, i)%interior .eq. 1) then
        allocate(ielem(n, i)%ineighn(ielem(n, i)%ifca))
        allocate(ielem(n, i)%ineighb(ielem(n, i)%ifca))
        allocate(ielem(n, i)%ineigh(ielem(n, i)%ifca))
        ielem(n, i)%ineigh(:) = 0; ielem(n, i)%ineighn(:) = 0
        ielem(n, i)%ineighb(:) = n
      else
        allocate(ielem(n, i)%ineighn(ielem(n, i)%ifca))
        allocate(ielem(n, i)%ineigh(ielem(n, i)%ifca))
        ielem(n, i)%ineighn(:) = 0
        ielem(n, i)%ineigh(:) = 0
      end if
    end do

    if (dimensiona .eq. 3) then
      ifdn = 4
    else
      ifdn = 2
    end if

    do i = 1, kmaxe
      if (ielem(n, i)%interior .eq. 0) then
        do j = 1, ielem(n, i)%ifca
          if (ielem(n, i)%ineighg(j) .gt. 0) then
            if (xmpie(ielem(n, i)%ineighg(j)) .eq. n) then
              k = xmpil(ielem(n, i)%ineighg(j))
              do ivt = 1, ielem(n, k)%ifca
                if (ielem(n, k)%ineighg(ivt) .eq. ielem(n, i)%ihexgl) then
                  ielem(n, i)%ineigh(j) = xmpil(ielem(n, i)%ineighg(j))
                  ielem(n, i)%ineighn(j) = ivt
                  go to 101
                end if
              end do
            end if
          end if
101       continue
        end do
      end if
    end do

    do i = 1, kmaxe
      if (ielem(n, i)%interior .eq. 1) then
        do j = 1, ielem(n, i)%ifca
          itax = 0
          if ((ielem(n, i)%ineighg(j) .gt. 0)) then
            if (xmpie(ielem(n, i)%ineighg(j)) .ne. n) then
              ielem(n, i)%ineighb(j) = xmpie(ielem(n, i)%ineighg(j))
            else
              k = xmpil(ielem(n, i)%ineighg(j))
              do ivt = 1, ielem(n, k)%ifca
                if (ielem(n, k)%ineighg(ivt) .eq. ielem(n, i)%ihexgl) then
                  ielem(n, i)%ineigh(j) = k
                  ielem(n, i)%ineighn(j) = ivt
                  ielem(n, i)%ineighb(j) = n
                end if
              end do
            end if
          end if
        end do
      end if
    end do

    jj1 = 0
    do k = 1, ineedt
      do ix1 = 1, tneedt
        do e = 1, iexchanger(k)%muchineed(1)
          do ix2 = 1, iexchanges(ix1)%muchtheyneed(1)
            if ((iexchanges1(ix1)%whattheyneed(ix2) .eq. iexchanger1(k)%sideineedn(e)) .and. &
            (iexchanges1(ix1)%sidetheyneedn(ix2) .eq. iexchanger1(k)%whatineed(e))) then
              i = xmpil(iexchanger1(k)%sideineedn(e))
              if (ielem(n, i)%interior .eq. 1) then
                do j = 1, ielem(n, i)%ifca
                  if (ielem(n, i)%ineighb(j) .ne. n) then
                    if (ielem(n, i)%ibounds(j) .gt. 0) then
                      if ((ibound(n, ielem(n, i)%ibounds(j))%icode .eq. 5) .or. (ibound(n, ielem(n, i)%ibounds(j))%icode .eq. 50)) then
                        if (dimensiona .eq. 3) then
                          if (ielem(n, i)%types_faces(j) .eq. 5) then
                            inum_points = qp_quad_n
                          else
                            inum_points = qp_triangle_n
                          end if
                        else
                          inum_points = qp_line_n
                        end if
                        if ((iexchanger(k)%procid .eq. ielem(n, i)%ineighb(j)) .and. (iexchanges(ix1)%procid .eq. ielem(n, i)%ineighb(j))) then
                          if ((iexchanges1(ix1)%sidetheyneedn(ix2) .eq. ielem(n, i)%ineighg(j)) .and. &
                          (ielem(n, i)%ineighg(j) .eq. iexchanger1(k)%whatineed(e))) then
                            do inn = 1, inum_points
                              if ((iexchanges1(ix1)%whattheyneed(ix2).eq.ielem(n,i)%ihexgl).and.(iexchanges1(ix1)%qtheyneed(ix2).eq.iexchanger1(k)%qineed(e)).and.(inn.eq.iexchanges1(ix1)%qtheyneed(ix2)))then
                                ielem(n, i)%q_face(j)%q_mapl(inn) = e
                                ielem(n, i)%ineighn(j) = k
                                iexchanges(ix1)%sidetheyneed(ix2) = j
                                jj1 = jj1 + 1
                              end if
                            end do
                          end if
                        end if
                      end if
                    else
                      if (dimensiona .eq. 3) then
                        if (ielem(n, i)%types_faces(j) .eq. 5) then
                          inum_points = qp_quad_n
                        else
                          inum_points = qp_triangle_n
                        end if
                      else
                        inum_points = qp_line_n
                      end if
                      if ((iexchanger(k)%procid .eq. ielem(n, i)%ineighb(j)) .and. (iexchanges(ix1)%procid .eq. ielem(n, i)%ineighb(j))) then
                        if ((iexchanges1(ix1)%sidetheyneedn(ix2) .eq. ielem(n, i)%ineighg(j)) .and. &
                        (ielem(n, i)%ineighg(j) .eq. iexchanger1(k)%whatineed(e))) then
                          do inn = 1, inum_points
                            if ((iexchanges1(ix1)%whattheyneed(ix2).eq.ielem(n,i)%ihexgl).and.(iexchanges1(ix1)%qtheyneed(ix2).eq.iexchanger1(k)%qineed(e)).and.(inn.eq.iexchanges1(ix1)%qtheyneed(ix2)))then
                              ielem(n, i)%q_face(j)%q_mapl(inn) = e
                              ielem(n, i)%ineighn(j) = k
                              iexchanges(ix1)%sidetheyneed(ix2) = j
                            end if
                          end do
                        end if
                      end if
                    end if
                  end if
                end do
              end if
            end if
          end do
        end do
      end do
    end do

  end subroutine renumber_neighbours

  subroutine solex_alloc(n)
    ! @brief
    ! this subroutine allocates the memory for the halo cells of the direct side neighbours only
    ! and is therefore used for lower-order schemes
    implicit none
    integer, intent(in)::n
    integer::i, j, k, l, m, e, kmaxe, ineedt, tneedt, icpuid
    real, dimension(1:dimensiona)::cords
    ineedt = iexchanger(1)%tot
    tneedt = iexchanges(1)%tot
    allocate(solchanger(ineedt))
    allocate(solchanges(tneedt))

    do i = 1, ineedt
      solchanger(i)%procid = iexchanger(i)%procid
      allocate(solchanger(i)%centres(iexchanger(i)%muchineed(1), dims))
      solchanger(i)%centres(:, :) = 0.d0
      allocate(solchanger(i)%sol(iexchanger(i)%muchineed(1), nof_variables + turbulenceequations + passivescalar))
      solchanger(i)%sol(:, :) = 0.0d0
    end do
    do i = 1, tneedt
      solchanges(i)%procid = iexchanges(i)%procid
      allocate(solchanges(i)%centres(iexchanges(i)%muchtheyneed(1), dims))
      solchanges(i)%centres(:, :) = 0.0d0
      allocate(solchanges(i)%sol(iexchanges(i)%muchtheyneed(1), nof_variables + turbulenceequations + passivescalar))
      solchanges(i)%sol(:, :) = 0.0d0
    end do

!-------------------for debugging only -----------------------------------------!
    do i = 1, tneedt
      do k = 1, iexchanges(i)%muchtheyneed(1)
        j = iexchanges(i)%localref(k)
        if (dimensiona .eq. 3) then
          call compute_centre3d(j, cords)
        else
          call compute_centre2d(j, cords)
        end if
        solchanges(i)%centres(k, 1:dims) = cords(1:dims)
      end do
    end do
    call mpi_barrier(mpi_comm_world, ierror)
    icpuid = n
    do i = 1, ineedt
      do k = 1, tneedt
        if (solchanger(i)%procid .eq. solchanges(k)%procid) then
          call mpi_sendrecv(solchanges(k)%centres(1:iexchanges(k)%muchtheyneed(1), 1:dims), &
                            iexchanges(k)%muchtheyneed(1)*dims, mpi_double_precision, solchanger(i)%procid, &
                            solchanges(k)%procid, solchanger(i)%centres(1:iexchanger(i)%muchineed(1), 1:dims), &
                            iexchanger(i)%muchineed(1)*dims, mpi_double_precision, &
                            solchanger(i)%procid, icpuid, mpi_comm_world, status, ierror)
        end if
      end do
    end do
  end subroutine solex_alloc

  subroutine estabexhange(n, ielem, imaxe, xmpie, xmpin, xmpielrank, ilocalstencil, iexchanger, &
                          iexchanges, irecexr, irecexs, numneighbours, ischeme, isize, iperiodicity, typesten, xmpil)
    ! @brief
    ! this subroutine is establishing the communication patterns for all the mpi processes and in particular for
    ! the halo cells of the reconstruction stencils and it must be noted that each cell can have a different number of stencils
    ! of different size

    implicit none
    type(element_number), allocatable, dimension(:, :), intent(inout)::ielem
    integer, intent(in)::n, imaxe, isize, iperiodicity, typesten
    integer, allocatable, dimension(:), intent(in)::xmpie, xmpin, xmpil
    integer, allocatable, dimension(:), intent(in)::xmpielrank
    type(exchange), allocatable, dimension(:), intent(inout)::iexchanges
    type(exchange), allocatable, dimension(:), intent(inout)::iexchanger
    integer, allocatable, dimension(:, :, :, :), intent(in)::ilocalstencil
    integer, intent(in)::numneighbours, ischeme
    type(recex), allocatable, dimension(:), intent(inout)::irecexr                !receive elements due to stencils
    type(recex), allocatable, dimension(:), intent(inout)::irecexs
    integer, allocatable, dimension(:)::totstenc, totstenc2
    integer, allocatable, dimension(:)::stenred
    integer, allocatable, dimension(:)::stenredproc
    type liststencils
      integer::procid
      integer, allocatable, dimension(:)::listsarray
    end type liststencils
    integer, dimension(0:isize - 1)::sendwh
    integer, dimension(0:isize - 1)::listofpr
    integer, dimension(0:isize - 1)::sendstwh
    type:: listgog
      integer::procid, imuch, iavc, iavt
      integer, allocatable, dimension(:)::globarray
      integer, allocatable, dimension(:, :)::nodex
    end type listgog
    type(liststencils), allocatable, dimension(:)::iliststen
    type(listgog), allocatable, dimension(:)::ilistgog, ilistside, ilistglo, ilistq
    integer::howmany, itee, source, tag, cpuidrec, icpuid, ixflag, iteedum, iavc, iavt
    integer, dimension(1:1)::rumts, sumts
    real, dimension(1:1)::dumts
    integer::stencounter, stencounter2, stensame, iiflag, tempsten2, iiflag2, kxk, ix, icx, icf, ifdn, ifdn2
    integer, dimension(0:isize - 1)::ihmste

    integer::i, j, ji, k, lm, kmaxn, kk, kmaxe, iaa, l, ingt, ik, tempint, temp2, itarget, isty5, i86, i87, i65, inum_points, iin
    kmaxe = xmpielrank(n)
    listofpr(:) = 0
    i65 = 0

    if (dimensiona .eq. 3) then
      ifdn = 4
    else
      ifdn = 2
    end if

    do k = 1, kmaxe
      if (ielem(n, k)%interior .eq. 1) then
        allocate(ielem(n, k)%q_face(ielem(n, k)%ifca))
        do l = 1, ielem(n, k)%ifca
          if (dimensiona .eq. 3) then
            if (ielem(n, k)%types_faces(l) .eq. 5) then
              inum_points = qp_quad_n
            else
              inum_points = qp_triangle_n
            end if
          else
            inum_points = qp_line_n
          end if
          allocate(ielem(n, k)%q_face(l)%q_mapl(inum_points))
          ielem(n, k)%q_face(l)%q_mapl(:) = 0
        end do
      end if
    end do

    do k = 1, kmaxe
      if (ielem(n, k)%interior .eq. 1) then
        do l = 1, ielem(n, k)%ifca
          if (ielem(n, k)%ineighg(l) .gt. 0) then
            if (xmpie(ielem(n, k)%ineighg(l)) .ne. n) then
              if (dimensiona .eq. 3) then
                if (ielem(n, k)%types_faces(l) .eq. 5) then
                  inum_points = qp_quad_n
                else
                  inum_points = qp_triangle_n
                end if
              else
                inum_points = qp_line_n
              end if
              listofpr(xmpie(ielem(n, k)%ineighg(l))) = listofpr(xmpie(ielem(n, k)%ineighg(l))) + inum_points
            end if
          end if
        end do
      end if
    end do

    howmany = 0
    do k = 0, isize - 1
      if (listofpr(k) .gt. 0) then
        howmany = howmany + 1
      end if
    end do
    allocate(ilistgog(howmany))
    allocate(ilistside(howmany))
    allocate(ilistq(howmany))
    allocate(ilistglo(howmany))

    kk = 0
    do k = 0, isize - 1
      if (listofpr(k) .gt. 0) then
        kk = kk + 1
        ilistgog(kk)%procid = k
        ilistgog(kk)%imuch = listofpr(k)
        ilistside(kk)%procid = k
        ilistside(kk)%imuch = listofpr(k)
        ilistq(kk)%procid = k
        ilistq(kk)%imuch = listofpr(k)
        allocate(ilistgog(kk)%globarray(listofpr(k)))
        allocate(ilistside(kk)%globarray(listofpr(k)))
        allocate(ilistside(kk)%nodex(listofpr(k), ifdn))
        allocate(ilistglo(kk)%globarray(listofpr(k)))
        allocate(ilistq(kk)%globarray(listofpr(k)))
        ilistgog(kk)%globarray(:) = 0
        ilistside(kk)%globarray(:) = 0
        ilistside(kk)%nodex(:, :) = 0
        ilistglo(kk)%globarray(:) = 0
        ilistq(kk)%globarray(:) = 0
      end if
    end do

    kk = 0
    do ingt = 0, isize - 1
      if (listofpr(ingt) .gt. 0) then
        kk = kk + 1
        iaa = 0
        do k = 1, kmaxe
          if (ielem(n, k)%interior .eq. 1) then
            do l = 1, ielem(n, k)%ifca
              if (ielem(n, k)%ineighg(l) .gt. 0) then
                if (xmpie(ielem(n, k)%ineighg(l)) .ne. n) then
                  if (xmpie(ielem(n, k)%ineighg(l)) .eq. ilistgog(kk)%procid) then
                    if (dimensiona .eq. 3) then
                      if (ielem(n, k)%types_faces(l) .eq. 5) then
                        inum_points = qp_quad_n
                        ifdn2 = 4
                      else
                        inum_points = qp_triangle_n
                        ifdn2 = 3
                      end if
                    else
                      inum_points = qp_line_n
                      ifdn2 = 2
                    end if
                    do iin = iaa + 1, iaa + inum_points
                      ilistgog(kk)%globarray(iin) = (ielem(n, k)%ineighg(l))
                      ilistglo(kk)%globarray(iin) = (ielem(n, k)%ihexgl)
                      ilistside(kk)%globarray(iin) = l
                      ilistq(kk)%globarray(iin) = iin - iaa
                      ilistside(kk)%nodex(iin, 1:ifdn2) = ielem(n, k)%nodes_faces(l, 1:ifdn2)
                    end do
                    iaa = iaa + inum_points
                  end if
                end if
              end if
            end do
          end if
        end do
      end if
    end do

    sendwh(:) = 0
    allocate(iexchanger(howmany))
    allocate(iexchanger1(howmany))
    iexchanger(:)%tot = howmany
    iexchanger1(:)%tot = howmany

    kk = 0
    do k = 0, isize - 1
      if (listofpr(k) .gt. 0) then
        kk = kk + 1
        iexchanger(kk)%procid = k
      end if
    end do

    do i = 1, howmany
      if (iexchanger(i)%procid .ne. n) then
        allocate(iexchanger(i)%muchineed(1))
        iexchanger(i)%muchineed(:) = 0
        iexchanger(i)%muchineed(1) = listofpr(iexchanger(i)%procid)
        allocate(iexchanger1(i)%whatineed(iexchanger(i)%muchineed(1)))        !element number global
        allocate(iexchanger1(i)%localref(iexchanger(i)%muchineed(1)))
        allocate(iexchanger1(i)%sideineed(iexchanger(i)%muchineed(1)))
        allocate(iexchanger1(i)%nodex(iexchanger(i)%muchineed(1), ifdn))
        allocate(iexchanger1(i)%qineed(iexchanger(i)%muchineed(1)))
        allocate(iexchanger1(i)%sideineedn(iexchanger(i)%muchineed(1)))
        iexchanger1(i)%whatineed(:) = 0        !element number global
        iexchanger1(i)%localref(:) = 0
        iexchanger1(i)%sideineed(:) = 0
        iexchanger1(i)%nodex(:, :) = 0
        iexchanger1(i)%sideineedn(:) = 0
        iexchanger1(i)%qineed(:) = 0
      end if
    end do
    tempint = 0
    do i = 1, howmany
      if (iexchanger(i)%procid .ne. n) then
        tempint = tempint + 1
        if (iexchanger(i)%procid .eq. ilistgog(tempint)%procid) then
          iexchanger1(i)%whatineed(:) = ilistgog(tempint)%globarray(:)
          iexchanger1(i)%sideineedn(:) = ilistglo(tempint)%globarray(:)
          iexchanger1(i)%nodex(1:iexchanger(i)%muchineed(1), 1:ifdn) = ilistside(tempint)%nodex(1:iexchanger(i)%muchineed(1), 1:ifdn)
          iexchanger1(i)%qineed(1:iexchanger(i)%muchineed(1)) = ilistq(tempint)%globarray(1:iexchanger(i)%muchineed(1))
        end if
      end if
    end do
    deallocate(ilistgog)
    call mpi_barrier(mpi_comm_world, ierror)
    icpuid = n
    do i = 1, howmany
      if (iexchanger(i)%procid .ne. n) then
        call mpi_sendrecv(iexchanger(i)%muchineed(1), 1, mpi_integer, iexchanger(i)%procid, iexchanger(i)%procid, &
                          itee, 1, mpi_integer, iexchanger(i)%procid, icpuid, mpi_comm_world, status, ierror)
        sendwh(iexchanger(i)%procid) = itee
      end if
    end do
    call mpi_barrier(mpi_comm_world, ierror)
    temp2 = 0
    do i = 0, isize - 1
      if (sendwh(i) .gt. 0) then
        temp2 = temp2 + 1
      end if
    end do

    allocate(iexchanges(temp2))
    allocate(iexchanges1(temp2))
    iexchanges(:)%tot = temp2
    kk = 0
    do i = 0, isize - 1
      if (sendwh(i) .gt. 0) then
        kk = kk + 1
        iexchanges(kk)%procid = i
        allocate(iexchanges(kk)%muchtheyneed(1))
        iexchanges(kk)%muchtheyneed(:) = 0
        iexchanges(kk)%muchtheyneed(1) = sendwh(i)
        allocate(iexchanges1(kk)%whattheyneed(iexchanges(kk)%muchtheyneed(1)))
        allocate(iexchanges(kk)%localref(iexchanges(kk)%muchtheyneed(1)))
        allocate(iexchanges(kk)%sidetheyneed(iexchanges(kk)%muchtheyneed(1)))
        allocate(iexchanges1(kk)%sidetheyneedn(iexchanges(kk)%muchtheyneed(1)))
        allocate(iexchanges1(kk)%nodex(iexchanges(kk)%muchtheyneed(1), ifdn))
        allocate(iexchanges1(kk)%qtheyneed(iexchanges(kk)%muchtheyneed(1)))
        allocate(iexchanges(kk)%qtheyneed(iexchanges(kk)%muchtheyneed(1)))
        iexchanges1(kk)%whattheyneed(:) = 0
        iexchanges(kk)%localref(:) = 0
        iexchanges(kk)%sidetheyneed(:) = 0
        iexchanges1(kk)%sidetheyneedn(:) = 0
        iexchanges1(kk)%nodex(:, :) = 0
        iexchanges1(kk)%qtheyneed(:) = 0
        iexchanges(kk)%qtheyneed(:) = 0
      end if
    end do

    do i = 1, howmany
      do k = 1, temp2
        if (iexchanges(k)%procid .eq. ilistside(i)%procid) then
          iexchanges(k)%sidetheyneed(:) = ilistside(i)%globarray(:)
        end if
      end do
    end do
    call mpi_barrier(mpi_comm_world, ierror)

    do i = 1, howmany
      do k = 1, temp2
        if (iexchanger(i)%procid .eq. iexchanges(k)%procid) then
          call mpi_sendrecv(iexchanger1(i)%whatineed(1:iexchanger(i)%muchineed(1)), &
                          iexchanger(i)%muchineed(1), mpi_integer, iexchanger(i)%procid, &
                          iexchanger(i)%procid, iexchanges1(k)%whattheyneed(1:iexchanges(k)%muchtheyneed(1)), &
                          iexchanges(k)%muchtheyneed(1), mpi_integer, iexchanger(i)%procid, &
                          icpuid, mpi_comm_world, status, ierror)
        end if
      end do
    end do
    do i = 1, howmany
      do k = 1, temp2
        if (iexchanger(i)%procid .eq. iexchanges(k)%procid) then
          call mpi_sendrecv(iexchanger1(i)%nodex(1:iexchanger(i)%muchineed(1), 1:ifdn), &
                            iexchanger(i)%muchineed(1)*ifdn, mpi_integer, iexchanger(i)%procid, &
                            iexchanger(i)%procid, iexchanges1(k)%nodex(1:iexchanges(k)%muchtheyneed(1), 1:ifdn), &
                            iexchanges(k)%muchtheyneed(1)*ifdn, mpi_integer, iexchanger(i)%procid, &
                            icpuid, mpi_comm_world, status, ierror)
        end if
      end do
    end do
    do i = 1, howmany
      do k = 1, temp2
        if (iexchanger(i)%procid .eq. iexchanges(k)%procid) then
          call mpi_sendrecv(iexchanger1(i)%qineed(1:iexchanger(i)%muchineed(1)), &
                            iexchanger(i)%muchineed(1), mpi_integer, iexchanger(i)%procid, &
                            iexchanger(i)%procid, iexchanges1(k)%qtheyneed(1:iexchanges(k)%muchtheyneed(1)), &
                            iexchanges(k)%muchtheyneed(1), mpi_integer, iexchanger(i)%procid, &
                            icpuid, mpi_comm_world, status, ierror)
          iexchanges(k)%qtheyneed(1:iexchanges(k)%muchtheyneed(1)) = iexchanges1(k)%qtheyneed(1:iexchanges(k)%muchtheyneed(1))
        end if
      end do
    end do

    call mpi_barrier(mpi_comm_world, ierror)
    do i = 1, howmany
      do k = 1, temp2
        if (iexchanger(i)%procid .eq. iexchanges(k)%procid) then
          call mpi_sendrecv(iexchanger1(i)%sideineedn(1:iexchanger(i)%muchineed(1)), &
                            iexchanger(i)%muchineed(1), mpi_integer, iexchanger(i)%procid, &
                            iexchanger(i)%procid, iexchanges1(k)%sidetheyneedn(1:iexchanges(k)%muchtheyneed(1)), &
                            iexchanges(k)%muchtheyneed(1), mpi_integer, iexchanger(i)%procid, &
                            icpuid, mpi_comm_world, status, ierror)
        end if
      end do
    end do

    call mpi_barrier(mpi_comm_world, ierror)
    do i = 1, howmany
      do k = 1, temp2
        if (iexchanger(i)%procid .eq. iexchanges(k)%procid) then
          call mpi_sendrecv(iexchanges(k)%sidetheyneed(1:iexchanges(k)%muchtheyneed(1)), &
                            iexchanges(k)%muchtheyneed(1), mpi_integer, iexchanges(k)%procid, &
                            iexchanges(k)%procid, iexchanger1(i)%sideineed(1:iexchanger(i)%muchineed(1)), &
                            iexchanger(i)%muchineed(1), mpi_integer, iexchanges(k)%procid, &
                            icpuid, mpi_comm_world, status, ierror)
        end if
      end do
    end do

    call mpi_barrier(mpi_comm_world, ierror)

    deallocate(ilistside, ilistq)
    do i = 1, temp2
      do k = 1, iexchanges(i)%muchtheyneed(1)
        if (xmpie(iexchanges1(i)%whattheyneed(k)) .eq. n) then
          iexchanges(i)%localref(k) = xmpil(iexchanges1(i)%whattheyneed(k))
        end if
      end do
    end do
    do i = 1, howmany
      do k = 1, temp2
        if (iexchanger(i)%procid .eq. iexchanges(k)%procid) then
          call mpi_sendrecv(iexchanges(k)%localref(1:iexchanges(k)%muchtheyneed(1)), &
                            iexchanges(k)%muchtheyneed(1), mpi_integer, iexchanges(k)%procid, &
                            iexchanges(k)%procid, iexchanger1(i)%localref(1:iexchanger(i)%muchineed(1)), &
                            iexchanger(i)%muchineed(1), mpi_integer, iexchanger(i)%procid, &
                            icpuid, mpi_comm_world, status, ierror)
        end if
      end do
    end do
    call mpi_barrier(mpi_comm_world, ierror)

    if (ischeme .gt. 1) then
      stencounter = 0
      allocate(totstenc2(imaxe))
      totstenc2 = 0
      do i = 1, kmaxe
        do k = 1, typesten
          isty5 = 0
          if ((k .eq. 1) .or. (ees .ne. 5)) then
            itarget = ielem(n, i)%inumneighbours
          else
            itarget = numneighbours2
          end if
          do kk = 1, itarget
            if (ilocalstencil(n, i, k, kk) .gt. 0) then
              isty5 = isty5 + 1
            end if
          end do
          if (isty5 .eq. itarget) then
            do kk = 1, itarget
              iiflag = 0
              if (xmpie(ilocalstencil(n, i, k, kk)) .eq. n) then
                iiflag = 1
              end if
              if (iiflag .eq. 0) then
                totstenc2((ilocalstencil(n, i, k, kk))) = ilocalstencil(n, i, k, kk)
              end if
            end do
          end if
        end do
      end do

      kxk = 0
      do i = 1, imaxe
        if (totstenc2(i) .gt. 0) then
          kxk = kxk + 1
        end if
      end do

      allocate(stenred(kxk))
      allocate(stenredproc(kxk))
      stencounter = kxk
      allocate(totstenc(stencounter))
      totstenc = 0
      kxk = 0
      do i = 1, imaxe
        if (totstenc2(i) .gt. 0) then
          kxk = kxk + 1
          totstenc(kxk) = totstenc2(i)
        end if
      end do

      deallocate(totstenc2)
      stenred = 0
      stenredproc = 0
      ihmste = 0
      kxk = 0
      do i = 1, stencounter
        if (totstenc(i) .gt. 0) then
          kxk = kxk + 1
          stenred(kxk) = totstenc(i)
        end if
      end do
      deallocate(totstenc)

      do i = 1, kxk
        stenredproc(i) = xmpie(stenred(i))
      end do
      do i = 0, isize - 1
        ihmste(i) = 0
        if (i .ne. n) then
          do ix = 1, kxk
            if (stenredproc(ix) .eq. i) then
            ihmste(i) = ihmste(i) + 1
            end if
          end do
        end if
      end do
      ix = 0
      do i = 0, isize - 1
        if (i .ne. n) then
          if (ihmste(i) .gt. 0) then
            ix = ix + 1
          end if
        end if
      end do
      call mpi_barrier(mpi_comm_world, ierror)
      allocate(iliststen(ix))
      ix = 0
      do i = 0, isize - 1
        if (i .ne. n) then
          if (ihmste(i) .gt. 0) then
            ix = ix + 1
            iliststen(ix)%procid = i
            allocate(iliststen(ix)%listsarray(ihmste(i)))
            iliststen(ix)%listsarray(:) = 0
            icf = 0
            do icx = 1, kxk
              if (stenredproc(icx) .eq. i) then
                icf = icf + 1
                iliststen(ix)%listsarray(icf) = stenred(icx)
              end if
            end do
          end if
        end if
      end do
      allocate(irecexr(ix))
      allocate(irecexr1(ix))
      irecexr(:)%tot = 0
      irecexr(:)%tot = ix
      icx = 0
      do i = 0, isize - 1
        if (ihmste(i) .gt. 0) then
          icx = icx + 1
          irecexr(icx)%procid = i
          allocate(irecexr(icx)%muchineed(1))
          irecexr(icx)%muchineed(:) = 0
          irecexr(icx)%muchineed(1) = ihmste(i)
          allocate(irecexr1(icx)%whatineed(irecexr(icx)%muchineed(1)))
          allocate(irecexr1(icx)%localref(irecexr(icx)%muchineed(1)))
          allocate(irecexr1(icx)%ishape(irecexr(icx)%muchineed(1)))
          irecexr1(icx)%whatineed(:) = 0
          irecexr1(icx)%localref(:) = 0
          irecexr1(icx)%ishape(:) = 0
          do icf = 1, ix
            if (iliststen(icf)%procid .eq. irecexr(icx)%procid) then
              irecexr1(icx)%whatineed(:) = iliststen(icf)%listsarray(:)
            end if
          end do
        end if
      end do

      sendstwh(:) = 0
      icpuid = n
      do k = 0, isize - 1
        ixflag = 0
        iteedum = 0
        itee = 0
        if (k .ne. n) then
          do i = 1, icx
            if (irecexr(i)%procid .eq. k) then
              call mpi_sendrecv(irecexr(i)%muchineed(1), 1, mpi_integer, irecexr(i)%procid, icpuid, &
                                itee, 1, mpi_integer, irecexr(i)%procid, irecexr(i)%procid, mpi_comm_world, status, ierror)
              sendstwh(irecexr(i)%procid) = itee
              ixflag = 1
            end if
          end do
          if (ixflag .eq. 0) then
            call mpi_sendrecv(iteedum, 1, mpi_integer, k, icpuid, &
                              itee, 1, mpi_integer, k, k, mpi_comm_world, status, ierror)
            sendstwh(k) = itee
          end if
        end if
      end do
      call mpi_barrier(mpi_comm_world, ierror)
      temp2 = 0
      do i = 0, isize - 1
        if (sendstwh(i) .gt. 0) then
          temp2 = temp2 + 1
        end if
      end do

      call mpi_barrier(mpi_comm_world, ierror)
      allocate(irecexs(temp2))
      allocate(irecexs1(temp2))

      irecexs(:)%tot = temp2
      kk = 0
      do i = 0, isize - 1
        if (sendstwh(i) .gt. 0) then
          kk = kk + 1
          irecexs(kk)%procid = i
          allocate(irecexs(kk)%muchtheyneed(1))
          irecexs(kk)%muchtheyneed(:) = 0
          irecexs(kk)%muchtheyneed(1) = sendstwh(i)
          allocate(irecexs1(kk)%whattheyneed(irecexs(kk)%muchtheyneed(1)))
          allocate(irecexs(kk)%localref(irecexs(kk)%muchtheyneed(1)))
          allocate(irecexs1(kk)%ishape(irecexs(kk)%muchtheyneed(1)))
          irecexs1(kk)%whattheyneed(:) = 0
          irecexs(kk)%localref(:) = 0
          irecexs1(kk)%ishape(:) = 0
        end if
      end do

      rumts(1:1) = 0
      call mpi_barrier(mpi_comm_world, ierror)
      do i = 0, isize - 1
        if (i .ne. n) then
          do j = 1, temp2
            iavt = 10000
            if (irecexs(j)%procid .eq. i) then
              iavt = j
              go to 7001
            end if
          end do
7001      continue
          do k = 1, icx
            iavc = 10000
            if (irecexr(k)%procid .eq. i) then
              iavc = k
              go to 8001
            end if
          end do
8001      continue
          if ((iavc .eq. 10000) .and. (iavt .ne. 10000)) then
            call mpi_sendrecv(sumts(1:1), 1, mpi_integer, i, icpuid, &
                              irecexs1(iavt)%whattheyneed(1:irecexs(iavt)%muchtheyneed(1)), irecexs(iavt)%muchtheyneed(1), &
                              mpi_integer, irecexs(iavt)%procid, irecexs(iavt)%procid, mpi_comm_world, status, ierror)
          end if
          if ((iavc .eq. 10000) .and. (iavt .eq. 10000)) then
            call mpi_sendrecv(sumts(1:1), 1, mpi_integer, i, icpuid, &
                              rumts(1:1), 1, mpi_integer, i, i, mpi_comm_world, status, ierror)
          end if
          if ((iavc .ne. 10000) .and. (iavt .eq. 10000)) then
            call mpi_sendrecv(irecexr1(iavc)%whatineed(1:irecexr(iavc)%muchineed(1)), irecexr(iavc)%muchineed(1), &
                        mpi_integer, irecexr(iavc)%procid, icpuid, rumts(1:1), 1, mpi_integer, i, i, mpi_comm_world, status, ierror)
          end if
          if ((iavc .ne. 10000) .and. (iavt .ne. 10000)) then
            call mpi_sendrecv(irecexr1(iavc)%whatineed(1:irecexr(iavc)%muchineed(1)), irecexr(iavc)%muchineed(1), &
                              mpi_integer, irecexr(iavc)%procid, icpuid, &
                         irecexs1(iavt)%whattheyneed(1:irecexs(iavt)%muchtheyneed(1)), irecexs(iavt)%muchtheyneed(1), mpi_integer, &
                              irecexs(iavt)%procid, irecexs(iavt)%procid, mpi_comm_world, status, ierror)
          end if
        end if        ! i.ne.n
      end do
      call mpi_barrier(mpi_comm_world, ierror)

      do i = 1, temp2
        do k = 1, irecexs(i)%muchtheyneed(1)
          if (xmpie(irecexs1(i)%whattheyneed(k)) .eq. n) then
            j = xmpil(irecexs1(i)%whattheyneed(k))
            irecexs(i)%localref(k) = j
            irecexs1(i)%ishape(k) = ielem(n, j)%ishape
          end if
        end do
      end do
      call mpi_barrier(mpi_comm_world, ierror)
      do i = 0, isize - 1
        rumts = 0
        sumts = 0
        if (i .ne. n) then
          do j = 1, temp2
            iavt = 10000
            if (irecexs(j)%procid .eq. i) then
              iavt = j
              go to 9001
            end if
          end do
9001      continue
          do k = 1, icx
            iavc = 10000
            if (irecexr(k)%procid .eq. i) then
              iavc = k
              go to 10001
            end if
          end do
10001     continue

          if ((iavt .eq. 10000) .and. (iavc .ne. 10000)) then
            call mpi_sendrecv(sumts(1:1), 1, mpi_integer, i, icpuid, &
                              irecexr1(iavc)%localref(1:irecexr(iavc)%muchineed(1)), irecexr(iavc)%muchineed(1), &
                              mpi_integer, irecexr(iavc)%procid, irecexr(iavc)%procid, mpi_comm_world, status, ierror)
            call mpi_sendrecv(sumts(1:1), 1, mpi_integer, i, icpuid, &
                              irecexr1(iavc)%ishape(1:irecexr(iavc)%muchineed(1)), irecexr(iavc)%muchineed(1), &
                              mpi_integer, irecexr(iavc)%procid, irecexr(iavc)%procid, mpi_comm_world, status, ierror)
          end if
          if ((iavt .eq. 10000) .and. (iavc .eq. 10000)) then
            call mpi_sendrecv(sumts(1:1), 1, mpi_integer, i, icpuid, &
                              rumts(1:1), 1, mpi_integer, i, i, mpi_comm_world, status, ierror)
            call mpi_sendrecv(sumts(1:1), 1, mpi_integer, i, icpuid, &
                              rumts(1:1), 1, mpi_integer, i, i, mpi_comm_world, status, ierror)
          end if
          if ((iavt .ne. 10000) .and. (iavc .eq. 10000)) then
            call mpi_sendrecv(irecexs(iavt)%localref(1:irecexs(iavt)%muchtheyneed(1)), &
                              irecexs(iavt)%muchtheyneed(1), mpi_integer, irecexs(iavt)%procid, icpuid, &
                              rumts(1:1), 1, mpi_integer, i, i, mpi_comm_world, status, ierror)
            call mpi_sendrecv(irecexs1(iavt)%ishape(1:irecexs(iavt)%muchtheyneed(1)), &
                              irecexs(iavt)%muchtheyneed(1), mpi_integer, irecexs(iavt)%procid, icpuid, &
                              rumts(1:1), 1, mpi_integer, i, i, mpi_comm_world, status, ierror)
          end if
          if ((iavt .ne. 10000) .and. (iavc .ne. 10000)) then
            call mpi_sendrecv(irecexs(iavt)%localref(1:irecexs(iavt)%muchtheyneed(1)), &
                              irecexs(iavt)%muchtheyneed(1), mpi_integer, irecexs(iavt)%procid, icpuid, irecexr1(iavc)%localref(1:irecexr(iavc)%muchineed(1)), &
                              irecexr(iavc)%muchineed(1), mpi_integer, irecexr(iavc)%procid, irecexr(iavc)%procid, mpi_comm_world, status, ierror)
            call mpi_sendrecv(irecexs1(iavt)%ishape(1:irecexs(iavt)%muchtheyneed(1)), irecexs(iavt)%muchtheyneed(1), &
                              mpi_integer, irecexs(iavt)%procid, icpuid, irecexr1(iavc)%ishape(1:irecexr(iavc)%muchineed(1)), irecexr(iavc)%muchineed(1), &
                              mpi_integer, irecexr(iavc)%procid, irecexr(iavc)%procid, mpi_comm_world, status, ierror)
          end if
        end if        ! i.ne.n
      end do
      call mpi_barrier(mpi_comm_world, ierror)
      deallocate(stenred)
      deallocate(stenredproc)
      deallocate(iliststen)
    end if
    deallocate(ilistglo)
  end subroutine estabexhange

  subroutine exchange_lower(n)
    ! @brief
    ! this subroutine is exchanging the variables of all the halo cells for direct side neighbours only
    ! and is primarily used by lower order schemes
    implicit none
    integer, intent(in)::n
    integer::i, k, ineedt, tneedt, icpuid, itest
    ineedt = iexchanger(1)%tot
    tneedt = iexchanges(1)%tot
    itest = nof_variables + turbulenceequations + passivescalar

    if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
      do i = 1, tneedt
        do k = 1, iexchanges(i)%muchtheyneed(1)
          solchanges(i)%sol(k, 1:nof_variables) = u_c(iexchanges(i)%localref(k))%val(1, 1:nof_variables)
          solchanges(i)%sol(k,nof_variables+1:nof_variables+turbulenceequations+passivescalar)=u_ct(iexchanges(i)%localref(k))%val(1,1:turbulenceequations+passivescalar)
        end do
      end do
    else
      do i = 1, tneedt
        do k = 1, iexchanges(i)%muchtheyneed(1)
          solchanges(i)%sol(k, 1:nof_variables) = u_c(iexchanges(i)%localref(k))%val(1, 1:nof_variables)
        end do
      end do
    end if

    icpuid = n
    do i = 1, ineedt
      do k = 1, tneedt
        if (solchanger(i)%procid .eq. solchanges(k)%procid) then
          call mpi_sendrecv(solchanges(k)%sol(1:iexchanges(k)%muchtheyneed(1), 1:itest), &
                            iexchanges(k)%muchtheyneed(1)*itest, mpi_double_precision, solchanger(i)%procid, &
                            solchanges(k)%procid, solchanger(i)%sol(1:iexchanger(i)%muchineed(1), 1:itest), &
                            iexchanger(i)%muchineed(1)*itest, mpi_double_precision, &
                            solchanger(i)%procid, icpuid, mpi_comm_world, status, ierror)
        end if
      end do
    end do

  end subroutine exchange_lower

  subroutine exchange_higher(n)
    ! @brief
    ! this subroutine is exchanging the variables of all the halo cells for all the reconstruction stencils
    implicit none
    integer, intent(in)::n
    integer::i, j, k, l, m, o, p, q, ineedt, tneedt, indl, tndl, icpuid, itest, itee, iteedum, itemp1, itemp2, iavc, iavt
    real, dimension(1:1)::dumts, rumts
    integer:: n_requests
    integer, dimension(:), allocatable:: requests
    itest = nof_variables + turbulenceequations + passivescalar
    ineedt = irecexr(1)%tot
    tneedt = irecexs(1)%tot

    if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
      do i = 1, tneedt
        do k = 1, irecexs(i)%muchtheyneed(1)
          iexsolhis(i)%sol(k, 1:nof_variables) = u_c(irecexs(i)%localref(k))%val(1, 1:nof_variables)
                  iexsolhis(i)%sol(k,nof_variables+1:nof_variables+turbulenceequations+passivescalar)=u_ct(irecexs(i)%localref(k))%val(1,1:turbulenceequations+passivescalar)
        end do
      end do
    else
      do i = 1, tneedt
        do k = 1, irecexs(i)%muchtheyneed(1)
          iexsolhis(i)%sol(k,1:nof_variables+turbulenceequations+passivescalar)=u_c(irecexs(i)%localref(k))%val(1,1:nof_variables+turbulenceequations+passivescalar)
        end do
      end do
    end if

    n_requests = 0
    allocate(requests(jtotal*2))
    requests(:) = 0
    icpuid = n

    do k = 1, jtotal
      if ((jtot(k, 1) .eq. -1) .and. (jtot(k, 2) .ne. -1)) then
        n_requests = n_requests + 1
        iavc = jtot(k, 2)
        call mpi_isend( dumts(1:1), & !sendbuf
                        1, mpi_double_precision, & !sendcount, sendtype
                        jtot(k, 3), 0, & !destination, tag
                        mpi_comm_world, requests(n_requests), ierror & !communicator, request handle, error
                      )
        n_requests = n_requests + 1
        iavc = jtot(k, 2)
        call mpi_irecv( iexsolhir(iavc)%sol(1:irecexr(iavc)%muchineed(1), 1:itest), & !recvbuf
                        irecexr(iavc)%muchineed(1)*itest, mpi_double_precision, & !recvcount, recvtype
                        iexsolhir(iavc)%procid, 0, & !source, tag
                        mpi_comm_world, requests(n_requests), ierror & !communicator, request handle, error
                      )
      end if

      if ((jtot(k, 1) .ne. -1) .and. (jtot(k, 2) .eq. -1)) then
        n_requests = n_requests + 1
        iavt = jtot(k, 1)
        call mpi_isend( iexsolhis(iavt)%sol(1:irecexs(iavt)%muchtheyneed(1), 1:itest), & !sendbuf
                        irecexs(iavt)%muchtheyneed(1)*itest, mpi_double_precision, & !sendcount, sendtype
                        iexsolhis(iavt)%procid, 0, & !destination, tag
                        mpi_comm_world, requests(n_requests), ierror & !communicator, request handle, error
                      )
        n_requests = n_requests + 1
        iavt = jtot(k, 1)
        call mpi_irecv( dumts(1:1), & !recvbuf
                        1, mpi_double_precision, & !recvcount, recvtype
                        jtot(k, 3), 0, & !source, tag
                        mpi_comm_world, requests(n_requests), ierror & !communicator, request handle, error
                      )
      end if

      if ((jtot(k, 1) .ne. -1) .and. (jtot(k, 2) .ne. -1)) then
        n_requests = n_requests + 1
        iavt = jtot(k, 1)
        iavc = jtot(k, 2)
        call mpi_isend( iexsolhis(iavt)%sol(1:irecexs(iavt)%muchtheyneed(1), 1:itest), & !sendbuf
                        irecexs(iavt)%muchtheyneed(1)*itest, mpi_double_precision, & !sendcount, sendtype
                        iexsolhis(iavt)%procid, 0, & !destination, tag
                        mpi_comm_world, requests(n_requests), ierror & !communicator, request handle, error
                      )
        n_requests = n_requests + 1
        iavc = jtot(k, 2)
        call mpi_irecv( iexsolhir(iavc)%sol(1:irecexr(iavc)%muchineed(1), 1:itest), & !recvbuf
                        irecexr(iavc)%muchineed(1)*itest, mpi_double_precision, & !recvcount, recvtype
                        iexsolhir(iavc)%procid, 0, & !source, tag
                        mpi_comm_world, requests(n_requests), ierror & !communicator, request handle, error
                      )
      end if
    end do
    call mpi_waitall(n_requests, requests, mpi_statuses_ignore, ierror)
    deallocate(requests)
  end subroutine exchange_higher

  subroutine exchange_adda_diss(n)
    ! @brief
    ! this subroutine is exchanging the variables of all the halo cells for all the reconstruction stencils
    implicit none
    integer, intent(in)::n
    integer::i, j, k, l, m, o, p, q, ineedt, tneedt, indl, tndl, icpuid, itest, itee, iteedum, itemp1, itemp2, iavc, iavt
    real, dimension(1:1)::dumts, rumts
    integer:: n_requests
    integer, dimension(:), allocatable:: requests
    itest = 1
    ineedt = irecexr(1)%tot
    tneedt = irecexs(1)%tot

    do i = 1, tneedt
      do k = 1, irecexs(i)%muchtheyneed(1)
        iexsolhisd(i)%sol(k, 1) = ielem(n, irecexs(i)%localref(k))%diss
      end do
    end do

    n_requests = 0
    allocate(requests(jtotal*2))
    requests(:) = 0
    icpuid = n

    do k = 1, jtotal
      if ((jtot(k, 1) .eq. -1) .and. (jtot(k, 2) .ne. -1)) then
        n_requests = n_requests + 1
        iavc = jtot(k, 2)
        call mpi_isend( dumts(1:1), & !sendbuf
                        1, mpi_double_precision, & !sendcount, sendtype
                        jtot(k, 3), 0, & !destination, tag
                        mpi_comm_world, requests(n_requests), ierror & !communicator, request handle, error
                      )
        n_requests = n_requests + 1
        iavc = jtot(k, 2)
        call mpi_irecv( iexsolhird(iavc)%sol(1:irecexr(iavc)%muchineed(1), 1), & !recvbuf
                        irecexr(iavc)%muchineed(1)*itest, mpi_double_precision, & !recvcount, recvtype
                        iexsolhir(iavc)%procid, 0, & !source, tag
                        mpi_comm_world, requests(n_requests), ierror & !communicator, request handle, error
                      )
      end if
      if ((jtot(k, 1) .ne. -1) .and. (jtot(k, 2) .eq. -1)) then
        n_requests = n_requests + 1
        iavt = jtot(k, 1)
        call mpi_isend( iexsolhisd(iavt)%sol(1:irecexs(iavt)%muchtheyneed(1), 1), & !sendbuf
                        irecexs(iavt)%muchtheyneed(1)*itest, mpi_double_precision, & !sendcount, sendtype
                        iexsolhis(iavt)%procid, 0, & !destination, tag
                        mpi_comm_world, requests(n_requests), ierror & !communicator, request handle, error
                      )
        n_requests = n_requests + 1
        iavt = jtot(k, 1)
        call mpi_irecv( dumts(1:1), & !recvbuf
                        1, mpi_double_precision, & !recvcount, recvtype
                        jtot(k, 3), 0, & !source, tag
                        mpi_comm_world, requests(n_requests), ierror & !communicator, request handle, error
                      )
      end if
      if ((jtot(k, 1) .ne. -1) .and. (jtot(k, 2) .ne. -1)) then
        n_requests = n_requests + 1
        iavt = jtot(k, 1)
        iavc = jtot(k, 2)
        call mpi_isend( iexsolhisd(iavt)%sol(1:irecexs(iavt)%muchtheyneed(1), 1:itest), & !sendbuf
                        irecexs(iavt)%muchtheyneed(1)*itest, mpi_double_precision, & !sendcount, sendtype
                        iexsolhis(iavt)%procid, 0, & !destination, tag
                        mpi_comm_world, requests(n_requests), ierror & !communicator, request handle, error
                      )
        n_requests = n_requests + 1
        iavc = jtot(k, 2)
        call mpi_irecv( iexsolhird(iavc)%sol(1:irecexr(iavc)%muchineed(1), 1:itest), & !recvbuf
                        irecexr(iavc)%muchineed(1)*itest, mpi_double_precision, & !recvcount, recvtype
                        iexsolhir(iavc)%procid, 0, & !source, tag
                        mpi_comm_world, requests(n_requests), ierror & !communicator, request handle, error
                      )
      end if
    end do
    call mpi_waitall(n_requests, requests, mpi_statuses_ignore, ierror)
    deallocate(requests)
  end subroutine exchange_adda_diss

  subroutine exchange_higher_av(n)
    ! @brief
    ! this subroutine is exchanging the averaged variables of all the halo cells for all the reconstruction stencils
    implicit none
    integer, intent(in)::n
    integer::i, j, k, l, m, o, p, q, ineedt, tneedt, indl, tndl, icpuid, itest, itee, iteedum, itemp1, itemp2, iavc, iavt
    real, dimension(1:1)::dumts, rumts
    integer:: n_requests, ind1
    integer, dimension(:), allocatable:: requests
    if (rungekutta .eq. 4) then
      ind1 = 7
    else
      ind1 = 5
    end if

    itest = nof_variables + turbulenceequations + passivescalar
    ineedt = irecexr(1)%tot
    tneedt = irecexs(1)%tot

    if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
      do i = 1, tneedt
        do k = 1, irecexs(i)%muchtheyneed(1)
          iexsolhis(i)%sol(k, 1:nof_variables) = u_c(irecexs(i)%localref(k))%val(ind1, 1:nof_variables)
                  iexsolhis(i)%sol(k,nof_variables+1:nof_variables+turbulenceequations+passivescalar)=u_ct(irecexs(i)%localref(k))%val(ind1,1:turbulenceequations+passivescalar)
        end do
      end do
    else
      do i = 1, tneedt
        do k = 1, irecexs(i)%muchtheyneed(1)
                  iexsolhis(i)%sol(k,1:nof_variables+turbulenceequations+passivescalar)=u_c(irecexs(i)%localref(k))%val(ind1,1:nof_variables+turbulenceequations+passivescalar)
        end do
      end do
    end if

    n_requests = 0
    allocate(requests(jtotal*2))
    requests(:) = 0
    icpuid = n

    do k = 1, jtotal
      if ((jtot(k, 1) .eq. -1) .and. (jtot(k, 2) .ne. -1)) then
        n_requests = n_requests + 1
        iavc = jtot(k, 2)
        call mpi_isend( dumts(1:1), & !sendbuf
                        1, mpi_double_precision, & !sendcount, sendtype
                        jtot(k, 3), 0, & !destination, tag
                        mpi_comm_world, requests(n_requests), ierror & !communicator, request handle, error
                      )
        n_requests = n_requests + 1
        iavc = jtot(k, 2)
        call mpi_irecv( iexsolhir(iavc)%sol(1:irecexr(iavc)%muchineed(1), 1:itest), & !recvbuf
                        irecexr(iavc)%muchineed(1)*itest, mpi_double_precision, & !recvcount, recvtype
                        iexsolhir(iavc)%procid, 0, & !source, tag
                        mpi_comm_world, requests(n_requests), ierror & !communicator, request handle, error
                      )
      end if

      if ((jtot(k, 1) .ne. -1) .and. (jtot(k, 2) .eq. -1)) then
        n_requests = n_requests + 1
        iavt = jtot(k, 1)
        call mpi_isend( iexsolhis(iavt)%sol(1:irecexs(iavt)%muchtheyneed(1), 1:itest), & !sendbuf
                        irecexs(iavt)%muchtheyneed(1)*itest, mpi_double_precision, & !sendcount, sendtype
                        iexsolhis(iavt)%procid, 0, & !destination, tag
                        mpi_comm_world, requests(n_requests), ierror & !communicator, request handle, error
                      )
        n_requests = n_requests + 1
        iavt = jtot(k, 1)
        call mpi_irecv( dumts(1:1), & !recvbuf
                        1, mpi_double_precision, & !recvcount, recvtype
                        jtot(k, 3), 0, & !source, tag
                        mpi_comm_world, requests(n_requests), ierror & !communicator, request handle, error
                      )
      end if

      if ((jtot(k, 1) .ne. -1) .and. (jtot(k, 2) .ne. -1)) then
        n_requests = n_requests + 1
        iavt = jtot(k, 1)
        iavc = jtot(k, 2)
        call mpi_isend( iexsolhis(iavt)%sol(1:irecexs(iavt)%muchtheyneed(1), 1:itest), & !sendbuf
                        irecexs(iavt)%muchtheyneed(1)*itest, mpi_double_precision, & !sendcount, sendtype
                        iexsolhis(iavt)%procid, 0, & !destination, tag
                        mpi_comm_world, requests(n_requests), ierror & !communicator, request handle, error
                      )
        n_requests = n_requests + 1
        iavc = jtot(k, 2)
        call mpi_irecv( iexsolhir(iavc)%sol(1:irecexr(iavc)%muchineed(1), 1:itest), & !recvbuf
                        irecexr(iavc)%muchineed(1)*itest, mpi_double_precision, & !recvcount, recvtype
                        iexsolhir(iavc)%procid, 0, & !source, tag
                        mpi_comm_world, requests(n_requests), ierror & !communicator, request handle, error
                      )
      end if
    end do

    call mpi_waitall(n_requests, requests, mpi_statuses_ignore, ierror)
    deallocate(requests)
  end subroutine exchange_higher_av

  subroutine exchange_higher_pre(n)
    ! @brief
    ! this subroutine is establishing the communication pattern and allocating the appropriate memory for
    ! exchanging the variables of all the halo cells for all the reconstruction stencils
    implicit none
    integer, intent(in)::n
    integer::i, j, k, l, m, o, p, q, ineedt, tneedt, indl, tndl, icpuid, itest, itee, iteedum, itemp1, itemp2, iavc, iavt
    real, dimension(1:1)::dumts, rumts
    itest = nof_variables + turbulenceequations + passivescalar
    ineedt = irecexr(1)%tot
    tneedt = irecexs(1)%tot

    if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
      do i = 1, tneedt
        do k = 1, irecexs(i)%muchtheyneed(1)
          iexsolhis(i)%sol(k, 1:nof_variables) = u_c(irecexs(i)%localref(k))%val(1, 1:nof_variables)
                  iexsolhis(i)%sol(k,nof_variables+1:nof_variables+turbulenceequations+passivescalar)=u_ct(irecexs(i)%localref(k))%val(1,1:turbulenceequations+passivescalar)
        end do
      end do
    else
      do i = 1, tneedt
        do k = 1, irecexs(i)%muchtheyneed(1)
                  iexsolhis(i)%sol(k,1:nof_variables+turbulenceequations+passivescalar)=u_c(irecexs(i)%localref(k))%val(1,1:nof_variables+turbulenceequations+passivescalar)
        end do
      end do
    end if

    jtotal = 0; jtotal1 = 0; jtotal2 = 0; jtotal3 = 0
    icpuid = n
    dumts = zero
    do i = 0, isize - 1
      if (i .ne. n) then
        do j = 1, tneedt
          iavt = 10000
          if (irecexs(j)%procid .eq. i) then
            itemp1 = irecexs(j)%muchtheyneed(1)*itest
            iavt = j
            go to 7001
          end if
        end do
7001    continue
        do k = 1, ineedt
          iavc = 10000
          if (irecexr(k)%procid .eq. i) then
            iavc = k
            itemp2 = irecexr(k)%muchineed(1)*itest
            go to 8001
          end if
        end do
8001    continue
        if ((iavt .eq. 10000) .and. (iavc .ne. 10000)) then
          !message 1
          jtotal1 = jtotal1 + 1
          call mpi_sendrecv(dumts(1:1), 1, mpi_double_precision, i, icpuid, &
                            iexsolhir(iavc)%sol(1:irecexr(iavc)%muchineed(1), 1:itest), &
                            itemp2, mpi_double_precision, iexsolhir(iavc)%procid, iexsolhir(iavc)%procid, mpi_comm_world, status, ierror)
        end if
        if ((iavt .ne. 10000) .and. (iavc .eq. 10000)) then
          !message 1
          jtotal2 = jtotal2 + 1
          call mpi_sendrecv(iexsolhis(iavt)%sol(1:irecexs(iavt)%muchtheyneed(1), 1:itest), itemp1, &
                            mpi_double_precision, iexsolhis(iavt)%procid, icpuid, &
                            dumts(1:1), 1, mpi_double_precision, i, i, mpi_comm_world, status, ierror)
        end if
        if ((iavt .ne. 10000) .and. (iavc .ne. 10000)) then
          !message 1
          jtotal3 = jtotal3 + 1
          call mpi_sendrecv(iexsolhis(iavt)%sol(1:irecexs(iavt)%muchtheyneed(1), 1:itest), itemp1, &
                            mpi_double_precision, iexsolhis(iavt)%procid, icpuid, &
                            iexsolhir(iavc)%sol(1:irecexr(iavc)%muchineed(1), 1:itest), itemp2, mpi_double_precision, &
                            iexsolhir(iavc)%procid, iexsolhir(iavc)%procid, mpi_comm_world, status, ierror)
        end if
      end if
    end do
    allocate(jtot1(jtotal1, 3), jtot2(jtotal2, 3), jtot3(jtotal3, 3), jtot(jtotal1 + jtotal2 + jtotal3, 3))
    jtot1(:, :) = -1; jtot2(:, :) = -1; jtot3(:, :) = -1; jtot(:, :) = -1
    jtotal = 0; jtotal1 = 0; jtotal2 = 0; jtotal3 = 0
    icpuid = n
    dumts = zero
    do i = 0, isize - 1
      if (i .ne. n) then
        do j = 1, tneedt
          iavt = 10000
          if (irecexs(j)%procid .eq. i) then
            itemp1 = irecexs(j)%muchtheyneed(1)*itest
            iavt = j
            go to 10001
          end if
        end do
10001   continue
        do k = 1, ineedt
          iavc = 10000
          if (irecexr(k)%procid .eq. i) then
            iavc = k
            itemp2 = irecexr(k)%muchineed(1)*itest
            go to 11001
          end if
        end do
11001   continue
        if ((iavt .eq. 10000) .and. (iavc .ne. 10000)) then
          !message 1
          jtotal = jtotal + 1
          jtot(jtotal, 1) = -1
          jtot(jtotal, 2) = k
          jtot(jtotal, 3) = i
          call mpi_sendrecv(dumts(1:1), 1, mpi_double_precision, i, icpuid, &
                            iexsolhir(iavc)%sol(1:irecexr(iavc)%muchineed(1), 1:itest), &
                            itemp2, mpi_double_precision, iexsolhir(iavc)%procid, iexsolhir(iavc)%procid, mpi_comm_world, status, ierror)
        end if
        if ((iavt .ne. 10000) .and. (iavc .eq. 10000)) then
          !message 1
          jtotal = jtotal + 1
          jtot(jtotal, 1) = j
          jtot(jtotal, 2) = -1
          jtot(jtotal, 3) = i
          call mpi_sendrecv(iexsolhis(iavt)%sol(1:irecexs(iavt)%muchtheyneed(1), 1:itest), itemp1, &
                            mpi_double_precision, iexsolhis(iavt)%procid, icpuid, &
                            dumts(1:1), 1, mpi_double_precision, i, i, mpi_comm_world, status, ierror)
        end if
        if ((iavt .ne. 10000) .and. (iavc .ne. 10000)) then
          !message 1
          jtotal = jtotal + 1
          jtot(jtotal, 1) = j
          jtot(jtotal, 2) = k
          call mpi_sendrecv(iexsolhis(iavt)%sol(1:irecexs(iavt)%muchtheyneed(1), 1:itest), itemp1, &
                            mpi_double_precision, iexsolhis(iavt)%procid, icpuid, &
                            iexsolhir(iavc)%sol(1:irecexr(iavc)%muchineed(1), 1:itest), itemp2, mpi_double_precision, &
                            iexsolhir(iavc)%procid, iexsolhir(iavc)%procid, mpi_comm_world, status, ierror)
        end if
      end if
    end do
  end subroutine exchange_higher_pre

  subroutine exhboundhigher(n)
    ! @brief
    ! this subroutine is communicating the boundary extrapolated values for the variables and their gradients
    ! for the gaussian quadrature points of direct-side neighbours between mpi processes
    implicit none
    integer, intent(in)::n
    integer::i, j, k, l, m, o, p, q, ineedt, tneedt, indl, tndl, icpuid, ittt, iex, imulti, k_cnt, nvar
    integer::itee, iteedum, jk, jjk, jjk4, jjk12, imulti2, icpe, jmnb, j76, j78, j79, j80, imulti3, i_cnt, cinout2
    integer:: n_requests
    integer, dimension(:), allocatable:: requests
    real::pr_t31, pr_t32, pr_t33, pr_t34, pr_t35, temp_prin, temp_prout
    cinout2 = 0
    indl = iexchanger(1)%tot
    tndl = iexchanges(1)%tot
    pr_t31 = zero
    pr_t32 = zero
    pr_t33 = zero
    pr_t34 = zero
    pr_t35 = zero
    temp_prin = zero
    temp_prout = zero

    if (indl .ne. tndl) then
      write (*, *) "exhbounhigher: indl and tndl are supposed to be equal; indl=", indl, "tndl=", tndl
      call mpi_abort(mpi_comm_world, 1, ierror)
    end if

    if (dimensiona .eq. 3) then
      if (itestcase .eq. 4) then
        i_cnt = (nof_variables + turbulenceequations + passivescalar) + ((4 + turbulenceequations + passivescalar)*3)
      else
        i_cnt = nof_variables
      end if
    else
      if (itestcase .eq. 4) then
        i_cnt = (nof_variables + turbulenceequations + passivescalar) + ((3 + turbulenceequations + passivescalar)*2)
      else
        i_cnt = nof_variables
      end if
    end if

    if (itestcase .le. 3) then
      do i = 1, tndl
        do k = 1, iexchanges(i)%muchtheyneed(1)
          iexboundhis(i)%facesol(k,1:nof_variables)=ilocal_recon3(iexchanges(i)%localref(k))%uleft(1:nof_variables,iexchanges(i)%sidetheyneed(k),iexchanges(i)%qtheyneed(k))
        end do
      end do
    end if

    if (itestcase .eq. 4) then
      if (turbulence .ne. 1) then
        do i = 1, tndl
          do k = 1, iexchanges(i)%muchtheyneed(1)
            iexboundhis(i)%facesol(k,1:nof_variables)=ilocal_recon3(iexchanges(i)%localref(k))%uleft(1:nof_variables,iexchanges(i)%sidetheyneed(k),iexchanges(i)%qtheyneed(k))
            ittt = 0
            do iex = 1, nof_variables - 1
              do nvar = 1, dims
                ittt = ittt + 1
                iexboundhis(i)%facesol(k,nof_variables+ittt)=ilocal_recon3(iexchanges(i)%localref(k))%uleftv(nvar,iex,iexchanges(i)%sidetheyneed(k),iexchanges(i)%qtheyneed(k))
              end do
            end do
          end do
        end do
      else
        do i = 1, tndl
          do k = 1, iexchanges(i)%muchtheyneed(1)
            iexboundhis(i)%facesol(k,1:nof_variables)=ilocal_recon3(iexchanges(i)%localref(k))%uleft(1:nof_variables,iexchanges(i)%sidetheyneed(k),iexchanges(i)%qtheyneed(k))
            iexboundhis(i)%facesol(k,nof_variables+1:nof_variables+turbulenceequations+passivescalar)=ilocal_recon3(iexchanges(i)%localref(k))%uleftturb(1:turbulenceequations+passivescalar,iexchanges(i)%sidetheyneed(k),iexchanges(i)%qtheyneed(k))
            ittt = 0
            do iex = 1, nof_variables - 1
              do nvar = 1, dims
                ittt = ittt + 1
                iexboundhis(i)%facesol(k,nof_variables+turbulenceequations+passivescalar+ittt)=ilocal_recon3(iexchanges(i)%localref(k))%uleftv(nvar,iex,iexchanges(i)%sidetheyneed(k),iexchanges(i)%qtheyneed(k))
              end do
            end do
            do iex = 1, turbulenceequations + passivescalar
              do nvar = 1, dims
                ittt = ittt + 1
                iexboundhis(i)%facesol(k,nof_variables+turbulenceequations+passivescalar+ittt)=ilocal_recon3(iexchanges(i)%localref(k))%uleftturbv(nvar,iex,iexchanges(i)%sidetheyneed(k),iexchanges(i)%qtheyneed(k))
              end do
            end do
          end do
        end do
      end if
    end if

!-------------------for debugging only -----------------------------------------!

    n_requests = 0
    allocate(requests(2*indl))
    requests(:) = 0
    icpuid = n

    do k = 1, indl
      j = 1
      do while (iexboundhir(k)%procid .ne. iexboundhis(j)%procid)
        j = j + 1
      end do
      n_requests = n_requests + 1
      call mpi_isend( iexboundhis(j)%facesol(1:iexchanges(j)%muchtheyneed(1), 1:i_cnt), & !sendbuf
                      iexchanges(j)%muchtheyneed(1)*i_cnt, mpi_double_precision, & !sendcount, sendtype
                      iexboundhis(j)%procid, 0, & !destination, tag
                      mpi_comm_world, requests(n_requests), ierror & !communicator, request handle, error
                    )
      n_requests = n_requests + 1
      call mpi_irecv( iexboundhir(k)%facesol(1:iexchanger(k)%muchineed(1), 1:i_cnt), & !recvbuf
                      iexchanger(k)%muchineed(1)*i_cnt, mpi_double_precision, & !recvcount, recvtype
                      iexboundhir(k)%procid, 0, & !source, tag
                      mpi_comm_world, requests(n_requests), ierror & !communicator, request handle, error
                    )
    end do

    call mpi_waitall(n_requests, requests, mpi_statuses_ignore, ierror)
    deallocate(requests)
  end subroutine exhboundhigher

  subroutine exhboundhigher2(n)
    ! @brief
    ! this subroutine is communicating the boundary extrapolated values for the variables and their gradients
    ! for the gaussian quadrature points of direct-side neighbours between mpi processes for the implicit time stepping
    implicit none
    integer, intent(in)::n
    integer::i, j, k, l, m, o, p, q, ineedt, tneedt, indl, tndl, icpuid, ittt, iex, imulti, k_cnt, nvar
    integer::itee, iteedum, jk, jjk, jjk4, jjk12, imulti2, icpe, jmnb, j76, j78, j79, j80, imulti3
    indl = iexchanger(1)%tot
    tndl = iexchanges(1)%tot

    if (itestcase .lt. 3) then
      iex = 1
      imulti = iex
      imulti2 = iex
    end if
    if (itestcase .eq. 3) then
      iex = nof_variables
      imulti2 = iex
    end if
    if (itestcase .eq. 4) then
      k_cnt = (nof_variables + turbulenceequations + passivescalar)
      iex = (nof_variables + turbulenceequations + passivescalar)
      imulti2 = k_cnt
      imulti3 = k_cnt
    end if

    imulti = iex
    if (itestcase .le. 3) then
      do i = 1, tndl
        do k = 1, iexchanges(i)%muchtheyneed(1)
          if (relax .eq. 3) then
            if (iscoun .eq. 1) then
              iexboundhisi(i)%facesol(k, 1:iex) = -rhs(iexchanges(i)%localref(k))%val(1:iex)/impdiag_mf(iexchanges(i)%localref(k))
            else
              iexboundhisi(i)%facesol(k,1:iex)=-(rhs(iexchanges(i)%localref(k))%val(1:nof_variables)+((((1.5*u_c(iexchanges(i)%localref(k))%val(1,1:nof_variables))-(2.0d0*u_c(iexchanges(i)%localref(k))%val(2,1:nof_variables))+(0.5d0*u_c(iexchanges(i)%localref(k))%val(3,1:nof_variables)))/(dt))*ielem(n,iexchanges(i)%localref(k))%totvolume))/impdiag_mf(iexchanges(i)%localref(k))
            end if
          else
            iexboundhisi(i)%facesol(k, 1:iex) = impdu(iexchanges(i)%localref(k), 1:iex)
          end if
        end do
      end do
    end if
    if (itestcase .eq. 4) then
      do i = 1, tndl
        do k = 1, iexchanges(i)%muchtheyneed(1)
          if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
            do jjk = 1, nof_variables
              if (relax .eq. 3) then
                if (iscoun .eq. 1) then
                  iexboundhisi(i)%facesol(k, jjk) = -rhs(iexchanges(i)%localref(k))%val(jjk)/impdiag_mf(iexchanges(i)%localref(k))
                else
                  iexboundhisi(i)%facesol(k,jjk)=-(rhs(iexchanges(i)%localref(k))%val(jjk)+((((1.5*u_c(iexchanges(i)%localref(k))%val(1,jjk))-(2.0d0*u_c(iexchanges(i)%localref(k))%val(2,jjk))+(0.5d0*u_c(iexchanges(i)%localref(k))%val(3,jjk)))/(dt))*ielem(n,iexchanges(i)%localref(k))%totvolume))/impdiag_mf(iexchanges(i)%localref(k))
                end if
              else
                iexboundhisi(i)%facesol(k, jjk) = &
                impdu(iexchanges(i)%localref(k), jjk)
              end if
            end do
            do nvar = 1, 0 + turbulenceequations + passivescalar
              if (relax .eq. 3) then
                if (iscoun .eq. 1) then
                  iexboundhisi(i)%facesol(k,nof_variables+nvar)=-rhst(iexchanges(i)%localref(k))%val(nvar)/impdiagt(iexchanges(i)%localref(k),nvar)
                else
                  iexboundhisi(i)%facesol(k,nof_variables+nvar)=-(rhst(iexchanges(i)%localref(k))%val(nvar)+((((1.5*u_ct(iexchanges(i)%localref(k))%val(1,nvar))-(2.0d0*u_ct(iexchanges(i)%localref(k))%val(2,nvar))+(0.5d0*u_ct(iexchanges(i)%localref(k))%val(3,nvar)))/(dt))*ielem(n,iexchanges(i)%localref(k))%totvolume))/impdiagt(iexchanges(i)%localref(k),nvar)
                end if
              else
                iexboundhisi(i)%facesol(k, nof_variables + nvar) = &
                impdu(iexchanges(i)%localref(k), nof_variables + nvar)
              end if
            end do
          end if
          if ((turbulence .eq. 0) .and. (passivescalar .eq. 0)) then
            do jjk = 1, iex
              if (relax .eq. 3) then
                if (iscoun .eq. 1) then
                  iexboundhisi(i)%facesol(k, jjk) = -rhs(iexchanges(i)%localref(k))%val(jjk)/impdiag_mf(iexchanges(i)%localref(k))
                else
                  iexboundhisi(i)%facesol(k,jjk)=-(rhs(iexchanges(i)%localref(k))%val(jjk)+((((1.5*u_c(iexchanges(i)%localref(k))%val(1,jjk))-(2.0d0*u_c(iexchanges(i)%localref(k))%val(2,jjk))+(0.5d0*u_c(iexchanges(i)%localref(k))%val(3,jjk)))/(dt))*ielem(n,iexchanges(i)%localref(k))%totvolume))/impdiag_mf(iexchanges(i)%localref(k))
                end if
              else
                iexboundhisi(i)%facesol(k, jjk) = &
                impdu(iexchanges(i)%localref(k), jjk)
              end if
            end do
          end if ! turbulence
        end do
      end do
    end if

    icpuid = n
    if (itestcase .le. 3) then
      do k = 1, indl
        do j = 1, tndl
          if (iexboundhiri(k)%procid .eq. iexboundhisi(j)%procid) then
            call mpi_sendrecv(iexboundhisi(j)%facesol(1:iexchanges(j)%muchtheyneed(1), 1:iex) &
                              , iexchanges(j)%muchtheyneed(1)*imulti2, mpi_double_precision, iexboundhisi(j)%procid, &
                              icpuid, iexboundhiri(k)%facesol(1:iexchanger(k)%muchineed(1), 1:iex), &
                              iexchanger(k)%muchineed(1)*imulti2, mpi_double_precision, iexboundhiri(k)%procid, &
                              iexboundhiri(k)%procid, mpi_comm_world, status, ierror)
          end if
        end do
      end do
    end if

    if (itestcase .eq. 4) then
      do k = 1, indl
        do j = 1, tndl
          if (iexboundhiri(k)%procid .eq. iexboundhisi(j)%procid) then
            call mpi_sendrecv(iexboundhisi(j)%facesol(1:iexchanges(j)%muchtheyneed(1), 1:imulti3), &
                              iexchanges(j)%muchtheyneed(1)*imulti2, mpi_double_precision, iexboundhisi(j)%procid, &
                              icpuid, iexboundhiri(k)%facesol(1:iexchanger(k)%muchineed(1), 1:imulti3), &
                              iexchanger(k)%muchineed(1)*imulti2, mpi_double_precision, iexboundhiri(k)%procid, &
                              iexboundhiri(k)%procid, mpi_comm_world, status, ierror)
          end if
        end do
      end do
    end if
  end subroutine exhboundhigher2

  subroutine exhboundhigherlu(n)
    ! @brief
    ! this subroutine is communicating the boundary extrapolated values for the variables and their gradients
    ! for the gaussian quadrature points of direct-side neighbours between mpi processes for the implicit time stepping

    implicit none
    integer, intent(in)::n
    integer::i, j, k, l, m, o, p, q, ineedt, tneedt, indl, tndl, icpuid, ittt, iex, imulti, k_cnt, nvar
    integer::itee, iteedum, jk, jjk, jjk4, jjk12, imulti2, icpe, jmnb, j76, j78, j79, j80, imulti3
    indl = iexchanger(1)%tot
    tndl = iexchanges(1)%tot

    if (itestcase .lt. 3) then
      iex = 1
      imulti = iex
      imulti2 = iex
    end if
    if (itestcase .eq. 3) then
      iex = nof_variables
      imulti2 = iex
    end if
    if (itestcase .eq. 4) then
      k_cnt = (nof_variables + turbulenceequations + passivescalar)
      iex = (nof_variables + turbulenceequations + passivescalar)
      imulti2 = k_cnt
      imulti3 = k_cnt
    end if

    imulti = iex
    if (itestcase .le. 3) then
      do i = 1, tndl
        do k = 1, iexchanges(i)%muchtheyneed(1)
          iexboundhisi(i)%facesol(k, 1:iex) = impdu(iexchanges(i)%localref(k), 1:iex)
        end do
      end do
    end if

    if (itestcase .eq. 4) then
      do i = 1, tndl
        do k = 1, iexchanges(i)%muchtheyneed(1)
          if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
            do jjk = 1, nof_variables
              iexboundhisi(i)%facesol(k, jjk) = &
              impdu(iexchanges(i)%localref(k), jjk)
            end do
            do nvar = 1, 0 + turbulenceequations + passivescalar
              iexboundhisi(i)%facesol(k, nof_variables + nvar) = &
              impdu(iexchanges(i)%localref(k), nof_variables + nvar)
            end do
          end if
          if ((turbulence .eq. 0) .and. (passivescalar .eq. 0)) then
            do jjk = 1, iex
              iexboundhisi(i)%facesol(k, jjk) = &
              impdu(iexchanges(i)%localref(k), jjk)
            end do
          end if ! turbulence
        end do
      end do
    end if

    icpuid = n
    if (itestcase .le. 3) then
      do k = 1, indl
        do j = 1, tndl
          if (iexboundhiri(k)%procid .eq. iexboundhisi(j)%procid) then
            call mpi_sendrecv(iexboundhisi(j)%facesol(1:iexchanges(j)%muchtheyneed(1), 1:iex) &
                              , iexchanges(j)%muchtheyneed(1)*imulti2, mpi_double_precision, iexboundhisi(j)%procid, &
                              icpuid, iexboundhiri(k)%facesol(1:iexchanger(k)%muchineed(1), 1:iex), &
                              iexchanger(k)%muchineed(1)*imulti2, mpi_double_precision, iexboundhiri(k)%procid, &
                              iexboundhiri(k)%procid, mpi_comm_world, status, ierror)
          end if
        end do
      end do
    end if

    if (itestcase .eq. 4) then
      do k = 1, indl
        do j = 1, tndl
          if (iexboundhiri(k)%procid .eq. iexboundhisi(j)%procid) then
            call mpi_sendrecv(iexboundhisi(j)%facesol(1:iexchanges(j)%muchtheyneed(1), 1:imulti3), &
                              iexchanges(j)%muchtheyneed(1)*imulti2, mpi_double_precision, iexboundhisi(j)%procid, &
                              icpuid, iexboundhiri(k)%facesol(1:iexchanger(k)%muchineed(1), 1:imulti3), &
                              iexchanger(k)%muchineed(1)*imulti2, mpi_double_precision, iexboundhiri(k)%procid, &
                              iexboundhiri(k)%procid, mpi_comm_world, status, ierror)
          end if
        end do
      end do
    end if

  end subroutine exhboundhigherlu

  subroutine exhboundhigher_mood(n)
    implicit none
    integer, intent(in)::n
    integer::i, j, k, l, m, o, p, q, ineedt, tneedt, indl, tndl, icpuid, ittt, iex, imulti, k_cnt, nvar
    integer::itee, iteedum, jk, jjk, jjk4, jjk12, imulti2, icpe, jmnb, j76, j78, j79, j80, imulti3, i_cnt, cinout2
    integer:: n_requests
    integer, dimension(:), allocatable:: requests
    real::pr_t31, pr_t32, pr_t33, pr_t34, pr_t35, temp_prin, temp_prout

    cinout2 = 0
    indl = iexchanger(1)%tot
    tndl = iexchanges(1)%tot
    pr_t31 = zero
    pr_t32 = zero
    pr_t33 = zero
    pr_t34 = zero
    pr_t35 = zero
    temp_prin = zero
    temp_prout = zero

    if (indl .ne. tndl) then
      write (*, *) "exhbounhigher: indl and tndl are supposed to be equal; indl=", indl, "tndl=", tndl
      call mpi_abort(mpi_comm_world, 1, ierror)
    end if

    if (dimensiona .eq. 3) then
      if (itestcase .eq. 4) then
        i_cnt = 1
      else
        i_cnt = 1
      end if
    else
      if (itestcase .eq. 4) then
        i_cnt = 1
      else
        i_cnt = 1
      end if
    end if

    if (itestcase .le. 3) then
      do i = 1, tndl
        do k = 1, iexchanges(i)%muchtheyneed(1)
          iexboundhis(i)%facesol_m(k, 1) = ielem(n, (iexchanges(i)%localref(k)))%mood
        end do
      end do
    end if

    if (itestcase .eq. 4) then
      if (turbulence .ne. 1) then
        do i = 1, tndl
          do k = 1, iexchanges(i)%muchtheyneed(1)
            iexboundhis(i)%facesol(k,1:nof_variables)=ilocal_recon3(iexchanges(i)%localref(k))%uleft(1:nof_variables,iexchanges(i)%sidetheyneed(k),iexchanges(i)%qtheyneed(k))
            ittt = 0
            do iex = 1, nof_variables - 1
              do nvar = 1, dims
                ittt = ittt + 1
                iexboundhis(i)%facesol(k,nof_variables+ittt)=ilocal_recon3(iexchanges(i)%localref(k))%uleftv(nvar,iex,iexchanges(i)%sidetheyneed(k),iexchanges(i)%qtheyneed(k))
              end do
            end do
          end do
        end do
      else
        do i = 1, tndl
          do k = 1, iexchanges(i)%muchtheyneed(1)
            iexboundhis(i)%facesol(k,1:nof_variables)=ilocal_recon3(iexchanges(i)%localref(k))%uleft(1:nof_variables,iexchanges(i)%sidetheyneed(k),iexchanges(i)%qtheyneed(k))
            iexboundhis(i)%facesol(k,nof_variables+1:nof_variables+turbulenceequations+passivescalar)=ilocal_recon3(iexchanges(i)%localref(k))%uleftturb(1:turbulenceequations+passivescalar,iexchanges(i)%sidetheyneed(k),iexchanges(i)%qtheyneed(k))
            ittt = 0
            do iex = 1, nof_variables - 1
              do nvar = 1, dims
                ittt = ittt + 1
                iexboundhis(i)%facesol(k,nof_variables+turbulenceequations+passivescalar+ittt)=ilocal_recon3(iexchanges(i)%localref(k))%uleftv(nvar,iex,iexchanges(i)%sidetheyneed(k),iexchanges(i)%qtheyneed(k))
              end do
            end do
            do iex = 1, turbulenceequations + passivescalar
              do nvar = 1, dims
                ittt = ittt + 1
                iexboundhis(i)%facesol(k,nof_variables+turbulenceequations+passivescalar+ittt)=ilocal_recon3(iexchanges(i)%localref(k))%uleftturbv(nvar,iex,iexchanges(i)%sidetheyneed(k),iexchanges(i)%qtheyneed(k))
              end do
            end do
          end do
        end do
      end if
    end if

    n_requests = 0
    allocate(requests(2*indl))
    icpuid = n

    do k = 1, indl
      j = 1
      do while (iexboundhir(k)%procid .ne. iexboundhis(j)%procid)
        j = j + 1
      end do
      n_requests = n_requests + 1
      call mpi_isend( &
                      iexboundhis(j)%facesol_m(1:iexchanges(j)%muchtheyneed(1), 1:i_cnt), & !sendbuf
                      iexchanges(j)%muchtheyneed(1)*i_cnt, mpi_double_precision, & !sendcount, sendtype
                      iexboundhis(j)%procid, 0, & !destination, tag
                      mpi_comm_world, requests(n_requests), ierror & !communicator, request handle, error
                      )
      n_requests = n_requests + 1
      call mpi_irecv( &
                      iexboundhir(k)%facesol_m(1:iexchanger(k)%muchineed(1), 1:i_cnt), & !recvbuf
                      iexchanger(k)%muchineed(1)*i_cnt, mpi_double_precision, & !recvcount, recvtype
                      iexboundhir(k)%procid, 0, & !source, tag
                      mpi_comm_world, requests(n_requests), ierror & !communicator, request handle, error
                    )
    end do

    call mpi_waitall(n_requests, requests, mpi_statuses_ignore, ierror)
    deallocate(requests)
  end subroutine exhboundhigher_mood

  subroutine exhboundhigher_dg(n)
    ! @brief
    ! this subroutine is communicating the boundary extrapolated values for the variables and their gradients
    ! for the gaussian quadrature points of direct-side neighbours between mpi processes
    implicit none
    integer, intent(in)::n
    integer::i, j, k, l, m, o, p, q, ineedt, tneedt, indl, tndl, icpuid, ittt, iex, imulti, k_cnt, nvar
    integer::itee, iteedum, jk, jjk, jjk4, jjk12, imulti2, icpe, jmnb, j76, j78, j79, j80, imulti3, i_cnt, cinout2
    integer:: n_requests
    integer, dimension(:), allocatable:: requests
    real::pr_t31, pr_t32, pr_t33, pr_t34, pr_t35, temp_prin, temp_prout

    cinout2 = 0
    indl = iexchanger(1)%tot
    tndl = iexchanges(1)%tot
    pr_t31 = zero
    pr_t32 = zero
    pr_t33 = zero
    pr_t34 = zero
    pr_t35 = zero
    temp_prin = zero
    temp_prout = zero

    if (indl .ne. tndl) then
      write (*, *) "exhbounhigher: indl and tndl are supposed to be equal; indl=", indl, "tndl=", tndl
      call mpi_abort(mpi_comm_world, 1, ierror)
    end if

    i_cnt = nof_variables
    do i = 1, tndl
      do k = 1, iexchanges(i)%muchtheyneed(1)
        iexboundhis(i)%facesol_dg(k,1:nof_variables)=ilocal_recon3(iexchanges(i)%localref(k))%uleft_dg(1:nof_variables,iexchanges(i)%sidetheyneed(k),iexchanges(i)%qtheyneed(k))
      end do
    end do

!-------------------for debugging only -----------------------------------------!

    n_requests = 0
    allocate(requests(2*indl))
    requests(:) = 0
    icpuid = n

    do k = 1, indl
      j = 1
      do while (iexboundhir(k)%procid .ne. iexboundhis(j)%procid)
        j = j + 1
      end do
      n_requests = n_requests + 1
      call mpi_isend( iexboundhis(j)%facesol_dg(1:iexchanges(j)%muchtheyneed(1), 1:i_cnt), & !sendbuf
                      iexchanges(j)%muchtheyneed(1)*i_cnt, mpi_double_precision, & !sendcount, sendtype
                      iexboundhis(j)%procid, 0, & !destination, tag
                      mpi_comm_world, requests(n_requests), ierror & !communicator, request handle, error
                    )
      n_requests = n_requests + 1
      call mpi_irecv( iexboundhir(k)%facesol_dg(1:iexchanger(k)%muchineed(1), 1:i_cnt), & !recvbuf
                      iexchanger(k)%muchineed(1)*i_cnt, mpi_double_precision, & !recvcount, recvtype
                      iexboundhir(k)%procid, 0, & !source, tag
                      mpi_comm_world, requests(n_requests), ierror & !communicator, request handle, error
                    )
    end do
    call mpi_waitall(n_requests, requests, mpi_statuses_ignore, ierror)
    deallocate(requests)
  end subroutine exhboundhigher_dg

  subroutine exhboundhigher_dg2(n)
    ! @brief
    ! this subroutine is communicating the boundary extrapolated values for the variables and their gradients
    ! for the gaussian quadrature points of direct-side neighbours between mpi processes
    implicit none
    integer, intent(in)::n
    integer::i, j, k, l, m, o, p, q, ineedt, tneedt, indl, tndl, icpuid, ittt, iex, imulti, k_cnt, nvar
    integer::itee, iteedum, jk, jjk, jjk4, jjk12, imulti2, icpe, jmnb, j76, j78, j79, j80, imulti3, i_cnt, cinout2
    integer:: n_requests
    integer, dimension(:), allocatable:: requests
    real::pr_t31, pr_t32, pr_t33, pr_t34, pr_t35, temp_prin, temp_prout
    cinout2 = 0
    indl = iexchanger(1)%tot
    tndl = iexchanges(1)%tot
    pr_t31 = zero
    pr_t32 = zero
    pr_t33 = zero
    pr_t34 = zero
    pr_t35 = zero
    temp_prin = zero
    temp_prout = zero

    if (indl .ne. tndl) then
      write (*, *) "exhbounhigher: indl and tndl are supposed to be equal; indl=", indl, "tndl=", tndl
      call mpi_abort(mpi_comm_world, 1, ierror)
    end if

    if (dimensiona .eq. 3) then
      if (itestcase .eq. 4) then
        i_cnt = (nof_variables + turbulenceequations + passivescalar) + ((4 + turbulenceequations + passivescalar)*3)
      else
        i_cnt = nof_variables
      end if
    else
      if (itestcase .eq. 4) then
        i_cnt = (nof_variables + turbulenceequations + passivescalar) + ((3 + turbulenceequations + passivescalar)*2)
      else
        i_cnt = nof_variables
      end if
    end if

    if (itestcase .le. 3) then
      if ((governingequations .eq. -1) .and. (br2_yn .eq. 1)) then
        do i = 1, tndl
          do k = 1, iexchanges(i)%muchtheyneed(1)
            iexboundhis(i)%facesol_dg(k,1:nof_variables)=ilocal_recon3(iexchanges(i)%localref(k))%uleft_dg(1:nof_variables,iexchanges(i)%sidetheyneed(k),iexchanges(i)%qtheyneed(k))
            ittt = 0
            do iex = 1, nof_variables - 4
              do nvar = 1, dims
                ittt = ittt + 1
                iexboundhis(i)%facesol_dg(k,nof_variables+ittt)=ilocal_recon3(iexchanges(i)%localref(k))%br2_aux_var(iex,nvar,iexchanges(i)%sidetheyneed(k),iexchanges(i)%qtheyneed(k))
              end do
            end do
          end do
        end do
      else
        do i = 1, tndl
          do k = 1, iexchanges(i)%muchtheyneed(1)
                iexboundhis(i)%facesol_dg(k,1:nof_variables)=ilocal_recon3(iexchanges(i)%localref(k))%uleft_dg(1:nof_variables,iexchanges(i)%sidetheyneed(k),iexchanges(i)%qtheyneed(k))
          end do
        end do
      end if
    end if

    if (itestcase .eq. 4) then
      if (turbulence .ne. 1) then
        do i = 1, tndl
          do k = 1, iexchanges(i)%muchtheyneed(1)
            iexboundhis(i)%facesol_dg(k,1:nof_variables)=ilocal_recon3(iexchanges(i)%localref(k))%uleft_dg(1:nof_variables,iexchanges(i)%sidetheyneed(k),iexchanges(i)%qtheyneed(k))
            ittt = 0
            do iex = 1, nof_variables - 1
              do nvar = 1, dims
                ittt = ittt + 1
                iexboundhis(i)%facesol_dg(k,nof_variables+ittt)=ilocal_recon3(iexchanges(i)%localref(k))%br2_aux_var(iex+1,nvar,iexchanges(i)%sidetheyneed(k),iexchanges(i)%qtheyneed(k))
              end do
            end do
          end do
        end do
      else
        do i = 1, tndl
          do k = 1, iexchanges(i)%muchtheyneed(1)
            iexboundhis(i)%facesol_dg(k,1:nof_variables)=ilocal_recon3(iexchanges(i)%localref(k))%uleft_dg(1:nof_variables,iexchanges(i)%sidetheyneed(k),iexchanges(i)%qtheyneed(k))
            iexboundhis(i)%facesol(k,nof_variables+1:nof_variables+turbulenceequations+passivescalar)=ilocal_recon3(iexchanges(i)%localref(k))%uleftturb(1:turbulenceequations+passivescalar,iexchanges(i)%sidetheyneed(k),iexchanges(i)%qtheyneed(k))
            ittt = 0
            do iex = 1, nof_variables - 1
              do nvar = 1, dims
                ittt = ittt + 1
                iexboundhis(i)%facesol(k,nof_variables+turbulenceequations+passivescalar+ittt)=ilocal_recon3(iexchanges(i)%localref(k))%uleftv(nvar,iex,iexchanges(i)%sidetheyneed(k),iexchanges(i)%qtheyneed(k))
              end do
            end do
            do iex = 1, turbulenceequations + passivescalar
              do nvar = 1, dims
                ittt = ittt + 1
                iexboundhis(i)%facesol(k,nof_variables+turbulenceequations+passivescalar+ittt)=ilocal_recon3(iexchanges(i)%localref(k))%uleftturbv(nvar,iex,iexchanges(i)%sidetheyneed(k),iexchanges(i)%qtheyneed(k))
              end do
            end do
          end do
        end do
      end if
    end if

!-------------------for debugging only -----------------------------------------!

    n_requests = 0
    allocate(requests(2*indl))
    requests(:) = 0
    icpuid = n

    do k = 1, indl
      j = 1
      do while (iexboundhir(k)%procid .ne. iexboundhis(j)%procid)
        j = j + 1
      end do
      n_requests = n_requests + 1
      call mpi_isend( &
                      iexboundhis(j)%facesol_dg(1:iexchanges(j)%muchtheyneed(1), 1:i_cnt), & !sendbuf
                      iexchanges(j)%muchtheyneed(1)*i_cnt, mpi_double_precision, & !sendcount, sendtype
                      iexboundhis(j)%procid, 0, & !destination, tag
                      mpi_comm_world, requests(n_requests), ierror & !communicator, request handle, error
                    )
      n_requests = n_requests + 1
      call mpi_irecv( &
                      iexboundhir(k)%facesol_dg(1:iexchanger(k)%muchineed(1), 1:i_cnt), & !recvbuf
                      iexchanger(k)%muchineed(1)*i_cnt, mpi_double_precision, & !recvcount, recvtype
                      iexboundhir(k)%procid, 0, & !source, tag
                      mpi_comm_world, requests(n_requests), ierror & !communicator, request handle, error
                    )
    end do
    call mpi_waitall(n_requests, requests, mpi_statuses_ignore, ierror)
    deallocate(requests)
  end subroutine exhboundhigher_dg2

end module communications
