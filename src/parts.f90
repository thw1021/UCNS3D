module partition
  use mpiinfo
  use declaration
  implicit none
contains
  subroutine partitioner5(n, imaxe, imaxn, xmpie, ieshape)
    !> @brief
    !> this subroutine partitions the mesh using metis
    integer, intent(in) :: n, imaxe, imaxn
    integer, allocatable, dimension(:), intent(in) :: ieshape
    integer, allocatable, dimension(:), intent(inout) :: xmpie
    integer :: i, j, k, ios, iox, ioz, i1, i2, i3, i4, i5, i6, i7, i8, i9, i10, i11, ioy, ihe, itri, idu
    character(len=12) :: metfile, celfile
    integer :: posa
    real :: xc, yc
    character(len=10) :: t, f, ss

    celfile = 'grid.cel'
    metfile = 'grid'
    open (9, file=metfile, form='formatted', status='new', action='write')
    open (8, file=celfile, form='formatted', status='old', action='read', iostat=ios)
    write (9, *) imaxe
    ihe = 0
    itri = 0

    if (dimensiona .eq. 3) then
      do j = 1, imaxe
        read (8, *) k, i1, i2, i3, i4, i5, i6, i7, i8
        write (9, "(8i10)") i1, i2, i3, i4, i5, i6, i7, i8
      end do
    else
      do j = 1, imaxe
        read (8, *) k, i1, i2, i3, i4
        if ((i3 .ne. i4)) then
          write (9, "(4i10)") i1, i2, i3, i4
        else
          write (9, "(3i10)") i1, i2, i3
        end if
      end do
    end if

    close (8)
    close (9)

    if (k .ge. 0) then
      do idu = 1, itold
        yc = xc + yc**2
        xc = xc*idu**2 + sqrt(xc**2)*sqrt(xc**4)*cos(xc)*sin(yc)*(1.0/(xc*sin(xc**4)))
        xc = yc*idu**2 + sqrt(xc**2)*sqrt(xc**4)*cos(xc)*sin(yc)*(1.0/(xc*sin(xc**4)))
        yc = xc*idu**2 + sqrt(xc**2)*sqrt(xc**4)*cos(xc)*sin(yc)*(1.0/(xc*sin(xc**4)))
      end do
    end if

    posa = isize
    write (t, fmt='(i10)') posa
    f = trim(adjustl(t))
    ss = trim(adjustr(f))
    call system('pwd')

    if (k .ge. 0) then
      do idu = 1, itold
        yc = xc + yc**2
        xc = xc*idu**2 + sqrt(xc**2)*sqrt(xc**4)*cos(xc)*sin(yc)*(1.0/(xc*sin(xc**4)))
        xc = yc*idu**2 + sqrt(xc**2)*sqrt(xc**4)*cos(xc)*sin(yc)*(1.0/(xc*sin(xc**4)))
        yc = xc*idu**2 + sqrt(xc**2)*sqrt(xc**4)*cos(xc)*sin(yc)*(1.0/(xc*sin(xc**4)))
      end do
    end if

    call system('./mpmetis star'//ss)

    if (k .ge. 0) then
      do idu = 1, itold
        yc = xc + yc**2
        xc = xc*idu**2 + sqrt(xc**2)*sqrt(xc**4)*cos(xc)*sin(yc)*(1.0/(xc*sin(xc**4)))
        xc = yc*idu**2 + sqrt(xc**2)*sqrt(xc**4)*cos(xc)*sin(yc)*(1.0/(xc*sin(xc**4)))
        yc = xc*idu**2 + sqrt(xc**2)*sqrt(xc**4)*cos(xc)*sin(yc)*(1.0/(xc*sin(xc**4)))
      end do
    end if

    call system('mv grid.epart.* grid.epart')

    if (k .ge. 0) then
      do idu = 1, itold
        yc = xc + yc**2
        xc = xc*idu**2 + sqrt(xc**2)*sqrt(xc**4)*cos(xc)*sin(yc)*(1.0/(xc*sin(xc**4)))
        xc = yc*idu**2 + sqrt(xc**2)*sqrt(xc**4)*cos(xc)*sin(yc)*(1.0/(xc*sin(xc**4)))
        yc = xc*idu**2 + sqrt(xc**2)*sqrt(xc**4)*cos(xc)*sin(yc)*(1.0/(xc*sin(xc**4)))
      end do
    end if

    open (10, file='grid.epart', form='formatted', status='old', action='read')
    do i = 1, imaxe
      read (10, *) k
      xmpie(i) = k
    end do
    close (10)

    if (k .ge. 0) then
      do idu = 1, itold
        yc = xc + yc**2
        xc = xc*idu**2 + sqrt(xc**2)*sqrt(xc**4)*cos(xc)*sin(yc)*(1.0/(xc*sin(xc**4)))
        xc = yc*idu**2 + sqrt(xc**2)*sqrt(xc**4)*cos(xc)*sin(yc)*(1.0/(xc*sin(xc**4)))
        yc = xc*idu**2 + sqrt(xc**2)*sqrt(xc**4)*cos(xc)*sin(yc)*(1.0/(xc*sin(xc**4)))
      end do
    end if

    call system('rm -rf grid.*part* grid')
  end subroutine partitioner5

  subroutine partitioner1(n, imaxe, imaxn, xmpie, ieshape)
    !> @brief
    !> this subroutine partitions the mesh using metis
    use iso_c_binding
    implicit none
    external metis_partmeshdual
    external metis_partmeshnodal
    integer, allocatable, dimension(:), intent(in) :: ieshape
    type :: elementglobal
      integer :: elementglid, nodeid1, nodeid2, nodeid3, nodeid4, nodeid5, nodeid6, nodeid7, nodeid8
    end type elementglobal
    real :: average
    integer :: maxi
    integer :: i, j, elementid, node1, node2, node3, node4, node5, node6, node7, node8, counternodes, k
    integer :: vweight
    integer, intent(in) :: n, imaxe, imaxn
    integer, allocatable, dimension(:), intent(inout) :: xmpie
    integer, allocatable, dimension(:) :: testar, pweight
    integer(c_int) :: imaxee, imaxnn, ncommon1, isizee, objval
    integer(c_int), allocatable, dimension(:) :: allnodesptr, allnodes
    integer(c_int), allocatable, dimension(:) :: xmpiee, xmpidumb, vwgt, vsize
    real(c_float), allocatable, dimension(:) :: tpwgts
    integer(c_int), dimension(0:39) :: options
    type(elementglobal), allocatable, dimension(:) :: elements

    imaxee = imaxe
    imaxnn = imaxn
    isizee = isize
    allocate(testar(imaxe))
    call metis_setdefaultoptions(options)
    counternodes = 0
    allocate(elements(imaxee))
    allocate(allnodes(allnodesgloball))
    allocate(allnodesptr(imaxe + 1))

    open(1112, file='grid.cel', form='formatted', status='old', action='read')
    allnodesptr(1) = 0
    do i = 1, imaxe
      read(1112, *) elementid, node1, node2, node3, node4, node5, node6, node7, node8
      elements(i)%elementglid = elementid
      elements(i)%nodeid1 = node1
      elements(i)%nodeid2 = node2
      elements(i)%nodeid3 = node3
      elements(i)%nodeid4 = node4
      elements(i)%nodeid5 = node5
      elements(i)%nodeid6 = node6
      elements(i)%nodeid7 = node7
      elements(i)%nodeid8 = node8

      if ((node3 .eq. node4) .and. (node5 .ne. node6) .and. (node7 .eq. node8)) then ! prism
        counternodes = counternodes + 1
        allnodes(counternodes) = node1
        counternodes = counternodes + 1
        allnodes(counternodes) = node2
        counternodes = counternodes + 1
        allnodes(counternodes) = node3
        counternodes = counternodes + 1
        allnodes(counternodes) = node5
        counternodes = counternodes + 1
        allnodes(counternodes) = node6
        counternodes = counternodes + 1
        allnodes(counternodes) = node7
      else if ((node5 .ne. node6) .and. (node6 .ne. node7) .and. (node7 .ne. node8) .and. (node3 .ne. node4)) then ! hexa
        counternodes = counternodes + 1
        allnodes(counternodes) = node1
        counternodes = counternodes + 1
        allnodes(counternodes) = node2
        counternodes = counternodes + 1
        allnodes(counternodes) = node3
        counternodes = counternodes + 1
        allnodes(counternodes) = node4
        counternodes = counternodes + 1
        allnodes(counternodes) = node5
        counternodes = counternodes + 1
        allnodes(counternodes) = node6
        counternodes = counternodes + 1
        allnodes(counternodes) = node7
        counternodes = counternodes + 1
        allnodes(counternodes) = node8
      else if ((node3 .eq. node4) .and. (node5 .eq. node6) .and. (node6 .eq. node7) .and. (node7 .eq. node8)) then ! tetra
        counternodes = counternodes + 1
        allnodes(counternodes) = node1
        counternodes = counternodes + 1
        allnodes(counternodes) = node2
        counternodes = counternodes + 1
        allnodes(counternodes) = node3
        counternodes = counternodes + 1
        allnodes(counternodes) = node5
      else
        counternodes = counternodes + 1
        allnodes(counternodes) = node1
        counternodes = counternodes + 1
        allnodes(counternodes) = node2
        counternodes = counternodes + 1
        allnodes(counternodes) = node3
        counternodes = counternodes + 1
        allnodes(counternodes) = node4
        counternodes = counternodes + 1
        allnodes(counternodes) = node5
      end if

      allnodesptr(i + 1) = counternodes
    end do
    close(1112)

    ncommon1 = 1
    allocate(vwgt(imaxee))
    allocate(vsize(imaxee))
    allocate(tpwgts(isizee))
    allocate(xmpidumb(imaxnn))
    allocate(xmpiee(imaxee))
    tpwgts(:) = 1.0/(1.0*isizee)
    vwgt(:) = 1
    vsize(:) = 1

    call metis_partmeshdual(imaxee, imaxnn, allnodesptr, allnodes, vwgt, vsize, ncommon1, isizee, tpwgts, options, objval, xmpiee, xmpidumb)

    do i = 1, imaxe
      xmpie(i) = xmpiee(i)
    end do

    allocate(pweight(0:isize - 1))
    pweight(:) = 0
    do i = 1, imaxe
      pweight(xmpiee(i)) = pweight(xmpiee(i)) + 1
    end do

    average = pweight(0)
    do i = 1, isize - 1
      average = pweight(i) + average
    end do
    average = average/isize
    maxi = maxval(pweight)

    open(63, file='history.txt', form='formatted', action='write', position='append')
    write(63, *) 'load imbalance of partioner', maxi/average
    close(63)

    deallocate(testar)
    deallocate(elements)
    deallocate(allnodes)
    deallocate(allnodesptr)
    deallocate(vwgt)
    deallocate(vsize)
    deallocate(tpwgts)
    deallocate(pweight)
    deallocate(xmpidumb)
    deallocate(xmpiee)
  end subroutine partitioner1

  subroutine partitioner2(n, imaxe, imaxn, xmpie, ieshape)
    !> @brief
    !> this subroutine partitions the mesh using metis
    use iso_c_binding
    implicit none
    external metis_partmeshdual
    external metis_partmeshnodal
    integer, allocatable, dimension(:), intent(in) :: ieshape
    type :: elementglobal
      integer :: elementglid, nodeid1, nodeid2, nodeid3, nodeid4, nodeid5, nodeid6, nodeid7, nodeid8
    end type elementglobal
    real :: average, average2
    integer :: maxi, maxi2
    integer :: i, j, elementid, node1, node2, node3, node4, node5, node6, node7, node8, counternodes, k
    integer :: vweight
    integer, allocatable, dimension(:, :, :) :: compweight, commweight
    integer, intent(in) :: n, imaxe, imaxn
    integer, allocatable, dimension(:), intent(inout) :: xmpie
    integer, allocatable, dimension(:) :: testar, pweight, pweight2
    integer(c_int) :: imaxee, imaxnn, ncommon1, isizee, objval
    integer(c_int), allocatable, dimension(:) :: allnodesptr, allnodes
    integer(c_int), allocatable, dimension(:) :: xmpiee, xmpidumb, vwgt, vsize
    real(c_float), allocatable, dimension(:) :: tpwgts
    integer(c_int), dimension(0:39) :: options
    type(elementglobal), allocatable, dimension(:) :: elements

    imaxee = imaxe
    imaxnn = imaxn
    isizee = isize
    allocate(testar(imaxe))
    call metis_setdefaultoptions(options)
    options(1) = 1 ! cut=0 or volume=1 objective
    counternodes = 0
    allocate(elements(imaxee))
    allocate(allnodes(allnodesgloball))
    allocate(allnodesptr(imaxe + 1))

    open(1112, file='grid.cel', form='formatted', status='old', action='read')
    allnodesptr(1) = 0
    do i = 1, imaxe
      read(1112, *) elementid, node1, node2, node3, node4, node5, node6, node7, node8
      elements(i)%elementglid = elementid
      elements(i)%nodeid1 = node1
      elements(i)%nodeid2 = node2
      elements(i)%nodeid3 = node3
      elements(i)%nodeid4 = node4
      elements(i)%nodeid5 = node5
      elements(i)%nodeid6 = node6
      elements(i)%nodeid7 = node7
      elements(i)%nodeid8 = node8

      if ((node3 .eq. node4) .and. (node5 .ne. node6) .and. (node7 .eq. node8)) then ! prism
        counternodes = counternodes + 1
        allnodes(counternodes) = node1
        counternodes = counternodes + 1
        allnodes(counternodes) = node2
        counternodes = counternodes + 1
        allnodes(counternodes) = node3
        counternodes = counternodes + 1
        allnodes(counternodes) = node5
        counternodes = counternodes + 1
        allnodes(counternodes) = node6
        counternodes = counternodes + 1
        allnodes(counternodes) = node7
      else if ((node5 .ne. node6) .and. (node6 .ne. node7) .and. (node7 .ne. node8) .and. (node3 .ne. node4)) then ! hexa
        counternodes = counternodes + 1
        allnodes(counternodes) = node1
        counternodes = counternodes + 1
        allnodes(counternodes) = node2
        counternodes = counternodes + 1
        allnodes(counternodes) = node3
        counternodes = counternodes + 1
        allnodes(counternodes) = node4
        counternodes = counternodes + 1
        allnodes(counternodes) = node5
        counternodes = counternodes + 1
        allnodes(counternodes) = node6
        counternodes = counternodes + 1
        allnodes(counternodes) = node7
        counternodes = counternodes + 1
        allnodes(counternodes) = node8
      else if ((node3 .eq. node4) .and. (node5 .eq. node6) .and. (node6 .eq. node7) .and. (node7 .eq. node8)) then ! tetra
        counternodes = counternodes + 1
        allnodes(counternodes) = node1
        counternodes = counternodes + 1
        allnodes(counternodes) = node2
        counternodes = counternodes + 1
        allnodes(counternodes) = node3
        counternodes = counternodes + 1
        allnodes(counternodes) = node5
      else
        counternodes = counternodes + 1
        allnodes(counternodes) = node1
        counternodes = counternodes + 1
        allnodes(counternodes) = node2
        counternodes = counternodes + 1
        allnodes(counternodes) = node3
        counternodes = counternodes + 1
        allnodes(counternodes) = node4
        counternodes = counternodes + 1
        allnodes(counternodes) = node5
      end if

      allnodesptr(i + 1) = counternodes
    end do
    close(1112)

    ncommon1 = 1
    allocate(vwgt(imaxee))
    allocate(vsize(imaxee))
    allocate(tpwgts(isizee))
    allocate(xmpidumb(imaxnn))
    allocate(xmpiee(imaxee))
    allocate(compweight(1:3, 1:10, 1:4))
    allocate(commweight(1:3, 1:10, 1:4))
    tpwgts(:) = 1.0/(1.0*isizee)
    compweight(:, :, :) = 1
    commweight(:, :, :) = 1

    compweight(:, :, 1) = 6
    compweight(:, :, 2) = 4
    compweight(:, :, 3) = 5
    compweight(:, :, 4) = 5

    compweight(3, :, 1) = 7
    compweight(3, :, 2) = 5
    compweight(3, :, 3) = 6
    compweight(3, :, 4) = 6

    commweight(:, :, 1) = 6
    commweight(:, :, 2) = 4
    commweight(:, :, 3) = 5
    commweight(:, :, 4) = 5

    do i = 1, imaxee
      vwgt(i) = compweight(spatiladiscret, spatialorder, ieshape(i))
    end do

    do i = 1, imaxee
      vsize(i) = commweight(1, 1, ieshape(i))
    end do

    deallocate(commweight)
    deallocate(compweight)

    call metis_partmeshdual(imaxee, imaxnn, allnodesptr, allnodes, vwgt, vsize, ncommon1, isizee, tpwgts, options, objval, xmpiee, xmpidumb)

    do i = 1, imaxe
      xmpie(i) = xmpiee(i)
    end do

    allocate(pweight(0:isize - 1))
    allocate(pweight2(0:isize - 1))
    pweight(:) = 0
    pweight2(:) = 0

    do i = 1, imaxe
      pweight(xmpiee(i)) = pweight(xmpiee(i)) + 1
      pweight2(xmpiee(i)) = pweight2(xmpiee(i)) + vwgt(i)
    end do

    average = pweight(0)
    average2 = pweight2(0)
    do i = 1, isize - 1
      average = pweight(i) + average
      average2 = pweight2(i) + average2
    end do
    average = average/isize
    average2 = average2/isize
    maxi = maxval(pweight)
    maxi2 = maxval(pweight2)

    open(63, file='history.txt', form='formatted', action='write', position='append')
    write(63, *) 'load imbalance of partioner', maxi/average
    write(63, *) 'load imbalance of partioner based on element weights', maxi2/average2
    close(63)

    deallocate(testar)
    deallocate(elements)
    deallocate(allnodes)
    deallocate(allnodesptr)
    deallocate(vwgt)
    deallocate(vsize)
    deallocate(tpwgts)
    deallocate(pweight)
    deallocate(pweight2)
    deallocate(xmpidumb)
    deallocate(xmpiee)
  end subroutine

  subroutine partitioner4(n, imaxe, imaxn, xmpie, ieshape)
    !> @brief
    !> this subroutine partitions the mesh using metis
    use iso_c_binding
    implicit none
    external metis_partmeshdual
    external metis_partmeshnodal
    integer, allocatable, dimension(:), intent(in) :: ieshape
    type :: elementglobal
      integer :: elementglid, nodeid1, nodeid2, nodeid3, nodeid4, nodeid5, nodeid6, nodeid7, nodeid8
    end type elementglobal
    real :: average, average2
    integer :: maxi, maxi2
    integer :: i, j, elementid, node1, node2, node3, node4, node5, node6, node7, node8, counternodes, k
    integer :: vweight
    integer, allocatable, dimension(:, :, :) :: compweight, commweight
    integer, intent(in) :: n, imaxe, imaxn
    integer, allocatable, dimension(:), intent(inout) :: xmpie
    integer, allocatable, dimension(:) :: testar, pweight, pweight2
    integer(c_int) :: imaxee, imaxnn, ncommon1, isizee, objval
    integer(c_int), allocatable, dimension(:) :: allnodesptr, allnodes
    integer(c_int), allocatable, dimension(:) :: xmpiee, xmpidumb, vwgt, vsize
    real(c_float), allocatable, dimension(:) :: tpwgts
    integer(c_int), dimension(0:39) :: options
    type(elementglobal), allocatable, dimension(:) :: elements

    imaxee = imaxe
    imaxnn = imaxn
    isizee = isize
    allocate(testar(imaxe))
    call metis_setdefaultoptions(options)
    options(1) = 1 ! cut=0 or volume=1 objective
    counternodes = 0
    allocate(elements(imaxee))
    allocate(allnodes(allnodesgloball))
    allocate(allnodesptr(imaxe + 1))

    open(1112, file='grid.cel', form='formatted', status='old', action='read')
    allnodesptr(1) = 0
    do i = 1, imaxe
      read(1112, *) elementid, node1, node2, node3, node4, node5, node6, node7, node8
      elements(i)%elementglid = elementid
      elements(i)%nodeid1 = node1
      elements(i)%nodeid2 = node2
      elements(i)%nodeid3 = node3
      elements(i)%nodeid4 = node4
      elements(i)%nodeid5 = node5
      elements(i)%nodeid6 = node6
      elements(i)%nodeid7 = node7
      elements(i)%nodeid8 = node8

      if ((node3 .eq. node4) .and. (node5 .ne. node6) .and. (node7 .eq. node8)) then ! prism
        counternodes = counternodes + 1
        allnodes(counternodes) = node1
        counternodes = counternodes + 1
        allnodes(counternodes) = node2
        counternodes = counternodes + 1
        allnodes(counternodes) = node3
        counternodes = counternodes + 1
        allnodes(counternodes) = node5
        counternodes = counternodes + 1
        allnodes(counternodes) = node6
        counternodes = counternodes + 1
        allnodes(counternodes) = node7
      else if ((node5 .ne. node6) .and. (node6 .ne. node7) .and. (node7 .ne. node8) .and. (node3 .ne. node4)) then ! hexa
        counternodes = counternodes + 1
        allnodes(counternodes) = node1
        counternodes = counternodes + 1
        allnodes(counternodes) = node2
        counternodes = counternodes + 1
        allnodes(counternodes) = node3
        counternodes = counternodes + 1
        allnodes(counternodes) = node4
        counternodes = counternodes + 1
        allnodes(counternodes) = node5
        counternodes = counternodes + 1
        allnodes(counternodes) = node6
        counternodes = counternodes + 1
        allnodes(counternodes) = node7
        counternodes = counternodes + 1
        allnodes(counternodes) = node8
      else if ((node3 .eq. node4) .and. (node5 .eq. node6) .and. (node6 .eq. node7) .and. (node7 .eq. node8)) then ! tetra
        counternodes = counternodes + 1
        allnodes(counternodes) = node1
        counternodes = counternodes + 1
        allnodes(counternodes) = node2
        counternodes = counternodes + 1
        allnodes(counternodes) = node3
        counternodes = counternodes + 1
        allnodes(counternodes) = node5
      else
        counternodes = counternodes + 1
        allnodes(counternodes) = node1
        counternodes = counternodes + 1
        allnodes(counternodes) = node2
        counternodes = counternodes + 1
        allnodes(counternodes) = node3
        counternodes = counternodes + 1
        allnodes(counternodes) = node4
        counternodes = counternodes + 1
        allnodes(counternodes) = node5
      end if

      allnodesptr(i + 1) = counternodes
    end do
    close(1112)

    ncommon1 = 1
    allocate(vwgt(imaxee))
    allocate(vsize(imaxee))
    allocate(tpwgts(isizee))
    allocate(xmpidumb(imaxnn))
    allocate(xmpiee(imaxee))
    allocate(compweight(1:3, 1:10, 1:4))
    allocate(commweight(1:3, 1:10, 1:4))
    tpwgts(:) = 1.0/(1.0*isizee)
    compweight(:, :, :) = 1
    commweight(:, :, :) = 1

    compweight(:, :, 1) = 6
    compweight(:, :, 2) = 4
    compweight(:, :, 3) = 5
    compweight(:, :, 4) = 5

    compweight(3, :, 1) = 7
    compweight(3, :, 2) = 5
    compweight(3, :, 3) = 6
    compweight(3, :, 4) = 6

    commweight(:, :, 1) = 6
    commweight(:, :, 2) = 4
    commweight(:, :, 3) = 5
    commweight(:, :, 4) = 5

    do i = 1, imaxee
      vwgt(i) = compweight(spatiladiscret, spatialorder, ieshape(i))
    end do

    do i = 1, imaxee
      vsize(i) = commweight(1, 1, ieshape(i))
    end do

    deallocate(commweight)
    deallocate(compweight)

    call metis_partmeshnodal(imaxee, imaxnn, allnodesptr, allnodes, vwgt, vsize, isizee, tpwgts, options, objval, xmpiee, xmpidumb)

    do i = 1, imaxe
      xmpie(i) = xmpiee(i)
    end do

    allocate(pweight(0:isize - 1))
    allocate(pweight2(0:isize - 1))
    pweight(:) = 0
    pweight2(:) = 0

    do i = 1, imaxe
      pweight(xmpiee(i)) = pweight(xmpiee(i)) + 1
      pweight2(xmpiee(i)) = pweight2(xmpiee(i)) + vwgt(i)
    end do

    average = pweight(0)
    average2 = pweight2(0)
    do i = 1, isize - 1
      average = pweight(i) + average
      average2 = pweight2(i) + average2
    end do
    average = average/isize
    average2 = average2/isize
    maxi = maxval(pweight)
    maxi2 = maxval(pweight2)

    open(63, file='history.txt', form='formatted', action='write', position='append')
    write(63, *) 'load imbalance of partioner', maxi/average
    write(63, *) 'load imbalance of partioner based on element weights', maxi2/average2
    close(63)

    deallocate(testar)
    deallocate(elements)
    deallocate(allnodes)
    deallocate(allnodesptr)
    deallocate(vwgt)
    deallocate(vsize)
    deallocate(tpwgts)
    deallocate(pweight)
    deallocate(pweight2)
    deallocate(xmpidumb)
    deallocate(xmpiee)
  end subroutine

  subroutine partitioner3(n, imaxe, imaxn, xmpie, ieshape)
    !> @brief
    !> this subroutine partitions the mesh using metis
    use iso_c_binding
    implicit none
    external metis_partmeshdual
    external metis_partmeshnodal
    integer, allocatable, dimension(:), intent(in) :: ieshape
    type :: elementglobal
      integer :: elementglid, nodeid1, nodeid2, nodeid3, nodeid4, nodeid5, nodeid6, nodeid7, nodeid8
    end type elementglobal
    real :: average
    integer :: maxi
    integer :: i, j, elementid, node1, node2, node3, node4, node5, node6, node7, node8, counternodes, k
    integer :: vweight
    integer, intent(in) :: n, imaxe, imaxn
    integer, allocatable, dimension(:), intent(inout) :: xmpie
    integer, allocatable, dimension(:) :: testar, pweight
    integer(c_int) :: imaxee, imaxnn, ncommon1, isizee, objval
    integer(c_int), allocatable, dimension(:) :: allnodesptr, allnodes
    integer(c_int), allocatable, dimension(:) :: xmpiee, xmpidumb, vwgt, vsize
    real(c_float), allocatable, dimension(:) :: tpwgts
    integer(c_int), dimension(0:39) :: options
    type(elementglobal), allocatable, dimension(:) :: elements

    imaxee = imaxe
    imaxnn = imaxn
    isizee = isize
    allocate(testar(imaxe))
    call metis_setdefaultoptions(options)
    counternodes = 0
    allocate(elements(imaxee))
    allnodesgloball = 8*imaxe
    allocate(allnodes(allnodesgloball))
    allocate(allnodesptr(imaxe + 1))

    open(1112, file='grid.cel', form='formatted', status='old', action='read')
    allnodesptr(1) = 0
    do i = 1, imaxe
      read(1112, *) elementid, node1, node2, node3, node4, node5, node6, node7, node8
      counternodes = counternodes + 1
      allnodes(counternodes) = node1
      counternodes = counternodes + 1
      allnodes(counternodes) = node2
      counternodes = counternodes + 1
      allnodes(counternodes) = node3
      counternodes = counternodes + 1
      allnodes(counternodes) = node4
      counternodes = counternodes + 1
      allnodes(counternodes) = node5
      counternodes = counternodes + 1
      allnodes(counternodes) = node6
      counternodes = counternodes + 1
      allnodes(counternodes) = node7
      counternodes = counternodes + 1
      allnodes(counternodes) = node8
      allnodesptr(i + 1) = counternodes
    end do
    close(1112)

    ncommon1 = 1
    allocate(vwgt(imaxee))
    allocate(vsize(imaxee))
    allocate(tpwgts(isizee))
    allocate(xmpidumb(imaxnn))
    allocate(xmpiee(imaxee))
    tpwgts(:) = 1.0/(1.0*isizee)
    vwgt(:) = 1
    vsize(:) = 1

    call metis_partmeshdual(imaxee, imaxnn, allnodesptr, allnodes, vwgt, vsize, ncommon1, isizee, tpwgts, options, objval, xmpiee, xmpidumb)

    do i = 1, imaxe
      xmpie(i) = xmpiee(i)
    end do

    allocate(pweight(0:isize - 1))
    pweight(:) = 0
    do i = 1, imaxe
      pweight(xmpiee(i)) = pweight(xmpiee(i)) + 1
    end do

    average = pweight(0)
    do i = 1, isize - 1
      average = pweight(i) + average
    end do
    average = average/isize
    maxi = maxval(pweight)

    open(63, file='history.txt', form='formatted', action='write', position='append')
    write(63, *) 'load imbalance of partioner', maxi/average
    close(63)

    deallocate(testar)
    deallocate(elements)
    deallocate(allnodes)
    deallocate(allnodesptr)
    deallocate(vwgt)
    deallocate(vsize)
    deallocate(tpwgts)
    deallocate(pweight)
    deallocate(xmpidumb)
    deallocate(xmpiee)
  end subroutine

  subroutine partitioner6(n, imaxe, imaxn, xmpie, ieshape)
    !> @brief
    !> this subroutine partitions the mesh using parmetis (preferred option large meshes can be partitioned)
    use iso_c_binding
    implicit none
    external parmetis_v3_partmeshkway
    integer, allocatable, dimension(:), intent(in) :: ieshape
    real :: average, average2, tsize
    integer :: maxi, maxi2
    integer :: i, j, elementid, node1, node2, node3, node4, node5, node6, node7, node8, counternodes, counternodes2, k
    integer :: vweight
    integer, allocatable, dimension(:, :, :) :: compweight, commweight
    integer, intent(in) :: n, imaxe, imaxn
    integer, allocatable, dimension(:), intent(inout) :: xmpie
    integer, allocatable, dimension(:) :: testar, pweight, pweight2
    integer(c_int) :: imaxee, imaxnn, ncommon1, isizee, objval, numflag, ncon, ncommonnodes, edgecut, wgtflag
    integer(c_int), allocatable, dimension(:) :: allnodesptr, allnodes, parts_el
    integer(c_int), allocatable, dimension(:) :: xmpiee, xmpidumb, vwgt, vsize
    integer(c_int), allocatable :: elmdist(:), eptr(:), eind(:), options(:)
    real(c_float), allocatable, dimension(:) :: tpwgts
    real(c_float), allocatable, dimension(:) :: ubvec

    wgtflag = 0
    imaxee = imaxe
    imaxnn = imaxn
    isizee = isize
    tsize = isizee
    allocate(xmpiee(1:isize))
    do i = 1, isize - 1
      xmpiee(i) = (imaxee/isize)
    end do
    xmpiee(isize) = imaxe - (xmpiee(isize - 1)*(isize - 1))

    allocate(elmdist(isize + 1))
    elmdist(1) = 1
    do i = 2, isize + 1
      elmdist(i) = elmdist(i - 1) + xmpiee(i - 1)
    end do

    if (dimensiona .eq. 3) then
      allnodesgloball = 8*xmpiee(n + 1)
    else
      allnodesgloball = 4*xmpiee(n + 1)
    end if

    allocate(allnodes(1:allnodesgloball))
    allocate(eptr(1:xmpiee(n + 1) + 1))
    counternodes = 0
    counternodes2 = 1
    eptr(counternodes2) = 1

    if (binio .eq. 0) then
      if (dimensiona .eq. 3) then
        open(1112, file='grid.cel', form='formatted', status='old', action='read')
        do i = 1, imaxe
          if ((i .ge. elmdist(n + 1)) .and. (i .lt. elmdist(n + 2))) then
            counternodes2 = counternodes2 + 1
            read(1112, *) elementid, node1, node2, node3, node4, node5, node6, node7, node8
            if ((node3 .eq. node4) .and. (node5 .ne. node6) .and. (node7 .eq. node8)) then ! prism
              counternodes = counternodes + 1
              allnodes(counternodes) = node1
              counternodes = counternodes + 1
              allnodes(counternodes) = node2
              counternodes = counternodes + 1
              allnodes(counternodes) = node3
              counternodes = counternodes + 1
              allnodes(counternodes) = node5
              counternodes = counternodes + 1
              allnodes(counternodes) = node6
              counternodes = counternodes + 1
              allnodes(counternodes) = node7
            else if ((node5 .ne. node6) .and. (node6 .ne. node7) .and. (node7 .ne. node8) .and. (node3 .ne. node4)) then ! hexa
              counternodes = counternodes + 1
              allnodes(counternodes) = node1
              counternodes = counternodes + 1
              allnodes(counternodes) = node2
              counternodes = counternodes + 1
              allnodes(counternodes) = node3
              counternodes = counternodes + 1
              allnodes(counternodes) = node4
              counternodes = counternodes + 1
              allnodes(counternodes) = node5
              counternodes = counternodes + 1
              allnodes(counternodes) = node6
              counternodes = counternodes + 1
              allnodes(counternodes) = node7
              counternodes = counternodes + 1
              allnodes(counternodes) = node8
            else if ((node3 .eq. node4) .and. (node5 .eq. node6) .and. (node6 .eq. node7) .and. (node7 .eq. node8)) then ! tetra
              counternodes = counternodes + 1
              allnodes(counternodes) = node1
              counternodes = counternodes + 1
              allnodes(counternodes) = node2
              counternodes = counternodes + 1
              allnodes(counternodes) = node3
              counternodes = counternodes + 1
              allnodes(counternodes) = node5
            else
              counternodes = counternodes + 1
              allnodes(counternodes) = node1
              counternodes = counternodes + 1
              allnodes(counternodes) = node2
              counternodes = counternodes + 1
              allnodes(counternodes) = node3
              counternodes = counternodes + 1
              allnodes(counternodes) = node4
              counternodes = counternodes + 1
              allnodes(counternodes) = node5
            end if
            eptr(counternodes2) = counternodes + 1
          else
            read(1112, *)
          end if
        end do
        close(1112)
      else
        open(1112, file='grid.cel', form='formatted', status='old', action='read')
        do i = 1, imaxe
          if ((i .ge. elmdist(n + 1)) .and. (i .lt. elmdist(n + 2))) then
            counternodes2 = counternodes2 + 1
            read(1112, *) elementid, node1, node2, node3, node4
            if ((node3 .eq. node4)) then ! triangle
              counternodes = counternodes + 1
              allnodes(counternodes) = node1
              counternodes = counternodes + 1
              allnodes(counternodes) = node2
              counternodes = counternodes + 1
              allnodes(counternodes) = node3
            else
              ! qudrilateral
              counternodes = counternodes + 1
              allnodes(counternodes) = node1
              counternodes = counternodes + 1
              allnodes(counternodes) = node2
              counternodes = counternodes + 1
              allnodes(counternodes) = node3
              counternodes = counternodes + 1
              allnodes(counternodes) = node4
            end if
            eptr(counternodes2) = counternodes + 1
          else
            read(1112, *)
          end if
        end do
        close(1112)
      end if
    else
      if (dimensiona .eq. 3) then
        open(1112, file='grid.cel', form='unformatted', status='old', action='read')
        do i = 1, imaxe
          if ((i .ge. elmdist(n + 1)) .and. (i .lt. elmdist(n + 2))) then
            counternodes2 = counternodes2 + 1
            read(1112) elementid, node1, node2, node3, node4, node5, node6, node7, node8
            if ((node3 .eq. node4) .and. (node5 .ne. node6) .and. (node7 .eq. node8)) then ! prism
              counternodes = counternodes + 1
              allnodes(counternodes) = node1
              counternodes = counternodes + 1
              allnodes(counternodes) = node2
              counternodes = counternodes + 1
              allnodes(counternodes) = node3
              counternodes = counternodes + 1
              allnodes(counternodes) = node5
              counternodes = counternodes + 1
              allnodes(counternodes) = node6
              counternodes = counternodes + 1
              allnodes(counternodes) = node7
            else if ((node5 .ne. node6) .and. (node6 .ne. node7) .and. (node7 .ne. node8) .and. (node3 .ne. node4)) then ! hexa
              counternodes = counternodes + 1
              allnodes(counternodes) = node1
              counternodes = counternodes + 1
              allnodes(counternodes) = node2
              counternodes = counternodes + 1
              allnodes(counternodes) = node3
              counternodes = counternodes + 1
              allnodes(counternodes) = node4
              counternodes = counternodes + 1
              allnodes(counternodes) = node5
              counternodes = counternodes + 1
              allnodes(counternodes) = node6
              counternodes = counternodes + 1
              allnodes(counternodes) = node7
              counternodes = counternodes + 1
              allnodes(counternodes) = node8
            else if ((node3 .eq. node4) .and. (node5 .eq. node6) .and. (node6 .eq. node7) .and. (node7 .eq. node8)) then ! tetra
              counternodes = counternodes + 1
              allnodes(counternodes) = node1
              counternodes = counternodes + 1
              allnodes(counternodes) = node2
              counternodes = counternodes + 1
              allnodes(counternodes) = node3
              counternodes = counternodes + 1
              allnodes(counternodes) = node5
            else
              counternodes = counternodes + 1
              allnodes(counternodes) = node1
              counternodes = counternodes + 1
              allnodes(counternodes) = node2
              counternodes = counternodes + 1
              allnodes(counternodes) = node3
              counternodes = counternodes + 1
              allnodes(counternodes) = node4
              counternodes = counternodes + 1
              allnodes(counternodes) = node5
            end if
            eptr(counternodes2) = counternodes + 1
          else
            read(1112)
          end if
        end do
        close(1112)
      else
        open(1112, file='grid.cel', form='unformatted', status='old', action='read')
        do i = 1, imaxe
          if ((i .ge. elmdist(n + 1)) .and. (i .lt. elmdist(n + 2))) then
            counternodes2 = counternodes2 + 1
            read(1112) elementid, node1, node2, node3, node4
            if ((node3 .eq. node4)) then ! triangle
              counternodes = counternodes + 1
              allnodes(counternodes) = node1
              counternodes = counternodes + 1
              allnodes(counternodes) = node2
              counternodes = counternodes + 1
              allnodes(counternodes) = node3
            else
              ! qudrilateral
              counternodes = counternodes + 1
              allnodes(counternodes) = node1
              counternodes = counternodes + 1
              allnodes(counternodes) = node2
              counternodes = counternodes + 1
              allnodes(counternodes) = node3
              counternodes = counternodes + 1
              allnodes(counternodes) = node4
            end if
            eptr(counternodes2) = counternodes + 1
          else
            read(1112)
          end if
        end do
        close(1112)
      end if
    end if

    allocate(eind(1:counternodes))
    eind(1:counternodes) = allnodes(1:counternodes)
    deallocate(allnodes)
    allocate(vwgt(xmpiee(n + 1)))
    vwgt = 0
    wgtflag = 0
    numflag = 1
    ncon = 1
    if (dimensiona .eq. 3) then
      ncommonnodes = 3
    else
      ncommonnodes = 2
    end if

    allocate(tpwgts(ncon*isize))
    tpwgts = 1./float(isize)
    allocate(ubvec(ncon))
    ubvec(1) = 1.050000
    allocate(options(3))
    options = [0, 0, 0]
    allocate(parts_el(xmpiee(n + 1)))
    parts_el = 0

    if (n .eq. 0) then
      write (*, *) '------------------------------------------------------------------------'
      write (*, *) '                       parmetis initiated                               '
      write (*, *) '------------------------------------------------------------------------'
    end if

    call mpi_barrier(mpi_comm_world, ierror)
    call parmetis_v3_partmeshkway(elmdist, eptr, eind, vwgt, wgtflag, numflag, ncon, ncommonnodes, isize, tpwgts, ubvec, options, edgecut, parts_el, mpi_comm_world)

    if (n .eq. 0) then
      write (*, *) '------------------------------------------------------------------------'
      write (*, *) '                       parmetis operation completed                     '
      write (*, *) '------------------------------------------------------------------------'
    end if

    deallocate(eptr)
    deallocate(eind)
    deallocate(vwgt)
    deallocate(tpwgts)
    deallocate(ubvec)
    deallocate(options)

    call mpi_gather(parts_el, xmpiee(1), mpi_integer, xmpie, xmpiee(1), mpi_integer, isize - 1, mpi_comm_world, ierror)
    if (n .eq. isize - 1) then
      xmpie(elmdist(n + 1):elmdist(n + 2) - 1) = parts_el(1:xmpiee(n + 1))
    end if

    call mpi_barrier(mpi_comm_world, ierror)
    call mpi_bcast(xmpie, imaxe, mpi_integer, isize - 1, mpi_comm_world, ierror)

    deallocate(parts_el)
    deallocate(xmpiee)
    deallocate(elmdist)
    xmpie = xmpie - 1
  end subroutine partitioner6
end module partition