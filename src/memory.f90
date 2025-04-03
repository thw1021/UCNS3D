module memory
  use mpi_info
  use declaration
  implicit none
contains
  subroutine allocate2
    ! @brief This subroutine allocates memory
    implicit none
    allocate(list(2000), ineb(6), iperb(6), nodelist(8))
  end subroutine allocate2

  subroutine allocate5
    ! @brief This subroutine allocates memory for the stencils
    implicit none
    allocate(ilocalallelg(n:n, xmpielrank(n), 1, iselemt(n)))
    ilocalallelg(:, :, :, :) = 0
    allocate(ilocalallelgper(n:n, xmpielrank(n), 1, iselemt(n)))
    ilocalallelgper(:, :, :, :) = 0
  end subroutine allocate5

  subroutine allocate3
    ! @brief This subroutine deallocates memory
    implicit none
    deallocate(list, ineb, iperb)
  end subroutine allocate3

  subroutine globaldea2(xmpil, xmpie)
    !> @brief This subroutine deallocates global lists
    implicit none
    integer, allocatable, dimension(:), intent(inout) :: xmpil, xmpie
    call mpi_barrier(mpi_comm_world, ierror)
    deallocate(xmpil)
    deallocate(xmpie)
    call mpi_barrier(mpi_comm_world, ierror)
  end subroutine globaldea2

  subroutine quad_alloc(numberofpoints, numberofpoints2)
    !> @brief This subroutine allocates memory for the quadrature points
    implicit none
    integer, intent(inout) :: numberofpoints, numberofpoints2
    integer :: i, kmaxe
    if (dimensiona .eq. 3) then
      numberofpoints = max(qp_hexa, qp_tetra, qp_pyra, qp_prism)
      if (dg .eq. 1) numberofpoints = max(qp_hexa, qp_tetra*6, qp_pyra, qp_prism)
      numberofpoints2 = max(qp_quad, qp_triangle, qp_triangle)
    else
      numberofpoints = max(qp_quad, qp_triangle)
      if (dg .eq. 1) numberofpoints = max(qp_quad, qp_triangle*2)
      numberofpoints2 = qp_line
    end if
    kmaxe = xmpielrank(n)
    do i = 1, kmaxe
      select case (ielem(n, i)%ishape)
      case (1) !hexa
        ielem(n, i)%itotalpoints = qp_tetra*6
      case (2) !tetra
        ielem(n, i)%itotalpoints = qp_tetra
      case (3) !pyramid
        ielem(n, i)%itotalpoints = qp_tetra*2
      case (4) !prism
        ielem(n, i)%itotalpoints = qp_tetra*3
      case (5) !quadrilateral
        ielem(n, i)%itotalpoints = qp_triangle*2
      case (6) !triangle
        ielem(n, i)%itotalpoints = qp_triangle
      end select
    end do
  end subroutine quad_alloc

  subroutine sumflux_allocation(n)
    !> @brief This subroutine allocates memory for the fluxes
    implicit none
    integer, intent(inout) :: n
    integer :: i, kmaxe
    kmaxe = xmpielrank(n)
    allocate(rhs(kmaxe))
    if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
      allocate(rhst(kmaxe))
    end if
    do i = 1, kmaxe
      if (dg .eq. 1) then
        allocate(rhs(i)%valdg(num_dg_dofs, nof_variables))
      end if
      allocate(rhs(i)%val(nof_variables))
      if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
        allocate(rhst(i)%val(turbulenceequations + passivescalar))
      end if
    end do
  end subroutine sumflux_allocation

  subroutine impallocate(n)
    !> @brief This subroutine allocates memory for implicit time stepping
    implicit none
    integer, intent(in) :: n
    integer :: kmaxe, interf
    kmaxe = xmpielrank(n)
    if (dimensiona .eq. 3) then
      interf = nof_variables
    else
      interf = nof_variables
    end if
    if (rungekutta .eq. 12) then
      allocate(impdu(kmaxe, 1:nof_variables + turbulenceequations + passivescalar))
      impdu(:, :) = zero
    else
      if (relax .eq. 3) then
        allocate(impdiag_mf(kmaxe))
        allocate(impoff_mf(kmaxe, interf))
        allocate(impdu(kmaxe, 1:nof_variables + turbulenceequations + passivescalar))
        if ((itestcase .eq. 4) .and. ((turbulence .gt. 0) .or. (passivescalar .gt. 0))) then
          allocate(impdiagt(kmaxe, turbulenceequations + passivescalar))
          allocate(impofft(kmaxe, interf, turbulenceequations + passivescalar))
          allocate(sht(kmaxe, turbulenceequations + passivescalar))
        end if
        impdiag_mf = zero
        impoff_mf = zero
        impdu = zero
      else
        if (dimensiona .eq. 3) then
          if (lowmemory .eq. 0) then
            allocate(impdiag(kmaxe, 1:nof_variables, 1:nof_variables))
            allocate(impoff(kmaxe, 6, 1:nof_variables, 1:nof_variables))
            allocate(impdu(kmaxe, 1:nof_variables + turbulenceequations + passivescalar))
            if ((itestcase .eq. 4) .and. ((turbulence .gt. 0) .or. (passivescalar .gt. 0))) then
              allocate(impofft(kmaxe, 6, turbulenceequations + passivescalar))
              allocate(impdiagt(kmaxe, turbulenceequations + passivescalar))
              allocate(sht(kmaxe, turbulenceequations + passivescalar))
            end if
            impdiag(:, :, :) = zero
            impoff(:, :, :, :) = zero
            impdu(:, :) = zero
          else
            allocate(impdu(kmaxe, 1:nof_variables + turbulenceequations + passivescalar))
            impdu(:, :) = zero
            if ((itestcase .eq. 4) .and. ((turbulence .gt. 0) .or. (passivescalar .gt. 0))) then
              allocate(sht(kmaxe, turbulenceequations + passivescalar))
            end if
          end if
        else
          if (lowmemory .eq. 0) then
            allocate(impdiag(kmaxe, 1:nof_variables, 1:nof_variables))
            allocate(impoff(kmaxe, 4, 1:nof_variables, 1:nof_variables))
            allocate(impdu(kmaxe, 1:nof_variables + turbulenceequations + passivescalar))
            if ((itestcase .eq. 4) .and. ((turbulence .gt. 0) .or. (passivescalar .gt. 0))) then
              allocate(impofft(kmaxe, 4, turbulenceequations + passivescalar))
              allocate(impdiagt(kmaxe, turbulenceequations + passivescalar))
              allocate(sht(kmaxe, turbulenceequations + passivescalar))
            end if
            impdiag(:, :, :) = zero
            impoff(:, :, :, :) = zero
            impdu(:, :) = zero
          else
            allocate(impdu(kmaxe, 1:nof_variables + turbulenceequations + passivescalar))
            if ((itestcase .eq. 4) .and. ((turbulence .gt. 0) .or. (passivescalar .gt. 0))) then
              allocate(sht(kmaxe, turbulenceequations + passivescalar))
            end if
            impdu(:, :) = zero
          end if
        end if
      end if
    end if
  end subroutine impallocate

  subroutine timing(n, cpux1, cpux2, cpux3, cpux4, cpux5, cpux6, timex1, timex2, timex3, timex4, timex5, timex6)
    !> @brief This subroutine allocates memory for the timers
    implicit none
    real, allocatable, dimension(:), intent(inout) :: cpux1, cpux2, cpux3, cpux4, cpux5, cpux6, timex1, timex2, timex3, timex4, timex5, timex6
    integer, intent(in) :: n
    allocate(cpux1(1))
    allocate(cpux2(1))
    allocate(cpux3(1))
    allocate(cpux4(1))
    allocate(cpux5(1))
    allocate(cpux6(1))
    allocate(timex1(1))
    allocate(timex2(1))
    allocate(timex3(1))
    allocate(timex4(1))
    allocate(timex5(1))
    allocate(timex6(1))
    cpux1(1) = 0.0; cpux2(1) = 0.0; cpux3(1) = 0.0; cpux4(1) = 0.0; cpux5(1) = 0.0; cpux6(1) = 0.0
    timex1(1) = 0.0; timex2(1) = 0.0; timex3(1) = 0.0; timex4(1) = 0.0; timex5(1) = 0.0; timex6(1) = 0.0
  end subroutine timing

  subroutine shallocation(ieshape, imaxe)
    !> @brief This subroutine allocates memory for the shapes
    implicit none
    integer, allocatable, dimension(:), intent(inout) :: ieshape
    integer, intent(inout) :: imaxe
    allocate(ieshape(imaxe))
    allocate(nodes_offset(imaxe), nodes_offset2(imaxe))
    ieshape = 0
    nodes_offset = 0
    nodes_offset2 = 0
  end subroutine shallocation

  subroutine shde_allocation(ieshape, imaxe)
    !> @brief This subroutine deallocates memory for the shapes
    implicit none
    integer, allocatable, dimension(:), intent(inout) :: ieshape
    integer, intent(inout) :: imaxe
    deallocate(ieshape, nodes_offset, nodes_offset2)
  end subroutine shde_allocation

  subroutine elallocation(n, xmpie, xmpielrank, ielem, imaxe, ieshape, itestcase, imaxb, ibound, xmin, xmax, ymin, ymax, zmin, zmax)
    !> @brief This subroutine allocates memory for the elements
    implicit none
    type(element_number), allocatable, dimension(:, :), intent(inout) :: ielem
    integer, allocatable, dimension(:), intent(in) :: ieshape
    integer, intent(in) :: n
    integer, intent(in) :: imaxe, itestcase, imaxb
    integer, allocatable, dimension(:), intent(in) :: xmpie
    integer, allocatable, dimension(:), intent(in) :: xmpielrank
    type(bound_number), allocatable, dimension(:, :) :: ibound
    integer :: i, j, k, lm, iex, kmaxe, kk
    real, allocatable, dimension(:), intent(inout) :: xmin, xmax, ymin, ymax, zmin, zmax
    allocate(xmin(n:n))
    allocate(ymin(n:n))
    allocate(zmin(n:n))
    allocate(xmax(n:n))
    allocate(ymax(n:n))
    allocate(zmax(n:n))
    kk = 0; i = 0; j = 0; lm = 0
    kmaxe = xmpielrank(n)
    allocate(ielem(n:n, xmpielrank(n)))
    if (tecplot .eq. 5) then
      allocate(nodes_offset_local(1:kmaxe)); nodes_offset_local = 0
      allocate(nodes_offset_local2(1:kmaxe)); nodes_offset_local2 = 0
    end if
    if (itestcase .lt. 3) then
      iex = 1
    end if
    if (itestcase .ge. 3) then
      if (dimensiona .eq. 3) then
        iex = 5
      else
        iex = 4
      end if
    end if
  end subroutine elallocation

  subroutine nodeallocation(n, inode, imaxn, xmpin, xmpinrank, inoden)
    !> @brief This subroutine allocates memory for the nodes
    implicit none
    integer, intent(in) :: n
    integer, allocatable, dimension(:), intent(in) :: xmpin
    integer, allocatable, dimension(:), intent(in) :: xmpinrank
    integer :: i, j, k, lm, iex, kmaxn, kk
    type(node_number), allocatable, dimension(:, :), intent(inout) :: inode
    type(node_ne), allocatable, dimension(:), intent(inout) :: inoden
    integer, intent(in) :: imaxn
  end subroutine nodeallocation

  subroutine nodedeallocation(n, inode, imaxn, xmpin, xmpinrank, inoden)
    !> @brief This subroutine deallocates memory from the nodes
    implicit none
    integer, intent(in) :: n
    integer, allocatable, dimension(:), intent(in) :: xmpin
    integer, allocatable, dimension(:), intent(in) :: xmpinrank
    integer :: i, j, k, lm, iex, kmaxn, kk
    type(node_number), allocatable, dimension(:, :), intent(inout) :: inode
    type(node_ne), allocatable, dimension(:), intent(inout) :: inoden
    integer, intent(in) :: imaxn
    deallocate(inode)
  end subroutine nodedeallocation

  subroutine xmpiallocate(xmpie, xmpil, xmpin, xmpielrank, xmpinrank, imaxe, imaxn, nproc)
    !> @brief This subroutine allocates memory for the global lists
    implicit none
    integer, allocatable, dimension(:), intent(inout) :: xmpie, xmpil
    integer, allocatable, dimension(:), intent(inout) :: xmpin
    integer, allocatable, dimension(:), intent(inout) :: xmpielrank
    integer, allocatable, dimension(:), intent(inout) :: xmpinrank
    integer, intent(inout) :: nproc, imaxe, imaxn
    allocate(xmpie(imaxe))
    allocate(xmpil(imaxe))
    allocate(xmpielrank(n:n))
    xmpie = 0
    xmpil = 0
    xmpielrank = 0
  end subroutine xmpiallocate

  subroutine deallocatempi1(n)
    !> @brief This subroutine deallocates memory for boundary exchange
    implicit none
    integer, intent(in) :: n
    deallocate(iexchanges1, iexchanger1)
  end subroutine deallocatempi1

  subroutine deallocatempi2(n)
    !> @brief This subroutine deallocates memory for the stencil exchange
    implicit none
    integer, intent(in) :: n
    deallocate(irecexr1, irecexs1)
  end subroutine deallocatempi2

  subroutine local_reconallocation3(n, ilocal_recon3)
    !> @brief This subroutine allocates memory for reconstruction prestoring
    implicit none
    type(local_recon3), allocatable, dimension(:), intent(inout) :: ilocal_recon3
    integer, intent(in) :: n
    integer :: i, j, k, m, ikg, itrue, kmaxe, idum
    integer :: inum_points, itarget
    integer :: imax, inum, ideg, imax2, inum2, ideg2, m2
    real :: perc, perde, perd, per1, per2, per3, per4, per5, per0, pef0, pef1, pef2, pef3, pef4, pef5, perv, per01, pef01
    kmaxe = xmpielrank(n)
    call mpi_barrier(mpi_comm_world, ierror)
    allocate(ilocal_recon3(kmaxe))
    if (dg .eq. 1) then
      allocate(ilocal_recon6(kmaxe))
      do i = 1, kmaxe        !for all elements
        allocate(ilocal_recon6(i)%dg2fv(1:idegfree, 1:nof_variables))
        ilocal_recon6(i)%dg2fv = 0.0d0
      end do
    end if
    perde = 0.0d0
    if (fastest .ne. 1) then
      do i = 1, kmaxe        !for all elements
        select case (ielem(n, i)%ishape)
        case (1, 2, 3, 4)
          imax = ielem(n, i)%inumneighbours - 1
          inum = ielem(n, i)%inumneighbours
          ideg = ielem(n, i)%idegfree
          m = ielem(n, i)%admis
          if (ees .eq. 5) then
            imax2 = numneighbours2 - 1
            inum2 = numneighbours2
            ideg2 = idegfree2
            m2 = ielem(n, i)%admis
          end if
          if (fastest .ne. 1) then
            allocate(ilocal_recon3(i)%invccjac(3, 3)); ilocal_recon3(i)%invccjac(:, :) = 0.0d0
            idum = 0
            if (ielem(n, i)%interior .eq. 1) then
              do j = 1, ielem(n, i)%ifca
                if (ielem(n, i)%ibounds(j) .gt. 0) then
                  if (ibound(n, ielem(n, i)%ibounds(j))%icode .eq. 4) then
                    idum = 1
                  end if
                end if
              end do
            end if
            if (idum .eq. 1) then
              allocate(ilocal_recon3(i)%volume(1, inum)); ilocal_recon3(i)%volume(:, :) = 0.0d0
            else
              allocate(ilocal_recon3(i)%volume(1, inum)); ilocal_recon3(i)%volume(:, :) = 0.0d0
            end if
            allocate(ilocal_recon3(i)%vext_ref(3)); ilocal_recon3(i)%vext_ref = 0.0d0
          end if
          if (firstorder .ne. 1) then
            if (greengo .eq. 0) then
              idum = 0
              if (ielem(n, i)%interior .eq. 1) then
                do j = 1, ielem(n, i)%ifca
                  if (ielem(n, i)%ibounds(j) .gt. 0) then
                    if (ibound(n, ielem(n, i)%ibounds(j))%icode .eq. 4) then
                      idum = 1
                    end if
                  end if
                end do
              end if
              if (idum .eq. 1) then
                allocate(ilocal_recon3(i)%stencils(m, imax, ideg)); ilocal_recon3(i)%stencils(:, :, :) = 0.0d0
                allocate(ilocal_recon3(i)%weightl(m, imax)); ilocal_recon3(i)%weightl(:, :) = 0.0d0!weightl
                if (ees .eq. 5) then
                  allocate(ilocal_recon3(i)%stencilsc(m2, imax2, ideg2)); ilocal_recon3(i)%stencilsc(:, :, :) = 0.0d0
                end if
              end if
            end if
            allocate(ilocal_recon3(i)%invmat_stencilt(ideg, imax, m)); ilocal_recon3(i)%invmat_stencilt(:, :, :) = 0.0d0
            if (ees .eq. 5) then
              allocate(ilocal_recon3(i)%invmat_stenciltc(ideg2, imax2, m2)); ilocal_recon3(i)%invmat_stenciltc(:, :, :) = 0.0d0
            end if
          end if
          allocate(ilocal_recon3(i)%ihexg(m, inum))
          allocate(ilocal_recon3(i)%ihexl(m, inum))
          allocate(ilocal_recon3(i)%periodicflag(m, inum))
          if (ees .eq. 5) then
            allocate(ilocal_recon3(i)%ihexgc(m2, inum2))
            allocate(ilocal_recon3(i)%ihexlc(m2, inum2))
          end if
          itrue = 0
          do j = 1, typesten
            ikg = 0
            if ((ees .ne. 5) .or. (j .eq. 1)) then
              itarget = inum
            else
              itarget = inum2
            end if
            do k = 1, itarget
              if (ilocalstencil(n, i, j, k) .gt. 0) then
                ikg = ikg + 1
                if (xmpie(ilocalstencil(n, i, j, k)) .ne. n) then
                  itrue = 1
                end if
              end if
            end do
          end do
          if (itrue .eq. 0) then
            ilocal_recon3(i)%local = 1
          else
            ilocal_recon3(i)%local = 0
          end if
          if (ilocal_recon3(i)%local .eq. 0) then
            allocate(ilocal_recon3(i)%ihexb(m, inum))
            allocate(ilocal_recon3(i)%ihexn(m, inum))
            if (ees .eq. 5) then
              allocate(ilocal_recon3(i)%ihexbc(m2, inum2))
              allocate(ilocal_recon3(i)%ihexnc(m2, inum2))
            end if
          end if
          if (iweno .eq. 1) then
            allocate(ilocal_recon3(i)%indicator(ideg, ideg)); ilocal_recon3(i)%indicator(:, :) = 0.0d0
            if (ees .eq. 5) then
              allocate(ilocal_recon3(i)%indicatorc(ideg2, ideg2)); ilocal_recon3(i)%indicatorc(:, :) = 0.0d0
            end if
          end if
        case (5, 6)
          imax = ielem(n, i)%inumneighbours - 1
          inum = ielem(n, i)%inumneighbours
          ideg = ielem(n, i)%idegfree
          m = ielem(n, i)%admis
          if (ees .eq. 5) then
            imax2 = numneighbours2 - 1
            inum2 = numneighbours2
            ideg2 = idegfree2
            m2 = ielem(n, i)%admis
          end if
          if (fastest .ne. 1) then
            allocate(ilocal_recon3(i)%invccjac(2, 2)); ilocal_recon3(i)%invccjac(:, :) = 0.0d0
            allocate(ilocal_recon3(i)%vext_ref(2)); ilocal_recon3(i)%vext_ref = 0.0d0
          end if
          idum = 0
          if (ielem(n, i)%interior .eq. 1) then
            do j = 1, ielem(n, i)%ifca
              if (ielem(n, i)%ibounds(j) .gt. 0) then
                if (ibound(n, ielem(n, i)%ibounds(j))%icode .eq. 4) then
                  idum = 1
                end if
              end if
            end do
          end if
          if (idum .eq. 1) then
            allocate(ilocal_recon3(i)%volume(1, inum)); ilocal_recon3(i)%volume(:, :) = 0.0d0
          else
            allocate(ilocal_recon3(i)%volume(1, inum)); ilocal_recon3(i)%volume(:, :) = 0.0d0
          end if
          if (firstorder .ne. 1) then
            if (greengo .eq. 0) then
              idum = 0
              if (ielem(n, i)%interior .eq. 1) then
                do j = 1, ielem(n, i)%ifca
                  if (ielem(n, i)%ibounds(j) .gt. 0) then
                    if (ibound(n, ielem(n, i)%ibounds(j))%icode .eq. 4) then
                      idum = 1
                    end if
                  end if
                end do
              end if
              if (idum .eq. 1) then
                allocate(ilocal_recon3(i)%stencils(m, imax, ideg)); ilocal_recon3(i)%stencils(:, :, :) = 0.0d0
                allocate(ilocal_recon3(i)%weightl(m, imax)); ilocal_recon3(i)%weightl(:, :) = 0.0d0
                if (ees .eq. 5) then
                  allocate(ilocal_recon3(i)%stencilsc(m2, imax2, ideg2)); ilocal_recon3(i)%stencilsc(:, :, :) = 0.0d0
                end if
              end if
            end if
            allocate(ilocal_recon3(i)%invmat_stencilt(ideg, imax, ielem(n, i)%admis)); ilocal_recon3(i)%invmat_stencilt(:, :, :) = 0.0d0
            if (ees .eq. 5) then
              allocate(ilocal_recon3(i)%invmat_stenciltc(ideg2, imax2, m2)); ilocal_recon3(i)%invmat_stenciltc(:, :, :) = 0.0d0
            end if
          end if
          allocate(ilocal_recon3(i)%ihexg(m, inum))
          allocate(ilocal_recon3(i)%ihexl(m, inum))
          if (initcond .eq. 0) then
            allocate(ilocal_recon3(i)%cond(7))
          end if
          if (ees .eq. 5) then
            allocate(ilocal_recon3(i)%ihexgc(m2, inum2))
            allocate(ilocal_recon3(i)%ihexlc(m2, inum2))
          end if
          itrue = 0
          do j = 1, typesten
            ikg = 0
            if ((ees .ne. 5) .or. (j .eq. 1)) then
              itarget = inum
            else
              itarget = inum2
            end if
            do k = 1, itarget
              if (ilocalstencil(n, i, j, k) .gt. 0) then
                ikg = ikg + 1
                if (xmpie(ilocalstencil(n, i, j, k)) .ne. n) then
                  itrue = 1
                end if
              end if
            end do
          end do
          if (itrue .eq. 0) then
            ilocal_recon3(i)%local = 1
          else
            ilocal_recon3(i)%local = 0
          end if
          if (ilocal_recon3(i)%local .eq. 0) then
            allocate(ilocal_recon3(i)%ihexb(m, inum))
            allocate(ilocal_recon3(i)%ihexn(m, inum))
            if (ees .eq. 5) then
              allocate(ilocal_recon3(i)%ihexbc(m2, inum2))
              allocate(ilocal_recon3(i)%ihexnc(m2, inum2))
            end if
          end if
          if (iweno .eq. 1) then
            allocate(ilocal_recon3(i)%indicator(ideg, ideg)); ilocal_recon3(i)%indicator(:, :) = 0.0d0
            if (ees .eq. 5) then
              allocate(ilocal_recon3(i)%indicatorc(ideg2, ideg2)); ilocal_recon3(i)%indicatorc(:, :) = 0.0d0
            end if
          end if
        end select
      end do
    end if
  end subroutine local_reconallocation3

  subroutine dealcordinates1(n, iexcordr, iexcords)
    !> @brief This subroutine deallocates memory for exchange of info between processes
    implicit none
    type(exchange_cord), allocatable, dimension(:), intent(inout) :: iexcordr
    type(exchange_cord), allocatable, dimension(:), intent(inout) :: iexcords
    integer, intent(in) :: n
    deallocate(iexcordr)
  end subroutine dealcordinates1

  subroutine dealcordinates2
    !> @brief This subroutine deallocates memory for exchange of info between processes
    implicit none
    if (allocated(iexcords)) deallocate(iexcords)
  end subroutine dealcordinates2

  subroutine allocate_basis_function(n, integ_basis, xmpielrank, idegfree)
    !> @brief This subroutine allocates memory for basis function integrals
    implicit none
    type(integralbasis), allocatable, dimension(:), intent(inout) :: integ_basis
    integer, allocatable, dimension(:), intent(in) :: xmpielrank
    integer, intent(in) :: idegfree, n
    integer :: kmaxe, i
    kmaxe = xmpielrank(n)
    allocate(integ_basis(kmaxe))
    if (dg .eq. 1) then
      allocate(integ_basis_dg(kmaxe))
      do i = 1, kmaxe
        allocate(integ_basis_dg(i)%value(1:idegfree)); integ_basis_dg(i)%value(:) = zero
      end do
    end if
    do i = 1, kmaxe
      allocate(integ_basis(i)%value(1:idegfree)); integ_basis(i)%value(:) = zero
      if (ees .eq. 5) then
        allocate(integ_basis(i)%valuec(1:idegfree2)); integ_basis(i)%valuec(:) = zero
      end if
    end do
  end subroutine allocate_basis_function

  subroutine deallocate_basis_function(n, integ_basis)
    !> @brief This subroutine deallocates memory for basis function integrals
    implicit none
    type(integralbasis), allocatable, dimension(:, :), intent(inout) :: integ_basis
    integer, intent(in) :: n
    deallocate(integ_basis)
  end subroutine deallocate_basis_function

  subroutine localsdeallocation(n, xmpielrank, ilocalstencil, ilocalstencilper, typesten, numneighbours)
    !> @brief This subroutine allocates memory for stencils
    implicit none
    integer, allocatable, dimension(:, :, :, :), intent(inout) :: ilocalstencil, ilocalstencilper
    integer, intent(in) :: numneighbours
    integer, intent(in) :: n
    integer, allocatable, dimension(:), intent(in) :: xmpielrank
    integer, intent(in) :: typesten
    deallocate(ilocalstencil)
    deallocate(ilocalstencilper)
  end subroutine localsdeallocation

  subroutine u_c_allocation(n, xmpielrank, u_c, u_e, itestcase, u_ct)
    !> @brief This subroutine allocates memory for solution vector
    implicit none
    type(u_centre), allocatable, dimension(:), intent(inout) :: u_c, u_ct
    type(u_exact), allocatable, dimension(:), intent(inout) :: u_e
    integer, allocatable, dimension(:), intent(in) :: xmpielrank
    integer, intent(in) :: itestcase, n
    integer :: i, kmaxe, istage, totalpoints, totalpointsvol
    kmaxe = xmpielrank(n)
    allocate(u_c(kmaxe))
    if (filtering .eq. 1) then
      allocate(u_cw(kmaxe))
      allocate(u_cs(kmaxe))
    end if
    if (itestcase .le. 4) then
      allocate(u_e(kmaxe))
    end if
    if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
      allocate(u_ct(kmaxe))
    end if
    select case (rungekutta)
    case (1)
      istage = 1
    case (2)
      istage = 2
    case (3)
      if (averaging .eq. 1) then
        istage = 5
      else
        istage = 3
        if (mood .eq. 1) then
          istage = 4
        end if
      end if
    case (4)
      if (averaging .eq. 1) then
        istage = 7
      else
        istage = 6
      end if
    case (5)
      istage = 2
    case (10)
      istage = 1
    case (11)
      if (averaging .eq. 1) then
        istage = 5
      else
        istage = 3
      end if
    case (12)
      if (averaging .eq. 1) then
        istage = 5
      else
        istage = 3
      end if
    end select
    if (dg .eq. 1) then
      allocate(m_1(kmaxe))
    end if
    do i = 1, kmaxe
      allocate(u_c(i)%val(istage, nof_variables)); u_c(i)%val = zero
      if (filtering .eq. 1) then
        allocate(u_cw(i)%val(1, nof_variables)); u_cw(i)%val = zero
        allocate(u_cs(i)%val(1, nof_variables)); u_cs(i)%val = zero
      end if
      if (dg .eq. 1) then
        allocate(u_c(i)%valdg(istage, nof_variables, ielem(n, i)%idegfree + 1)); u_c(i)%valdg = zero
        if (filtering .eq. 1) then
          allocate(u_cw(i)%valdg(1, nof_variables, ielem(n, i)%idegfree + 1)); u_cw(i)%valdg = zero
          allocate(u_cs(i)%valdg(1, nof_variables, ielem(n, i)%idegfree + 1)); u_cs(i)%valdg = zero
        end if
        allocate(m_1(i)%val(1:idegfree + 1, 1:idegfree + 1)); m_1(i)%val = zero
        if (itestcase .eq. 4) allocate(u_c(i)%br2_aux_var(ielem(n, i)%idegfree + 1, nof_variables, dimensiona)) ! ns
      end if
      if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
        allocate(u_ct(i)%val(istage, turbulenceequations + passivescalar)); u_ct(i)%val = zero
      end if
      if (averaging .eq. 1) then
        allocate(u_c(i)%rms(7))
        u_c(i)%rms(:) = zero
      end if
      if (itestcase .le. 4) then
        allocate(u_e(i)%val(1, nof_variables)); u_e(i)%val = zero
      end if
    end do
    if (mood .eq. 1) then
      do i = 1, kmaxe
        ielem(n, i)%recalc = 0
      end do
    else
      do i = 1, kmaxe
        ielem(n, i)%recalc = 1
      end do
    end if
  end subroutine u_c_allocation

  subroutine local_reconallocation5(n)
    !> @brief This subroutine allocates memory for reconstruction (one per process since these are destroyed after each element)
    implicit none
    integer, intent(in) :: n
    integer :: kmaxe, i
    kmaxe = xmpielrank(n)
    allocate(ilocal_recon5(1:kmaxe))
    do i = 1, kmaxe
      if (dimensiona .eq. 3) then
        allocate(ilocal_recon5(i)%gradients(typesten, idegfree, nof_variables)) !1000
        if (ees .eq. 5) then
          allocate(ilocal_recon5(i)%gradientsc(typesten, idegfree2, nof_variables)) !1000
        end if
        if ((turbulenceequations .gt. 0) .or. (passivescalar .gt. 0)) then
          allocate(ilocal_recon5(i)%gradients2(typesten, idegfree, 0 + turbulenceequations + passivescalar))!20
          if (ees .eq. 5) then
            allocate(ilocal_recon5(i)%gradientsc2(typesten, idegfree2, 0 + turbulenceequations + passivescalar))!20
          end if
        end if
        if (itestcase .eq. 4) then
          if ((turbulenceequations .gt. 0) .or. (passivescalar .gt. 0)) then
            allocate(ilocal_recon5(i)%gradientsturb(1, idegfree, 0 + turbulenceequations + passivescalar)) !10
          end if
          allocate(ilocal_recon5(i)%gradientstemp(idegfree))!10
          allocate(ilocal_recon5(i)%velocitydof(3, idegfree))!30
        end if
      else
        allocate(ilocal_recon5(i)%gradients(typesten, idegfree, nof_variables)) !1000
        if (ees .eq. 5) then
          allocate(ilocal_recon5(i)%gradientsc(typesten, idegfree2, nof_variables)) !1000
        end if
        if ((turbulenceequations .gt. 0) .or. (passivescalar .gt. 0)) then
          allocate(ilocal_recon5(i)%gradients2(typesten, idegfree, 0 + turbulenceequations + passivescalar))!20
          if (ees .eq. 5) then
            allocate(ilocal_recon5(i)%gradientsc2(typesten, idegfree2, 0 + turbulenceequations + passivescalar))!20
          end if
        end if
        if (itestcase .eq. 4) then
          if ((turbulenceequations .gt. 0) .or. (passivescalar .gt. 0)) then
            allocate(ilocal_recon5(i)%gradientsturb(1, idegfree, 0 + turbulenceequations + passivescalar)) !10
          end if
          allocate(ilocal_recon5(i)%gradientstemp(idegfree))!10
          allocate(ilocal_recon5(i)%velocitydof(2, idegfree))!30
        end if
      end if
    end do
  end subroutine local_reconallocation5

  subroutine local_reconallocation4(n)
    !> @brief This subroutine allocates memory for reconstruction
    implicit none
    integer, intent(in) :: n
    integer :: k, i, j, l, m, it, kmaxe, idum, iccf, decomf, svg, points, ii
    integer :: q5, q4, q3, q2, q1, q0, q01, icnn, iconsidered
    real :: perc, perde, perdi, per1, per2, per3, per4, per5, per0, pef0, pef1, pef2, pef3, pef4, pef5, perv, per01, pef01
    kmaxe = xmpielrank(n)
    if (itestcase .ge. 3) then
      it = nof_variables
    end if
    if (itestcase .lt. 3) then
      it = 1
    end if
    call mpi_barrier(mpi_comm_world, ierror)
    do ii = 1, nof_interior        !for all the interior elements
      i = el_int(ii)
      iconsidered = i
      select case (ielem(n, i)%ishape)
      case (1, 3, 4)
        points = qp_quad_n
      case (2)
        points = qp_triangle_n
      end select
      if (itestcase .eq. 4) then
        if (fastest .ne. 1) then
          allocate(ilocal_recon3(i)%uleftv(dims, it - 1, ielem(n, i)%ifca, points))
        else
          allocate(ilocal_recon3(i)%uleftv(dims, it - 1, ielem(n, i)%ifca, 1))
        end if
        ilocal_recon3(i)%uleftv = zero
        if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
          svg = (turbulenceequations + passivescalar)
          if (fastest .ne. 1) then
            allocate(ilocal_recon3(i)%uleftturbv(dims, svg, ielem(n, i)%ifca, points))        ! the derivatives of the turbulence model
            allocate(ilocal_recon3(i)%uleftturb(turbulenceequations + passivescalar, ielem(n, i)%ifca, points))
          else
            allocate(ilocal_recon3(i)%uleftturbv(dims, svg, ielem(n, i)%ifca, 1))        ! the derivatives of the turbulence model
            allocate(ilocal_recon3(i)%uleftturb(turbulenceequations + passivescalar, ielem(n, i)%ifca, 1))
          end if
          ilocal_recon3(i)%uleftturb = zero; ilocal_recon3(i)%uleftturbv = zero
        end if
      end if
      allocate(ilocal_recon3(i)%grads(4 + turbulenceequations + passivescalar + (qsas_model*3), 3))
      if ((outsurf .eq. 1) .and. (averaging .eq. 1)) then
        allocate(ilocal_recon3(i)%gradsav(4, 3))
      end if
      if (fastest .ne. 1) then
        allocate(ilocal_recon3(i)%uleft(it, ielem(n, i)%ifca, points))
      else
        allocate(ilocal_recon3(i)%uleft(it, ielem(n, i)%ifca, 1))
      end if
      ilocal_recon3(i)%uleft = zero
    end do
    do ii = 1, nof_bounded
      i = el_bnd(ii)
      iconsidered = i
      select case (ielem(n, i)%ishape)
      case (1, 3, 4)
        points = qp_quad_n
      case (2)
        points = qp_triangle_n
      end select
      if (itestcase .eq. 4) then
        if (fastest .ne. 1) then
          allocate(ilocal_recon3(i)%uleftv(dims, it - 1, ielem(n, i)%ifca, points))
        else
          allocate(ilocal_recon3(i)%uleftv(dims, it - 1, ielem(n, i)%ifca, 1))
        end if
        ilocal_recon3(i)%uleftv = zero
        if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
          svg = (turbulenceequations + passivescalar)
          if (fastest .ne. 1) then
            allocate(ilocal_recon3(i)%uleftturbv(dims, svg, ielem(n, i)%ifca, points))        ! the derivatives of the turbulence model
            allocate(ilocal_recon3(i)%uleftturb(turbulenceequations + passivescalar, ielem(n, i)%ifca, points))
          else
            allocate(ilocal_recon3(i)%uleftturbv(dims, svg, ielem(n, i)%ifca, 1))        ! the derivatives of the turbulence model
            allocate(ilocal_recon3(i)%uleftturb(turbulenceequations + passivescalar, ielem(n, i)%ifca, 1))
          end if
          ilocal_recon3(i)%uleftturb = zero; ilocal_recon3(i)%uleftturbv = zero
        end if
      end if
      allocate(ilocal_recon3(i)%grads(4 + turbulenceequations + passivescalar + (qsas_model*3), 3))
      if ((outsurf .eq. 1) .and. (averaging .eq. 1)) then
        allocate(ilocal_recon3(i)%gradsav(4, 3))
      end if
      if (fastest .ne. 1) then
        allocate(ilocal_recon3(i)%uleft(it, ielem(n, i)%ifca, points))
      else
        allocate(ilocal_recon3(i)%uleft(it, ielem(n, i)%ifca, 1))
      end if
      ilocal_recon3(i)%uleft = zero
    end do
  end subroutine local_reconallocation4

  subroutine local_reconallocation42d(n)
    !> @brief This subroutine allocates memory for reconstruction in 2d
    implicit none
    integer, intent(in) :: n
    integer :: k, i, j, l, m, it, kmaxe, idum, iccf, decomf, svg, points
    integer :: q5, q4, q3, q2, q1, q0, q01, icnn
    real :: perc, perde, perdi, per1, per2, per3, per4, per5, per0, pef0, pef1, pef2, pef3, pef4, pef5, perv, per01, pef01
    kmaxe = xmpielrank(n)
    if (itestcase .ge. 3) then
      it = nof_variables
    end if
    if (itestcase .lt. 3) then
      it = 1
    end if
    call mpi_barrier(mpi_comm_world, ierror)
    do i = 1, kmaxe
      points = qp_line_n
      if (itestcase .eq. 4) then !linear step?
        if (fastest .ne. 1) then
          allocate(ilocal_recon3(i)%uleftv(dims, it - 1, ielem(n, i)%ifca, points))
        else
          allocate(ilocal_recon3(i)%uleftv(dims, it - 1, ielem(n, i)%ifca, 1))
        end if
        ilocal_recon3(i)%uleftv = zero
        if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
          svg = (turbulenceequations + passivescalar)
          if (fastest .ne. 1) then
            allocate(ilocal_recon3(i)%uleftturbv(dims, svg, ielem(n, i)%ifca, points))        ! the derivatives of the turbulence model
            allocate(ilocal_recon3(i)%uleftturb(turbulenceequations + passivescalar, ielem(n, i)%ifca, points))
          else
            allocate(ilocal_recon3(i)%uleftturbv(dims, svg, ielem(n, i)%ifca, 1))        ! the derivatives of the turbulence model
            allocate(ilocal_recon3(i)%uleftturb(turbulenceequations + passivescalar, ielem(n, i)%ifca, 1))
          end if
          ilocal_recon3(i)%uleftturb = zero; ilocal_recon3(i)%uleftturbv = zero
        end if
        allocate(ilocal_recon3(i)%grads(3 + turbulenceequations + passivescalar + (qsas_model*2), 2))
      end if
      if (fastest .ne. 1) then
        allocate(ilocal_recon3(i)%uleft(it, ielem(n, i)%ifca, points))
        if (code_profile .eq. 2) then
          allocate(ilocal_recon3(i)%uleftx(it, ielem(n, i)%ifca, points))
        end if
      else
        allocate(ilocal_recon3(i)%uleft(it, ielem(n, i)%ifca, 1))
      end if
      ilocal_recon3(i)%uleft = zero
    end do
  end subroutine local_reconallocation42d
end module memory