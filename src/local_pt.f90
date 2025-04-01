module local
  use library
  use declaration
  use transform
  implicit none
contains
  subroutine exch_cords(n)
!> @brief
!> this subroutine establishes and communicates the exchange of coordinates for the boundary extrapolated values
    implicit none
    integer, intent(in)::n
    integer::i,j,k,l,ineedt,tneedt,indl,tndl,icpuid,ixflag,itee,iteedum,iavc,iavt,i_cnt,i_cnt2,i_cnt3,i_cnt4,ixf4,kmaxe,ixfv,i_cnt5
    real, dimension(1)::dumts, rumts
    integer, dimension(4)::icfv1, icfv2
    real::rcfv1, rcfv2
!-------------------for debugging only -----------------------------------------!

!-------------------for debugging only -----------------------------------------!
    kmaxe = xmpielrank(n)
    indl = iexchanger(1)%tot
    tndl = iexchanges(1)%tot
    if (fastest .ne. 1) then
      ineedt = irecexr(1)%tot
      tneedt = irecexs(1)%tot
      allocate (iexcordr(ineedt))
      allocate (iexcords(tneedt))
      allocate (iexsolhir(ineedt))
      allocate (iexsolhis(tneedt))
      if (adda .eq. 1) then
        allocate (iexsolhird(ineedt))
        allocate (iexsolhisd(tneedt))
      end if
    end if
    allocate (iexboundhir(indl))
    allocate (iexboundhis(tndl))
    allocate (iexboundhirr(indl))
    allocate (iexboundhiss(tndl))
    call mpi_barrier(mpi_comm_world, ierror)
    if (dimensiona .eq. 3) then
      i_cnt2 = 4; i_cnt3 = 3; i_cnt4 = 8
    else
      i_cnt2 = 2; i_cnt3 = 2; i_cnt4 = 4
    end if
    i_cnt5 = i_cnt3*i_cnt4
    if (fastest .ne. 1) then
    do i = 1, ineedt
      iexsolhir(i)%procid = irecexr(i)%procid
      iexcordr(i)%procid = irecexr(i)%procid
      allocate (iexcordr(i)%nodecord(irecexr(i)%muchineed(1), i_cnt4, i_cnt3))
      iexcordr(i)%nodecord(1:irecexr(i)%muchineed(1), i_cnt4, i_cnt3) = -tolbig
      allocate (iexsolhir(i)%sol(irecexr(i)%muchineed(1), nof_variables + turbulenceequations + passivescalar))
      if (adda .eq. 1) then
        allocate (iexsolhird(i)%sol(irecexr(i)%muchineed(1), 1))
      end if
      iexsolhir(i)%sol(:, :) = 0.0d0
    end do
    end if
    do i = 1, indl
      iexboundhir(i)%procid = iexchanger(i)%procid
      iexboundhirr(i)%procid = iexchanger(i)%procid
      if (itestcase .le. 3) then
        allocate (iexboundhirr(i)%vertpp(iexchanger(i)%muchineed(1), i_cnt2))
      else
        if (dimensiona .eq. 3) then
          i_cnt = (nof_variables + turbulenceequations + passivescalar) + ((4 + turbulenceequations + passivescalar)*3)
        else
          i_cnt = (nof_variables + turbulenceequations + passivescalar) + ((3 + turbulenceequations + passivescalar)*2)
        end if
        allocate (iexboundhirr(i)%vertpp(iexchanger(i)%muchineed(1), i_cnt2))
      end if
      iexboundhirr(i)%vertpp(:, :) = 0
    end do
    if (fastest .ne. 1) then
    do i = 1, tneedt
      iexsolhis(i)%procid = irecexs(i)%procid
      iexcords(i)%procid = irecexs(i)%procid
      allocate (iexcords(i)%nodecord(irecexs(i)%muchtheyneed(1), i_cnt4, i_cnt3))
      iexcords(i)%nodecord(1:irecexs(i)%muchtheyneed(1), 1:i_cnt4, 1:i_cnt3) = -tolbig
      allocate (iexsolhis(i)%sol(irecexs(i)%muchtheyneed(1), nof_variables + turbulenceequations + passivescalar))
      if (adda .eq. 1) then
        allocate (iexsolhisd(i)%sol(irecexs(i)%muchtheyneed(1), 1))
      end if
      iexsolhis(i)%sol(:, :) = 0.0d0
    end do
    end if
    do i = 1, tndl
      iexboundhis(i)%procid = iexchanges(i)%procid
      iexboundhiss(i)%procid = iexchanges(i)%procid
      if (itestcase .le. 3) then
        allocate (iexboundhiss(i)%vertpp(iexchanges(i)%muchtheyneed(1), i_cnt2))
      else
        if (dimensiona .eq. 3) then
          i_cnt = (nof_variables + turbulenceequations + passivescalar) + ((4 + turbulenceequations + passivescalar)*3)
        else
          i_cnt = (nof_variables + turbulenceequations + passivescalar) + ((3 + turbulenceequations + passivescalar)*2)
        end if
        allocate (iexboundhiss(i)%vertpp(iexchanges(i)%muchtheyneed(1), i_cnt2))
      end if
      iexboundhiss(i)%vertpp(:, :) = 0
    end do
    call exchange_cordx(n, ineedt, tneedt, i_cnt4, i_cnt3)
  end subroutine exch_cords
  subroutine exchange_cordx(n, ineedt, tneedt, i_cnt4, i_cnt3)
    integer, intent(in)::n, ineedt, tneedt, i_cnt4, i_cnt3
    integer::i, j, k, l, icpuid, ixflag, itee, iteedum, iavc, iavt
    real, dimension(1)::dumts, rumts
    if (fastest .ne. 1) then
      do i = 1, tneedt
          end do
        end do
      end do
      call mpi_barrier(mpi_comm_world, ierror)
      icpuid = n
      dumts(1:1) = tolsmall
      do i = 0, isize - 1
        if (i .ne. n) then
          do j = 1, tneedt
            iavt = 10000
            if (irecexs(j)%procid .eq. i) then
              iavt = j
              go to 7001
            end if
          end do
7001      continue
          do k = 1, ineedt
            iavc = 10000
            if (irecexr(k)%procid .eq. i) then
              iavc = k
              go to 8001
            end if
          end do
8001      continue
          if ((iavt .eq. 10000) .and. (iavc .ne. 10000)) then
            call mpi_sendrecv(dumts(1:1), 1, mpi_double_precision, i, icpuid, &
                              iexcordr(iavc)%nodecord(1:irecexr(iavc)%muchineed(1), 1:i_cnt4, 1:i_cnt3), &
irecexr(iavc)%muchineed(1)*i_cnt4*i_cnt3,mpi_double_precision,iexcordr(iavc)%procid,iexcordr(iavc)%procid,mpi_comm_world,status,ierror)
          end if
          if ((iavt .ne. 10000) .and. (iavc .eq. 10000)) then
            !message 1
            call mpi_sendrecv(iexcords(iavt)%nodecord(1:irecexs(iavt)%muchtheyneed(1), 1:i_cnt4, 1:i_cnt3), &
                              irecexs(iavt)%muchtheyneed(1)*i_cnt4*i_cnt3, mpi_double_precision, iexcords(iavt)%procid, icpuid, &
                              dumts(1:1), 1, mpi_double_precision, i, i, mpi_comm_world, status, ierror)
          end if
          if ((iavt .ne. 10000) .and. (iavc .ne. 10000)) then
            !message 1
            call mpi_sendrecv(iexcords(iavt)%nodecord(1:irecexs(iavt)%muchtheyneed(1), 1:i_cnt4, 1:i_cnt3), &
                              irecexs(iavt)%muchtheyneed(1)*i_cnt4*i_cnt3, mpi_double_precision, iexcords(iavt)%procid, icpuid, &
              iexcordr(iavc)%nodecord(1:irecexr(iavc)%muchineed(1), 1:i_cnt4, 1:i_cnt3), irecexr(iavc)%muchineed(1)*i_cnt4*i_cnt3, &
                              mpi_double_precision, iexcordr(iavc)%procid, iexcordr(iavc)%procid, mpi_comm_world, status, ierror)
          end if
          !end do! j
        end if        ! i.ne.n
      end do
    end if
    call mpi_barrier(mpi_comm_world, ierror)
  end subroutine exchange_cordx
  subroutine exch_cords_opt(n)
!> @brief
!> this subroutine establishes and communicates the exchange of coordinates for the boundary extrapolated values
    implicit none
    integer, intent(in)::n
    integer::i,j,k,l,ineedt,tneedt,indl,tndl,icpuid,ixflag,itee,iteedum,iavc,iavt,i_cnt,i_cnt2,i_cnt3,i_cnt4,ixf4,kmaxe,ixfv,i_cnt5
    real, dimension(1)::dumts, rumts
    integer, dimension(4)::icfv1, icfv2
    real::rcfv1, rcfv2
!-------------------for debugging only -----------------------------------------!
!-------------------for debugging only -----------------------------------------!
    kmaxe = xmpielrank(n)
    indl = iexchanger(1)%tot
    tndl = iexchanges(1)%tot
    if (dimensiona .eq. 3) then
      i_cnt2 = 4; i_cnt3 = 3; i_cnt4 = 8
    else
      i_cnt2 = 2; i_cnt3 = 2; i_cnt4 = 4
    end if
    i_cnt5 = i_cnt3*i_cnt4
    do i = 1, indl
      if (itestcase .le. 3) then
        allocate (iexboundhir(i)%facesol(iexchanger(i)%muchineed(1), nof_variables))
        if (dg .eq. 1) then
          allocate (iexboundhir(i)%facesol_dg(iexchanger(i)%muchineed(1), nof_variables))
        end if
        if (mood .eq. 1) then
          allocate (iexboundhir(i)%facesol_m(iexchanger(i)%muchineed(1), 1))
        end if
      else
        if (dimensiona .eq. 3) then
          i_cnt = (nof_variables + turbulenceequations + passivescalar) + ((4 + turbulenceequations + passivescalar)*3)
        else
          i_cnt = (nof_variables + turbulenceequations + passivescalar) + ((3 + turbulenceequations + passivescalar)*2)
        end if
        allocate (iexboundhir(i)%facesol(iexchanger(i)%muchineed(1), i_cnt))
        if (dg .eq. 1) then
          allocate (iexboundhir(i)%facesol_dg(iexchanger(i)%muchineed(1), i_cnt))
        end if
      end if
      iexboundhir(i)%facesol(:, :) = 0.0d0
      if (dg .eq. 1) then
        iexboundhir(i)%facesol_dg(:, :) = 0.0d0
      end if
    end do
    do i = 1, tndl
      if (itestcase .le. 3) then
        allocate (iexboundhis(i)%facesol(iexchanges(i)%muchtheyneed(1), nof_variables))
        if (dg .eq. 1) then
          allocate (iexboundhis(i)%facesol_dg(iexchanges(i)%muchtheyneed(1), nof_variables))
        end if
        if (mood .eq. 1) then
          allocate (iexboundhis(i)%facesol_m(iexchanges(i)%muchtheyneed(1), 1))
        end if
      else
        if (dimensiona .eq. 3) then
          i_cnt = (nof_variables + turbulenceequations + passivescalar) + ((4 + turbulenceequations + passivescalar)*3)
        else
          i_cnt = (nof_variables + turbulenceequations + passivescalar) + ((3 + turbulenceequations + passivescalar)*2)
        end if
        allocate (iexboundhis(i)%facesol(iexchanges(i)%muchtheyneed(1), i_cnt))
        if (dg .eq. 1) then
          allocate (iexboundhis(i)%facesol_dg(iexchanges(i)%muchtheyneed(1), i_cnt))
        end if
      end if
      iexboundhis(i)%facesol(:, :) = 0.0d0
      if (dg .eq. 1) then
        iexboundhis(i)%facesol_dg(:, :) = 0.0d0
      end if
    end do
    call mpi_barrier(mpi_comm_world, ierror)
  end subroutine exch_cords_opt
  subroutine exch_cord3(n)
!> @brief
!> this subroutine establishes and communicates the exchange of coordinates for the boundary extrapolated values and mapping of gaussian quadrature points
    implicit none
    integer, intent(in)::n
    integer::i,j,k,l,ineedt,tneedt,indl,tndl,icpuid,ixflag,itee,iteedum,iavc,iavt,i_cnt,i_cnt2,i_cnt3,i_cnt4,ixf4,kmaxe,ixfv,i_cnt5
    real, dimension(1)::dumts, rumts
    integer, dimension(4)::icfv1, icfv2
    real::rcfv1, rcfv2
    icpuid = n
    kmaxe = xmpielrank(n)
    indl = iexchanger(1)%tot
    tndl = iexchanges(1)%tot
    if (dimensiona .eq. 3) then
      i_cnt2 = 4; i_cnt3 = 3; i_cnt4 = 8
    else
      i_cnt2 = 2; i_cnt3 = 2; i_cnt4 = 4
    end if
    i_cnt5 = i_cnt3*i_cnt4
    if (dimensiona .eq. 3) then
      do j = 1, indl
      do i = 1, tndl
      if (iexchanger(j)%procid .eq. iexchanges(i)%procid) then
      do k = 1, iexchanges(i)%muchtheyneed(1)
        if (ielem(n, (iexchanges(i)%localref(k)))%types_faces(iexchanges(i)%sidetheyneed(k)) .eq. 5) then
          ixf4 = 4
        else
          ixf4 = 3
        end if
        iexboundhiss(i)%vertpp(k, 1:ixf4) = ielem(n, (iexchanges(i)%localref(k)))%nodes_faces(iexchanges(i)%sidetheyneed(k), 1:ixf4)
      end do
      end if
      end do
      end do
      call mpi_barrier(mpi_comm_world, ierror)
      do k = 1, indl
      do j = 1, tndl
        if (iexchanger(k)%procid .eq. iexchanges(j)%procid) then
          call mpi_sendrecv(iexboundhiss(j)%vertpp(1:iexchanges(j)%muchtheyneed(1), 1:i_cnt2) &
                            , iexchanges(j)%muchtheyneed(1)*i_cnt2, mpi_integer, iexchanges(j)%procid, &
                            iexchanges(j)%procid, iexboundhirr(k)%vertpp(1:iexchanger(k)%muchineed(1), 1:i_cnt2), &
                            iexchanger(k)%muchineed(1)*i_cnt2, mpi_integer, iexchanges(j)%procid, &
                            icpuid, mpi_comm_world, status, ierror)
        end if
      end do
      end do
      do i = 1, kmaxe
        if (ielem(n, i)%interior .eq. 1) then
          do k = 1, ielem(n, i)%ifca
            if (ielem(n, i)%types_faces(k) .eq. 5) then
              ixf4 = 4
            else
              ixf4 = 3
            end if
            if (ielem(n, i)%ineighg(k) .gt. 0) then
              if (ielem(n, i)%ineighb(k) .ne. n) then
                if (ielem(n, i)%ineighg(k) .gt. ielem(n, i)%ihexgl) then
                  if (ielem(n, i)%ibounds(k) .gt. 0) then
                    if ((ibound(n, ielem(n, i)%ibounds(k))%icode .eq. 5) .or. (ibound(n, ielem(n, i)%ibounds(k))%icode .eq. 50)) then
                      do ixfv = 1, ixf4
                  inoder(iexboundhirr(ielem(n,i)%ineighn(k))%vertpp(ielem(n,i)%q_face(k)%q_mapl(1),ixfv))%itor=iexboundhirr(ielem(n,i)%ineighn(k))%vertpp(ielem(n,i)%q_face(k)%q_mapl(1),ixfv)
                      end do
                    end if
                  end if
                end if
              end if
            end if
          end do
        end if
      end do
    else
!now try to remap the local and global neighbours and mapping of the gaussian quadrature points
!first the global
      do j = 1, indl
      do i = 1, tndl
      if (iexchanger(j)%procid .eq. iexchanges(i)%procid) then
      do k = 1, iexchanges(i)%muchtheyneed(1)
        ixf4 = 2
        iexboundhiss(i)%vertpp(k, 1:ixf4) = ielem(n, (iexchanges(i)%localref(k)))%nodes_faces(iexchanges(i)%sidetheyneed(k), 1:ixf4)
      end do
      end if
      end do
      end do
      do k = 1, indl
      do j = 1, tndl
        if (iexboundhirr(k)%procid .eq. iexboundhiss(j)%procid) then
          call mpi_sendrecv(iexboundhiss(j)%vertpp(1:iexchanges(j)%muchtheyneed(1), 1:i_cnt2) &
                            , iexchanges(j)%muchtheyneed(1)*i_cnt2, mpi_integer, iexboundhiss(j)%procid, &
                            icpuid, iexboundhirr(k)%vertpp(1:iexchanger(k)%muchineed(1), 1:i_cnt2), &
                            iexchanger(k)%muchineed(1)*i_cnt2, mpi_integer, iexboundhirr(k)%procid, &
                            iexboundhirr(k)%procid, mpi_comm_world, status, ierror)
        end if
      end do
      end do
      do i = 1, kmaxe
        if (ielem(n, i)%interior .eq. 1) then
          do k = 1, ielem(n, i)%ifca
            ixf4 = 2
            if (ielem(n, i)%ineighg(k) .gt. 0) then
              if (ielem(n, i)%ineighb(k) .ne. n) then
                if (ielem(n, i)%ineighg(k) .gt. ielem(n, i)%ihexgl) then
                  if (ielem(n, i)%ibounds(k) .gt. 0) then
                    if ((ibound(n, ielem(n, i)%ibounds(k))%icode .eq. 5) .or. (ibound(n, ielem(n, i)%ibounds(k))%icode .eq. 50)) then
                      do ixfv = 1, ixf4
                  inoder(iexboundhirr(ielem(n,i)%ineighn(k))%vertpp(ielem(n,i)%q_face(k)%q_mapl(1),ixfv))%itor=iexboundhirr(ielem(n,i)%ineighn(k))%vertpp(ielem(n,i)%q_face(k)%q_mapl(1),ixfv)
                      end do
                    end if
                  end if
                end if
              else
              end if
            end if
          end do
        else
        end if
      end do
    end if
  end subroutine exch_cord3
  subroutine exch_cords2(n, isize, iexboundhiri, iexboundhisi, &
                         itestcase, numberofpoints2, iexchanger, iexchanges)
!> @brief
!> this subroutine establishes and communicates the exchange of coordinates for the boundary extrapolated values in 2d
    implicit none
    type(exchange_boundhi), allocatable, dimension(:), intent(inout)::iexboundhiri
    type(exchange_boundhi), allocatable, dimension(:), intent(inout)::iexboundhisi
    type(exchange), allocatable, dimension(:), intent(in)::iexchanger
    type(exchange), allocatable, dimension(:), intent(in)::iexchanges
    integer, intent(in)::isize, n, itestcase, numberofpoints2
    integer::i, j, k, l, ineedt, tneedt, indl, tndl, icpuid, ixflag, itee, iteedum, iavc, iavt, i_cnt
    real, dimension(1:1)::dumts, rumts
    indl = iexchanger(1)%tot
    tndl = iexchanges(1)%tot
    allocate (iexboundhiri(indl))
    allocate (iexboundhisi(tndl))
    i_cnt = (nof_variables + turbulenceequations + passivescalar)
    do i = 1, indl
      iexboundhiri(i)%procid = iexchanger(i)%procid
      allocate (iexboundhiri(i)%facesol(iexchanger(i)%muchineed(1), i_cnt))
      iexboundhiri(i)%facesol(:, :) = 0.0d0
    end do
    do i = 1, tndl
      iexboundhisi(i)%procid = iexchanges(i)%procid
      allocate (iexboundhisi(i)%facesol(iexchanges(i)%muchtheyneed(1), i_cnt))
      iexboundhisi(i)%facesol(:, :) = 0.0d0
    end do
    call mpi_barrier(mpi_comm_world, ierror)
  end subroutine exch_cords2
subroutine localise_stencil(n,iconsi,ilox_ihexg,ilox_ihexl,ilox_ihexb,ilox_ihexn,ilox_ishape,ilox_xxc,ilox_yyc,ilox_zzc,ilox_volume,ilox_periodicflag,ilon_nodcount,ilon_x,ilon_y,ilon_z)
!> @brief
!> this subroutine starts expressing all the stencil elements coordinates and volumes with respect to the considered cell
    implicit none
    integer, intent(in)::n, iconsi
    integer::i, j, l, ixff, ixsst, ikg, ismp, inv, ineedt, ikg2, in_sten, k, itarget
    integer, allocatable, dimension(:, :), intent(inout)::ilox_ihexg  !global index of cells
    integer, allocatable, dimension(:, :), intent(inout)::ilox_ihexl  !local index of cells
    integer, allocatable, dimension(:, :), intent(inout)::ilox_ihexb  !cpu that that each cell belongs to
    integer, allocatable, dimension(:, :), intent(inout)::ilox_ihexn  !internal index from where to take the values from communicated messages
    integer, allocatable, dimension(:, :), intent(inout)::ilox_ishape !shape of each element
    real, allocatable, dimension(:, :), intent(inout)::ilox_xxc       !cell centre coordinates in x
    real, allocatable, dimension(:, :), intent(inout)::ilox_yyc       !cell centre coordinates in y
    real, allocatable, dimension(:, :), intent(inout)::ilox_zzc      !cell centre coordinates in z
    real, allocatable, dimension(:, :), intent(inout)::ilox_volume    !cell volume
    integer, allocatable, dimension(:, :), intent(inout)::ilox_periodicflag
    integer, allocatable, dimension(:, :, :), intent(inout)::ilon_nodcount  !number of nodes
    real, allocatable, dimension(:, :, :), intent(inout)::ilon_x           !coordinates of each node in x
    real, allocatable, dimension(:, :, :), intent(inout)::ilon_y           !coordinates of each node in x
    real, allocatable, dimension(:, :, :), intent(inout)::ilon_z           !coordinates of each node in x
    real::rin_sten
    real, dimension(1:dimensiona)::cords
    ineedt = irecexr(1)%tot
    i = iconsi
    ikg2 = 0
    do ismp = 1, typesten
      ikg = 0
      if ((ees .ne. 5) .or. (ismp .eq. 1)) then
        itarget = ielem(n, i)%inumneighbours
      else
        itarget = numneighbours2
      end if
      do l = 1, itarget
        if (ilocalstencil(n, i, ismp, l) .gt. 0) then
          ikg = ikg + 1
        end if
      end do
      if (ikg .eq. itarget) then
        if (ilocal_recon3(i)%local .eq. 1) then
          ikg2 = ikg2 + 1
          do l = 1, itarget
            ilox_ihexg(ikg2, l) = ilocalstencil(n, i, ismp, l)
            ilox_periodicflag(ikg2, l) = ilocalstencilper(n, i, ismp, l)
            j = (xmpil(ilox_ihexg(ikg2, l)))
            ilox_ihexl(ikg2, l) = j
            ilox_ishape(ikg2, l) = ielem(n, j)%ishape
            call compute_centre3d(j, cords)
            ilox_xxc(ikg2, l) = cords(1)
            ilox_yyc(ikg2, l) = cords(2)
            ilox_zzc(ikg2, l) = cords(3)
            ilon_nodcount(ikg2, l, 1:ielem(n, j)%nonodes) = ielem(n, j)%nodes(1:ielem(n, j)%nonodes)
            do k = 1, ielem(n, j)%nonodes
              ilon_x(ikg2, l, k) = inoder(ielem(n, j)%nodes(k))%cord(1)
              ilon_y(ikg2, l, k) = inoder(ielem(n, j)%nodes(k))%cord(2)
              ilon_z(ikg2, l, k) = inoder(ielem(n, j)%nodes(k))%cord(3)
            end do
          end do
        else
          ikg2 = ikg2 + 1
          do l = 1, itarget
            ilox_ihexg(ikg2, l) = ilocalstencil(n, i, ismp, l)
            ilox_periodicflag(ikg2, l) = ilocalstencilper(n, i, ismp, l)
            ilox_ihexb(ikg2, l) = xmpie(ilocalstencil(n, i, ismp, l))
            if (ilox_ihexb(ikg2, l) .eq. n) then
              j = (xmpil(ilox_ihexg(ikg2, l)))
              if (ilox_ihexg(ikg2, l) .eq. ielem(n, j)%ihexgl) then
                ilox_ihexl(ikg2, l) = j
                ilox_ishape(ikg2, l) = ielem(n, j)%ishape
              end if
              call compute_centre3d(j, cords)
              ilox_xxc(ikg2, l) = cords(1)
              ilox_yyc(ikg2, l) = cords(2)
              ilox_zzc(ikg2, l) = cords(3)
              ilon_nodcount(ikg2, l, 1:ielem(n, j)%nonodes) = ielem(n, j)%nodes(1:ielem(n, j)%nonodes)
              do k = 1, ielem(n, j)%nonodes
                ilon_x(ikg2, l, k) = inoder(ielem(n, j)%nodes(k))%cord(1)
                ilon_y(ikg2, l, k) = inoder(ielem(n, j)%nodes(k))%cord(2)
                ilon_z(ikg2, l, k) = inoder(ielem(n, j)%nodes(k))%cord(3)
              end do
            else
              do ixff = 1, ineedt
                if (iexcordr(ixff)%procid .eq. ilox_ihexb(ikg2, l)) then
                  do ixsst = 1, irecexr(ixff)%muchineed(1)
                    if ((ilox_ihexg(ikg2, l)) .eq. irecexr1(ixff)%whatineed(ixsst)) then
                      ilox_ihexl(ikg2, l) = ixsst
                      ilox_ihexn(ikg2, l) = ixff
                      ilox_ishape(ikg2, l) = irecexr1(ixff)%ishape(ixsst)
                      exit
                    end if
                  end do
                end if
              end do
              select case (ilox_ishape(ikg2, l))
              case (1)
                in_sten = 8
                rin_sten = 8.0d0
              case (2)
                in_sten = 4
                rin_sten = 4.0d0

              case (3)
                in_sten = 5
                rin_sten = 5.0d0

              case (4)
                in_sten = 6
                rin_sten = 6.0d0
              end select
              ilon_x(ikg2, l, 1:in_sten) = iexcordr(ilox_ihexn(ikg2, l))%nodecord(ilox_ihexl(ikg2, l), 1:in_sten, 1)
              ilon_y(ikg2, l, 1:in_sten) = iexcordr(ilox_ihexn(ikg2, l))%nodecord(ilox_ihexl(ikg2, l), 1:in_sten, 2)
              ilon_z(ikg2, l, 1:in_sten) = iexcordr(ilox_ihexn(ikg2, l))%nodecord(ilox_ihexl(ikg2, l), 1:in_sten, 3)
              ilox_xxc(ikg2, l) = sum(ilon_x(ikg2, l, 1:in_sten))/rin_sten
              ilox_yyc(ikg2, l) = sum(ilon_y(ikg2, l, 1:in_sten))/rin_sten
              ilox_zzc(ikg2, l) = sum(ilon_z(ikg2, l, 1:in_sten))/rin_sten
            end if
          end do
        end if
      end if
    end do
  end subroutine localise_stencil
subroutine localise_stencil2d(n,iconsi,ilox_ihexg,ilox_ihexl,ilox_ihexb,ilox_ihexn,ilox_ishape,ilox_xxc,ilox_yyc,ilox_zzc,ilox_volume,ilox_periodicflag,ilon_nodcount,ilon_x,ilon_y,ilon_z)
!> @brief
!> this subroutine starts expressing all the stencil elements coordinates and volumes with respect to the considered cell in 2d
    implicit none
    integer, intent(in)::n, iconsi
    integer, allocatable, dimension(:, :), intent(inout)::ilox_ihexg  !global index of cells
    integer, allocatable, dimension(:, :), intent(inout)::ilox_ihexl  !local index of cells
    integer, allocatable, dimension(:, :), intent(inout)::ilox_ihexb  !cpu that that each cell belongs to
    integer, allocatable, dimension(:, :), intent(inout)::ilox_ihexn  !internal index from where to take the values from communicated messages
    integer, allocatable, dimension(:, :), intent(inout)::ilox_ishape !shape of each element
    real, allocatable, dimension(:, :), intent(inout)::ilox_xxc       !cell centre coordinates in x
    real, allocatable, dimension(:, :), intent(inout)::ilox_yyc       !cell centre coordinates in y
    real, allocatable, dimension(:, :), intent(inout)::ilox_zzc      !cell centre coordinates in z
    real, allocatable, dimension(:, :), intent(inout)::ilox_volume    !cell volume
    integer, allocatable, dimension(:, :), intent(inout)::ilox_periodicflag
    integer, allocatable, dimension(:, :, :), intent(inout)::ilon_nodcount  !number of nodes
    real, allocatable, dimension(:, :, :), intent(inout)::ilon_x           !coordinates of each node in x
    real, allocatable, dimension(:, :, :), intent(inout)::ilon_y           !coordinates of each node in x
    real, allocatable, dimension(:, :, :), intent(inout)::ilon_z           !coordinates of each node in x
    real, dimension(1:dimensiona)::cords
    integer::i, j, l, ixff, ixsst, ikg, ismp, inv, ineedt, ikg2, in_sten, k, itarget, nojcount
    real::rin_sten
    ineedt = irecexr(1)%tot
    i = iconsi
    ikg2 = 0
    do ismp = 1, typesten
      ikg = 0
      if ((ees .ne. 5) .or. (ismp .eq. 1)) then
        itarget = ielem(n, i)%inumneighbours
      else
        itarget = numneighbours2
      end if
      do l = 1, itarget
        if (ilocalstencil(n, i, ismp, l) .gt. 0) then
          ikg = ikg + 1
        end if
      end do
      if (ikg .eq. itarget) then
        if (ilocal_recon3(i)%local .eq. 1) then
          ikg2 = ikg2 + 1
          do l = 1, itarget
            ilox_ihexg(ikg2, l) = ilocalstencil(n, i, ismp, l)
            j = (xmpil(ilox_ihexg(ikg2, l)))
            ilox_ihexl(ikg2, l) = j
            ilox_ishape(ikg2, l) = ielem(n, j)%ishape
            call compute_centre2d(j, cords)
            ilox_xxc(ikg2, l) = cords(1)
            ilox_yyc(ikg2, l) = cords(2)
            ilon_nodcount(ikg2, l, 1:ielem(n, j)%nonodes) = ielem(n, j)%nodes(1:ielem(n, j)%nonodes)
            do k = 1, ielem(n, j)%nonodes
              ilon_x(ikg2, l, k) = inoder(ielem(n, j)%nodes(k))%cord(1)
              ilon_y(ikg2, l, k) = inoder(ielem(n, j)%nodes(k))%cord(2)
            end do
          end do
        else
          ikg2 = ikg2 + 1
          do l = 1, itarget
            ilox_ihexg(ikg2, l) = ilocalstencil(n, i, ismp, l)
            ilox_ihexb(ikg2, l) = xmpie(ilocalstencil(n, i, ismp, l))
            if (ilox_ihexb(ikg2, l) .eq. n) then
              j = (xmpil(ilox_ihexg(ikg2, l)))
              if (ilox_ihexg(ikg2, l) .eq. ielem(n, j)%ihexgl) then
                ilox_ihexl(ikg2, l) = j
                ilox_ishape(ikg2, l) = ielem(n, j)%ishape
              end if
              call compute_centre2d(j, cords)
              ilox_xxc(ikg2, l) = cords(1)
              ilox_yyc(ikg2, l) = cords(2)
              ilon_nodcount(ikg2, l, 1:ielem(n, j)%nonodes) = ielem(n, j)%nodes(1:ielem(n, j)%nonodes)
              do k = 1, ielem(n, j)%nonodes
                ilon_x(ikg2, l, k) = inoder(ielem(n, j)%nodes(k))%cord(1)
                ilon_y(ikg2, l, k) = inoder(ielem(n, j)%nodes(k))%cord(2)
              end do
            else
              do ixff = 1, ineedt
                if (iexcordr(ixff)%procid .eq. ilox_ihexb(ikg2, l)) then
                  do ixsst = 1, irecexr(ixff)%muchineed(1)
                    if ((ilox_ihexg(ikg2, l)) .eq. irecexr1(ixff)%whatineed(ixsst)) then
                      ilox_ihexl(ikg2, l) = ixsst
                      ilox_ihexn(ikg2, l) = ixff
                      ilox_ishape(ikg2, l) = irecexr1(ixff)%ishape(ixsst)
                      exit
                    end if
                  end do
                end if
              end do
              select case (ilox_ishape(ikg2, l))
              case (5)
                in_sten = 4
                rin_sten = 4.d0
              case (6)
                in_sten = 3
                rin_sten = 3.d0
              end select
              ilon_x(ikg2, l, 1:in_sten) = iexcordr(ilox_ihexn(ikg2, l))%nodecord(ilox_ihexl(ikg2, l), 1:in_sten, 1)
              ilon_y(ikg2, l, 1:in_sten) = iexcordr(ilox_ihexn(ikg2, l))%nodecord(ilox_ihexl(ikg2, l), 1:in_sten, 2)
              ilox_xxc(ikg2, l) = sum(ilon_x(ikg2, l, 1:in_sten))/rin_sten
              ilox_yyc(ikg2, l) = sum(ilon_y(ikg2, l, 1:in_sten))/rin_sten
            end if
          end do
        end if
      end if
    end do
  end subroutine localise_stencil2d
subroutine  localise_sten2(n,iconsi,ilox_ihexg,ilox_ihexl,ilox_ihexb,ilox_ihexn,ilox_ishape,ilox_xxc,ilox_yyc,ilox_zzc,ilox_volume,ilox_periodicflag,ilon_nodcount,ilon_x,ilon_y,ilon_z)
!> @brief
!> this subroutine continues expressing all the stencil elements coordinates and volumes with respect to the considered cell
    implicit none
    integer, intent(in)::iconsi, n
 integer::i, j, k, l, kk, prk, jj, kmaxe, ineedt, jx2, jx, iivd, iivd3, facexx, ixxfff, in1, itarget, idum, eltype, n_node, elem_dec
    integer, allocatable, dimension(:, :), intent(inout)::ilox_ihexg  !global index of cells
    integer, allocatable, dimension(:, :), intent(inout)::ilox_ihexl  !local index of cells
    integer, allocatable, dimension(:, :), intent(inout)::ilox_ihexb  !cpu that that each cell belongs to
    integer, allocatable, dimension(:, :), intent(inout)::ilox_ihexn  !internal index from where to take the values from communicated messages
    integer, allocatable, dimension(:, :), intent(inout)::ilox_ishape !shape of each element
    real, allocatable, dimension(:, :), intent(inout)::ilox_xxc       !cell centre coordinates in x
    real, allocatable, dimension(:, :), intent(inout)::ilox_yyc       !cell centre coordinates in y
    real, allocatable, dimension(:, :), intent(inout)::ilox_zzc      !cell centre coordinates in z
    real, allocatable, dimension(:, :), intent(inout)::ilox_volume    !cell volume
    integer, allocatable, dimension(:, :), intent(inout)::ilox_periodicflag
    integer, allocatable, dimension(:, :, :), intent(inout)::ilon_nodcount  !number of nodes
    real, allocatable, dimension(:, :, :), intent(inout)::ilon_x           !coordinates of each node in x
    real, allocatable, dimension(:, :, :), intent(inout)::ilon_y           !coordinates of each node in x
    real, allocatable, dimension(:, :, :), intent(inout)::ilon_z           !coordinates of each node in x
    real, dimension(3)::tempcentres, temp_cg
    real::dumv1, dumv2, detjc, dist1, ma, mb, mc, md, me, mf, mg, mh, mi, mdd, tempxx
    real, dimension(1:dimensiona)::cords
    real, dimension(1:8, 1:dimensiona)::vext, nodes_list
    real, dimension(1:6, 1:4, 1:dimensiona)::elem_listd
    real, dimension(1:dimensiona, 1:numberofpoints)::qpoints
    real, dimension(1:numberofpoints)::wequa3d
    real, dimension(1:dimensiona, 1:dimensiona)::vva1
    real, dimension(1)::deta
    ineedt = irecexr(1)%tot
    i = iconsi
    if (iperiodicity .eq. 1) then
      do jj = 1, ielem(n, i)%admis
        if ((ees .ne. 5) .or. (jj .eq. 1)) then
          itarget = ielem(n, i)%inumneighbours
        else
          itarget = numneighbours2
        end if
        do j = 2, itarget
        if (per_rot .eq. 0) then
        if (abs(ilox_xxc(jj, j) - ilox_xxc(jj, 1)) .gt. xper*oo2) then
          ilox_xxc(jj, j) = ilox_xxc(jj, j) + (xper*sign(1.0, ilox_xxc(jj, 1) - xper*oo2))
          do kk = 1, 8
            ilon_x(jj, j, kk) = ilon_x(jj, j, kk) + (xper*sign(1.0, ilox_xxc(jj, 1) - xper*oo2))
          end do
        end if
        if (abs(ilox_yyc(jj, j) - ilox_yyc(jj, 1)) .gt. yper*oo2) then
          ilox_yyc(jj, j) = ilox_yyc(jj, j) + (yper*sign(1.0, ilox_yyc(jj, 1) - yper*oo2))
          do kk = 1, 8
            ilon_y(jj, j, kk) = ilon_y(jj, j, kk) + (yper*sign(1.0, ilox_yyc(jj, 1) - yper*oo2))
          end do
        end if
        if (abs(ilox_zzc(jj, j) - ilox_zzc(jj, 1)) .gt. zper*oo2) then
          ilox_zzc(jj, j) = ilox_zzc(jj, j) + (zper*sign(1.0, ilox_zzc(jj, 1) - zper*oo2))
          do kk = 1, 8
            ilon_z(jj, j, kk) = ilon_z(jj, j, kk) + (zper*sign(1.0, ilox_zzc(jj, 1) - zper*oo2))
          end do
        end if
        else
        if (ilox_periodicflag(jj, j) .eq. 2) then
          tempxx = ilox_xxc(jj, j)
          ilox_xxc(jj, j) = tempxx*cos(angle_per) - ilox_yyc(jj, j)*sin(angle_per)
          ilox_yyc(jj, j) = tempxx*sin(angle_per) + ilox_yyc(jj, j)*cos(angle_per)
          !write(3300+n,'(6es14.6,i5)'),ilox_xxc(jj,1),ilox_yyc(jj,1),ilox_zzc(jj,1),ilox_xxc(jj,j),ilox_yyc(jj,j),ilox_zzc(jj,j),ilox_periodicflag(jj,j)
          do kk = 1, 8
            tempxx = ilon_x(jj, j, kk)
            ilon_x(jj, j, kk) = tempxx*cos(angle_per) - ilon_y(jj, j, kk)*sin(angle_per)
            ilon_y(jj, j, kk) = tempxx*sin(angle_per) + ilon_y(jj, j, kk)*cos(angle_per)
          end do
        end if
        if (ilox_periodicflag(jj, j) .eq. 1) then
          tempxx = ilox_xxc(jj, j)
          ilox_xxc(jj, j) = tempxx*cos(-angle_per) - ilox_yyc(jj, j)*sin(-angle_per)
          ilox_yyc(jj, j) = tempxx*sin(-angle_per) + ilox_yyc(jj, j)*cos(-angle_per)
          !write(3300+n,'(6es14.6,i5)'),ilox_xxc(jj,1),ilox_yyc(jj,1),ilox_zzc(jj,1),ilox_xxc(jj,j),ilox_yyc(jj,j),ilox_zzc(jj,j),ilox_periodicflag(jj,j)
          do kk = 1, 8
            tempxx = ilon_x(jj, j, kk)
            ilon_x(jj, j, kk) = tempxx*cos(-angle_per) - ilon_y(jj, j, kk)*sin(-angle_per)
            ilon_y(jj, j, kk) = tempxx*sin(-angle_per) + ilon_y(jj, j, kk)*cos(-angle_per)
          end do
        end if
        end if
        end do
      end do
    end if
    vext = 0.0d0
    nodes_list = 0.0d0
    eltype = ielem(n, i)%ishape
    elem_dec = ielem(n, i)%vdec
    elem_listd = 0.0d0
    jx = ielem(n, i)%nonodes
    do k = 1, jx
      jx2 = ielem(n, i)%nodes(k)
      nodes_list(k, 1) = ilon_x(1, 1, k)
      nodes_list(k, 2) = ilon_y(1, 1, k)
      nodes_list(k, 3) = ilon_z(1, 1, k)
      vext(k, :) = nodes_list(k, :)
    end do
    call decompose3(n, eltype, nodes_list, elem_listd)
    select case (ielem(n, i)%ishape)
    case (1)
      if (ielem(n, i)%mode .eq. 0) then
        call quadraturehexa(n, igqrules, vext, qpoints, wequa3d)
        dumv1 = hexavolume(n, vext, qpoints, wequa3d)
        call compute_centre3d(i, cords)
        vext(1, 1:dims) = cords(1:dims)
      else
        vext(1:4, 1:3) = elem_listd(1, 1:4, 1:3)
        call computejacobians(n, vext, vva1, deta)
      end if
    case (2)
      vext(1:4, 1:3) = elem_listd(1, 1:4, 1:3)
      call computejacobians(n, vext, vva1, deta)
    case (3)
      vext(1:4, 1:3) = elem_listd(1, 1:4, 1:3)
      call computejacobians(n, vext, vva1, deta)
    case (4)
      if (ielem(n, i)%mode .eq. 0) then
        call quadratureprism(n, igqrules, vext, qpoints, wequa3d)
        dumv1 = prismvolume(n, vext, qpoints, wequa3d)
        vext(1, 1:3) = vext(6, 1:3)
      else
        vext(1:4, 1:3) = elem_listd(1, 1:4, 1:3)
        call computejacobians(n, vext, vva1, deta)
      end if
    end select
    temp_cg(1:3) = vext(1, 1:3)
    ilocal_recon3(i)%invccjac(1:3, 1:3) = vva1(1:3, 1:3)
    detjc = deta(1)
    ilocal_recon3(i)%vext_ref(1:3) = temp_cg(1:3)
    if ((poly .eq. 4) .or. (dg .eq. 1)) then
      ilocal_recon3(i)%invccjac(1:3, 1:3) = zero
      ilocal_recon3(i)%invccjac(1, 1) = 1.0d0
      ilocal_recon3(i)%invccjac(2, 2) = 1.0d0
      ilocal_recon3(i)%invccjac(3, 3) = 1.0d0
      detjc = 1.0
      temp_cg(1) = ielem(n, i)%xxc
      temp_cg(2) = ielem(n, i)%yyc
      temp_cg(3) = ielem(n, i)%zzc
      ilocal_recon3(i)%vext_ref(1:3) = temp_cg(1:3)
    end if
    do jj = 1, ielem(n, i)%admis
      if ((ees .ne. 5) .or. (jj .eq. 1)) then
        itarget = ielem(n, i)%inumneighbours
      else
        itarget = numneighbours2
      end if
      do l = 1, itarget
        if (fastest .ne. 1) then
          do kk = 1, 8
            tempcentres = zero
            tempcentres(1) = ilon_x(jj, l, kk)
            tempcentres(2) = ilon_y(jj, l, kk)
            tempcentres(3) = ilon_z(jj, l, kk)
            tempcentres(:) = matmul(ilocal_recon3(i)%invccjac(:, :), tempcentres(:) - temp_cg(:))
            ilon_x(jj, l, kk) = tempcentres(1)
            ilon_y(jj, l, kk) = tempcentres(2)
            ilon_z(jj, l, kk) = tempcentres(3)
          end do
        end if
      end do
    end do
    idum = 0
    if (ielem(n, i)%interior .eq. 1) then
      do j = 1, ielem(n, i)%ifca
      if (ielem(n, i)%ibounds(j) .gt. 0) then
        if ((ibound(n, ielem(n, i)%ibounds(j))%icode .eq. 4) .or. (ibound(n, ielem(n, i)%ibounds(j))%icode .eq. 99)) then
          idum = 1
        end if
      end if
      end do
    end if
    do jj = 1, ielem(n, i)%admis
      if ((ees .ne. 5) .or. (jj .eq. 1)) then
        itarget = ielem(n, i)%inumneighbours
      else
        itarget = numneighbours2
      end if
      do l = 1, itarget
        if (ilocal_recon3(i)%local .eq. 1) then
          ilox_volume(jj, l) = (ielem(n, ilox_ihexl(jj, l))%totvolume)/abs(detjc)
          if ((ees .ne. 5) .or. (jj .eq. 1)) then
          if (idum .eq. 1) then
            ilocal_recon3(i)%volume(1, l) = ilox_volume(1, l)
          else
            ilocal_recon3(i)%volume(1, 1) = ilox_volume(1, 1)
          end if
          ilocal_recon3(i)%ihexg(jj, l) = ilox_ihexg(jj, l)
          ilocal_recon3(i)%ihexl(jj, l) = ilox_ihexl(jj, l)
          ilocal_recon3(i)%periodicflag(jj, l) = ilox_periodicflag(jj, l)
          else
          ilocal_recon3(i)%ihexgc(jj, l) = ilox_ihexg(jj, l)
          ilocal_recon3(i)%ihexlc(jj, l) = ilox_ihexl(jj, l)
          end if
        else
          if (ilox_ihexb(jj, l) .eq. n) then
            ilox_volume(jj, l) = (ielem(n, ilox_ihexl(jj, l))%totvolume)/abs(detjc)
            if ((ees .ne. 5) .or. (jj .eq. 1)) then
            if (idum .eq. 1) then
              ilocal_recon3(i)%volume(1, l) = ilox_volume(1, l)
            else
              ilocal_recon3(i)%volume(1, 1) = ilox_volume(1, 1)
            end if
            ilocal_recon3(i)%ihexg(jj, l) = ilox_ihexg(jj, l)
            ilocal_recon3(i)%ihexl(jj, l) = ilox_ihexl(jj, l)
            ilocal_recon3(i)%periodicflag(jj, l) = ilox_periodicflag(jj, l)
            ilocal_recon3(i)%ihexb(jj, l) = ilox_ihexb(jj, l)
            else
            ilocal_recon3(i)%volume(1, 1) = ilox_volume(1, 1)
            ilocal_recon3(i)%ihexgc(jj, l) = ilox_ihexg(jj, l)
            ilocal_recon3(i)%ihexlc(jj, l) = ilox_ihexl(jj, l)
            ilocal_recon3(i)%ihexbc(jj, l) = ilox_ihexb(jj, l)
            end if
          else
            if ((ees .ne. 5) .or. (jj .eq. 1)) then
              ilocal_recon3(i)%ihexg(jj, l) = ilox_ihexg(jj, l)
              ilocal_recon3(i)%ihexl(jj, l) = ilox_ihexl(jj, l)
              ilocal_recon3(i)%periodicflag(jj, l) = ilox_periodicflag(jj, l)
              ilocal_recon3(i)%ihexb(jj, l) = ilox_ihexb(jj, l)
              ilocal_recon3(i)%ihexn(jj, l) = ilox_ihexn(jj, l)
            else
              ilocal_recon3(i)%ihexgc(jj, l) = ilox_ihexg(jj, l)
              ilocal_recon3(i)%ihexlc(jj, l) = ilox_ihexl(jj, l)
              ilocal_recon3(i)%ihexbc(jj, l) = ilox_ihexb(jj, l)
              ilocal_recon3(i)%ihexnc(jj, l) = ilox_ihexn(jj, l)
            end if
            eltype = ilox_ishape(jj, l)
            elem_listd = 0.0d0; vext = 0.0d0; nodes_list = 0.0d0

            select case (eltype)
            case (1)
              elem_dec = 6; jx = 8
            case (2)
              elem_dec = 1; jx = 4
            case (3)
              elem_dec = 2; jx = 5
            case (4)
              elem_dec = 3; jx = 6
            end select
            do k = 1, jx
              nodes_list(k, 1) = ilon_x(jj, l, k)
              nodes_list(k, 2) = ilon_y(jj, l, k)
              nodes_list(k, 3) = ilon_z(jj, l, k)
              vext(k, :) = nodes_list(k, :)
            end do
            call decompose3(n, eltype, nodes_list, elem_listd)
            dumv2 = 0.0d0
            do k = 1, elem_dec
              vext(1:4, 1:3) = elem_listd(k, 1:4, 1:3)
              dumv2 = dumv2 + tetravolume(n, vext)
            end do
            ilox_volume(jj, l) = dumv2
            if ((ees .ne. 5) .or. (jj .eq. 1)) then
            if (idum .eq. 1) then
              ilocal_recon3(i)%volume(1, l) = ilox_volume(1, l)
            else
              ilocal_recon3(i)%volume(1, 1) = ilox_volume(1, 1)
            end if
            !ilocal_recon3(i)%volume(jj,l)=ilox_volume(jj,l)
            else
            ilocal_recon3(i)%volume(1, 1) = ilox_volume(1, 1)
            !ilocal_recon3(i)%volumec(jj,l)=ilox_volume(jj,l)
            end if
          end if
        end if
      end do
    end do
    if (ielem(n, i)%interior .eq. 0) then
      call compute_centre3d(i, cords)
      vext(1, 1:dims) = cords(1:dims)
      do k = 1, ielem(n, i)%ifca
        j = ielem(n, i)%ineigh(k)
        call compute_centre3d(j, cords)
        vext(2, 1:dims) = cords(1:dims)
        dist1 = distance3(n, vext)
        if (rungekutta .ge. 2) then
          ielem(n, i)%dih(k) = dist1
        end if
      end do
    else
      call compute_centre3d(i, cords)
      vext(1, 1:dims) = cords(1:dims)
      do k = 1, ielem(n, i)%ifca
        if (ielem(n, i)%ineighg(k) .eq. 0) then        !boundaries except other cpus and periodics
          facexx = k
          select case (ielem(n, i)%types_faces(k))
          case (5)
            ixxfff = 4
          case (6)
            ixxfff = 3
          end select
          call compute_centre3df(n, i, k, ixxfff, cords)
          vext(2, 1:dims) = cords(1:dims)
          dist1 = distance3(n, vext)
          if (rungekutta .ge. 2) then
            ielem(n, i)%dih(k) = dist1*2.0d0
          end if
        end if
        if ((ielem(n, i)%ineighg(k) .gt. 0) .and. (ielem(n, i)%ibounds(k) .eq. 0)) then        !non periodic boundaries
        if (ielem(n, i)%ineighb(k) .eq. n) then                !within my cpu
          j = ielem(n, i)%ineigh(k)
          call compute_centre3d(j, cords)
          vext(2, 1:dims) = cords(1:dims)
          dist1 = distance3(n, vext)
          if (rungekutta .ge. 2) then
            ielem(n, i)%dih(k) = dist1
          end if
        else                                                !from another cpu
          do in1 = 1, ielem(n, i)%inumneighbours
            if (ielem(n, i)%ineighg(k) .eq. ilocal_recon3(i)%ihexg(1, in1)) then
              ielem(n, i)%indexi(k) = in1
              if (rungekutta .ge. 2) then
                vext(2, 1) = ilox_xxc(1, in1); vext(2, 2) = ilox_yyc(1, in1); vext(2, 3) = ilox_zzc(1, in1)
                dist1 = distance3(n, vext)
                ielem(n, i)%dih(k) = dist1
              end if
            end if
          end do
        end if
        end if
        if ((ielem(n, i)%ineighg(k) .gt. 0) .and. (ielem(n, i)%ibounds(k) .gt. 0)) then        !periodic boundaries within my cpu
        if (ielem(n, i)%ineighb(k) .eq. n) then
          j = ielem(n, i)%ineigh(k)
          call compute_centre3d(j, cords)
          vext(2, 1:dims) = cords(1:dims)
          if (per_rot .eq. 0) then
          if (abs(vext(2, 1) - vext(1, 1)) .gt. xper*oo2) then
            vext(2, 1) = vext(2, 1) + (xper*sign(1.0, vext(1, 1) - xper*oo2))
          end if
          if (abs(vext(2, 2) - vext(1, 2)) .gt. yper*oo2) then
            vext(2, 2) = vext(2, 2) + (yper*sign(1.0, vext(1, 2) - yper*oo2))
          end if
          if (abs(vext(2, 3) - vext(1, 3)) .gt. zper*oo2) then
            vext(2, 3) = vext(2, 3) + (zper*sign(1.0, vext(1, 3) - zper*oo2))
          end if
          else
          if (ibound(n, ielem(n, i)%ibounds(k))%icode .eq. 5) then
            do kk = 1, n_node
              tempxx = vext(kk, 1)
              vext(kk, 1) = tempxx*cos(-angle_per) - sin(-angle_per)*vext(kk, 2)
              vext(kk, 2) = tempxx*sin(-angle_per) + cos(-angle_per)*vext(kk, 2)
            end do
          else
            do kk = 1, n_node
              tempxx = vext(kk, 1)
              vext(kk, 1) = tempxx*cos(angle_per) - sin(angle_per)*vext(kk, 2)
              vext(kk, 2) = tempxx*sin(angle_per) + cos(angle_per)*vext(kk, 2)
            end do
          end if
          end if
          dist1 = distance3(n, vext)
          if (rungekutta .ge. 2) then
            ielem(n, i)%dih(k) = dist1
          end if
        else        !periodic boundaries from another cpu
          do in1 = 1, ielem(n, i)%inumneighbours
            if (ielem(n, i)%ineighg(k) .eq. ilocal_recon3(i)%ihexg(1, in1)) then
              ielem(n, i)%indexi(k) = in1
              if (rungekutta .ge. 2) then
                vext(2, 1) = ilox_xxc(1, in1); vext(2, 2) = ilox_yyc(1, in1); vext(2, 3) = ilox_zzc(1, in1)
                if (per_rot .eq. 0) then
                if (abs(vext(2, 1) - vext(1, 1)) .gt. xper*oo2) then
                  vext(2, 1) = vext(2, 1) + (xper*sign(1.0, vext(1, 1) - xper*oo2))
                end if
                if (abs(vext(2, 2) - vext(1, 2)) .gt. yper*oo2) then
                  vext(2, 2) = vext(2, 2) + (yper*sign(1.0, vext(1, 2) - yper*oo2))
                end if
                if (abs(vext(2, 3) - vext(1, 3)) .gt. zper*oo2) then
                  vext(2, 3) = vext(2, 3) + (zper*sign(1.0, vext(1, 3) - zper*oo2))
                end if
                else
                if (ibound(n, ielem(n, i)%ibounds(k))%icode .eq. 5) then
                  do kk = 1, n_node
                    tempxx = vext(kk, 1)
                    vext(kk, 1) = tempxx*cos(-angle_per) - sin(-angle_per)*vext(kk, 2)
                    vext(kk, 2) = tempxx*sin(-angle_per) + cos(-angle_per)*vext(kk, 2)
                  end do
                else
                  do kk = 1, n_node
                    tempxx = vext(kk, 1)
                    vext(kk, 1) = tempxx*cos(angle_per) - sin(angle_per)*vext(kk, 2)
                    vext(kk, 2) = tempxx*sin(angle_per) + cos(angle_per)*vext(kk, 2)
                  end do
                end if
                end if
                dist1 = distance3(n, vext)
                ielem(n, i)%dih(k) = dist1
              end if
            end if
          end do

        end if
        end if
      end do
    end if
    do jj = 1, ielem(n, i)%admis
      if ((ees .ne. 5) .or. (jj .eq. 1)) then
        itarget = ielem(n, i)%inumneighbours
      else
        itarget = numneighbours2
      end if
      do l = 1, itarget
        if (fastest .ne. 1) then
          tempcentres = zero
          tempcentres(1) = ilox_xxc(jj, l)
          tempcentres(2) = ilox_yyc(jj, l)
          tempcentres(3) = ilox_zzc(jj, l)
          tempcentres(:) = matmul(ilocal_recon3(i)%invccjac(:, :), tempcentres(:) - temp_cg(:))
          ilox_xxc(jj, l) = tempcentres(1)
          ilox_yyc(jj, l) = tempcentres(2)
          ilox_zzc(jj, l) = tempcentres(3)
        end if
      end do
    end do
  end subroutine localise_sten2

subroutine  localise_sten2d(n,iconsi,ilox_ihexg,ilox_ihexl,ilox_ihexb,ilox_ihexn,ilox_ishape,ilox_xxc,ilox_yyc,ilox_zzc,ilox_volume,ilox_periodicflag,ilon_nodcount,ilon_x,ilon_y,ilon_z)
!> @brief
!> this subroutine continues expressing all the stencil elements coordinates and volumes with respect to the considered cell in 2d
    implicit none
    integer::i,j,k,l,kk,prk,jj,kmaxe,ineedt,jx2,jx,in1,facexx,ixxfff,ihgt,ihgj,itarget,idum,in_sten,nj,elem_dec,eltype
    integer, allocatable, dimension(:, :), intent(inout)::ilox_ihexg  !global index of cells
    integer, allocatable, dimension(:, :), intent(inout)::ilox_ihexl  !local index of cells
    integer, allocatable, dimension(:, :), intent(inout)::ilox_ihexb  !cpu that that each cell belongs to
    integer, allocatable, dimension(:, :), intent(inout)::ilox_ihexn  !internal index from where to take the values from communicated messages
    integer, allocatable, dimension(:, :), intent(inout)::ilox_ishape !shape of each element
    real, allocatable, dimension(:, :), intent(inout)::ilox_xxc       !cell centre coordinates in x
    real, allocatable, dimension(:, :), intent(inout)::ilox_yyc       !cell centre coordinates in y
    real, allocatable, dimension(:, :), intent(inout)::ilox_zzc      !cell centre coordinates in z
    real, allocatable, dimension(:, :), intent(inout)::ilox_volume    !cell volume
    integer, allocatable, dimension(:, :), intent(inout)::ilox_periodicflag
    integer, allocatable, dimension(:, :, :), intent(inout)::ilon_nodcount  !number of nodes
    real, allocatable, dimension(:, :, :), intent(inout)::ilon_x           !coordinates of each node in x
    real, allocatable, dimension(:, :, :), intent(inout)::ilon_y           !coordinates of each node in x
    real, allocatable, dimension(:, :, :), intent(inout)::ilon_z           !coordinates of each node in x
    real, dimension(2)::tempcentres, rel2
    real::dumv1, dumv2, detjc, dist1, distfd, x1x, y1y, x2x, y2y
    real, dimension(1:dimensiona)::cords
    integer, intent(in)::iconsi, n
    integer, dimension(4)::nojcount
    real, dimension(1:8, 1:dimensiona)::vext, nodes_list
    real, dimension(1:6, 1:4, 1:dimensiona)::elem_listd
    real, dimension(1:dimensiona, 1:numberofpoints)::qpoints
    real, dimension(1:numberofpoints)::wequa3d
    real, dimension(1:dimensiona, 1:dimensiona)::vva1
    real, dimension(1)::deta
    kmaxe = xmpielrank(n)
    ineedt = irecexr(1)%tot
    i = iconsi
    if (iperiodicity .eq. 1) then
      do jj = 1, ielem(n, i)%admis
      if ((ees .ne. 5) .or. (jj .eq. 1)) then
        itarget = ielem(n, i)%inumneighbours
      else
        itarget = numneighbours
      end if
      do j = 2, itarget
      if (abs(ilox_xxc(jj, j) - ilox_xxc(jj, 1)) .gt. xper*oo2) then
        ilox_xxc(jj, j) = ilox_xxc(jj, j) + (xper*sign(1.0, ilox_xxc(jj, 1) - xper*oo2))
        do kk = 1, 4
          ilon_x(jj, j, kk) = ilon_x(jj, j, kk) + (xper*sign(1.0, ilox_xxc(jj, 1) - xper*oo2))
        end do
      end if
      if (abs(ilox_yyc(jj, j) - ilox_yyc(jj, 1)) .gt. yper*oo2) then
        ilox_yyc(jj, j) = ilox_yyc(jj, j) + (yper*sign(1.0, ilox_yyc(jj, 1) - yper*oo2))
        do kk = 1, 4
          ilon_y(jj, j, kk) = ilon_y(jj, j, kk) + (yper*sign(1.0, ilox_yyc(jj, 1) - yper*oo2))
        end do
      end if
      end do
      end do
    end if
    if (code_profile .eq. 30) then
      do nj = 1, ielem(n, i)%nonodes
        x1x = ilon_x(1, 1, nj)
        y1y = ilon_y(1, 1, nj)
        nojcount(nj) = 0
        if (ilocal_recon3(i)%local .eq. 1) then
          do l = 2, itarget
            j = (xmpil(ilox_ihexg(1, l)))
            do k = 1, ielem(n, j)%nonodes
              x2x = ilon_x(1, l, k)
              y2y = ilon_y(1, l, k)
              distfd = sqrt(((x1x - x2x)**2) + ((y1y - y2y)**2))
              if (distfd .lt. tolsmall) then
                nojcount(nj) = nojcount(nj) + 1
                ielem(n, i)%nodes_neighbours(nj, nojcount(nj)) = l
              end if
            end do
          end do
        end if
        !---------------------- mixed----------!
        if (ilocal_recon3(i)%local .ne. 1) then
          do l = 2, itarget
            if (ilox_ihexb(1, l) .eq. n) then
              j = (xmpil(ilox_ihexg(1, l)))
              do k = 1, ielem(n, j)%nonodes
                x2x = ilon_x(1, l, k)
                y2y = ilon_y(1, l, k)
                distfd = sqrt(((x1x - x2x)**2) + ((y1y - y2y)**2))
                if (distfd .lt. tolsmall) then
                  nojcount(nj) = nojcount(nj) + 1
                  ielem(n, i)%nodes_neighbours(nj, nojcount(nj)) = l
                end if
              end do
            else
              select case (ilox_ishape(1, l))
              case (5)
                in_sten = 4
              case (6)
                in_sten = 3
              end select
              do k = 1, in_sten
                x2x = ilon_x(1, l, k)
                y2y = ilon_y(1, l, k)
                distfd = sqrt(((x1x - x2x)**2) + ((y1y - y2y)**2))
                if (distfd .lt. tolsmall) then
                  nojcount(nj) = nojcount(nj) + 1
                  ielem(n, i)%nodes_neighbours(nj, nojcount(nj)) = l
                end if
              end do
            end if
          end do
        end if
      end do
      write (630 + n, *) "element number global", ielem(n, i)%ihexgl
      allocate (ielem(n, i)%nojecount(ielem(n, i)%nonodes))
      do nj = 1, ielem(n, i)%nonodes
        write (630 + n, *) "node number", nj
        ielem(n, i)%nojecount(nj) = nojcount(nj)
        do j = 1, nojcount(nj)
          write (630 + n, *) j, ielem(n, i)%nodes_neighbours(nj, j)
          !end if
        end do
      end do
    end if
    vext = 0.0d0
    nodes_list = 0.0d0
    eltype = ielem(n, i)%ishape
    elem_dec = ielem(n, i)%vdec
    elem_listd = 0.0d0
    jx = ielem(n, i)%nonodes
    do k = 1, jx
      jx2 = ielem(n, i)%nodes(k)
      nodes_list(k, 1) = ilon_x(1, 1, k)
      nodes_list(k, 2) = ilon_y(1, 1, k)
      vext(k, :) = nodes_list(k, :)
    end do
    call decompose2(n, eltype, nodes_list, elem_listd)
    select case (ielem(n, i)%ishape)
    case (5)
      if (ielem(n, i)%mode .eq. 0) then
        call quadraturequad(n, igqrules, vext, qpoints, wequa3d)
        dumv1 = quadvolume(n, vext, qpoints, wequa3d)
        call compute_centre2d(i, cords)
        vext(1, 1:dims) = cords(1:dims)
      else
        vext(1:3, 1:2) = elem_listd(1, 1:3, 1:2)
        call computejacobians2(n, vext, vva1, deta)
      end if
    case (6)
      vext(1:3, 1:2) = elem_listd(1, 1:3, 1:2)
      call computejacobians2(n, vext, vva1, deta)
    end select
    ilocal_recon3(i)%invccjac(1:2, 1:2) = vva1(1:2, 1:2)
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
    detjc = deta(1)
    ilocal_recon3(i)%vext_ref(1:2) = vext(1, 1:2)
    if ((poly .eq. 4) .or. (dg .eq. 1)) then
      ilocal_recon3(i)%invccjac(1:2, 1:2) = 0.0d0
      ilocal_recon3(i)%invccjac(1, 1) = 1.0d0
      ilocal_recon3(i)%invccjac(2, 2) = 1.0d0
      detjc = 1.0
      vext(1, 1) = ielem(n, i)%xxc
      vext(1, 2) = ielem(n, i)%yyc
      ilocal_recon3(i)%vext_ref(1:2) = vext(1, 1:2)
    end if
    do jj = 1, ielem(n, i)%admis
      if ((ees .ne. 5) .or. (jj .eq. 1)) then
        itarget = ielem(n, i)%inumneighbours
      else
        itarget = numneighbours2
      end if
      do l = 1, itarget
        do kk = 1, 4
          tempcentres = zero
          tempcentres(1) = ilon_x(jj, l, kk)
          tempcentres(2) = ilon_y(jj, l, kk)
          tempcentres(:) = matmul(ilocal_recon3(i)%invccjac(:, :), tempcentres(:) - vext(1, :))
          ilon_x(jj, l, kk) = tempcentres(1)
          ilon_y(jj, l, kk) = tempcentres(2)
        end do
      end do
    end do
    do jj = 1, ielem(n, i)%admis
      if ((ees .ne. 5) .or. (jj .eq. 1)) then
        itarget = ielem(n, i)%inumneighbours
      else
        itarget = numneighbours2
      end if
      do l = 1, itarget
        if (ilocal_recon3(i)%local .eq. 1) then
          ilox_volume(jj, l) = (ielem(n, ilox_ihexl(jj, l))%totvolume)/abs(detjc)
          if ((ees .ne. 5) .or. (jj .eq. 1)) then
          if (idum .eq. 1) then
            ilocal_recon3(i)%volume(1, l) = ilox_volume(1, l)
          else
            ilocal_recon3(i)%volume(1, 1) = ilox_volume(1, 1)
          end if
          ilocal_recon3(i)%ihexg(jj, l) = ilox_ihexg(jj, l)
          ilocal_recon3(i)%ihexl(jj, l) = ilox_ihexl(jj, l)
          else
          ilocal_recon3(i)%volume(1, 1) = ilox_volume(1, 1)
          ilocal_recon3(i)%ihexgc(jj, l) = ilox_ihexg(jj, l)
          ilocal_recon3(i)%ihexlc(jj, l) = ilox_ihexl(jj, l)
          end if
        else
          if (ilox_ihexb(jj, l) .eq. n) then
            ilox_volume(jj, l) = (ielem(n, ilox_ihexl(jj, l))%totvolume)/abs(detjc)
            if ((ees .ne. 5) .or. (jj .eq. 1)) then
              if (idum .eq. 1) then
                ilocal_recon3(i)%volume(1, l) = ilox_volume(1, l)
              else
                ilocal_recon3(i)%volume(1, 1) = ilox_volume(1, 1)
              end if
              ilocal_recon3(i)%ihexg(jj, l) = ilox_ihexg(jj, l)
              ilocal_recon3(i)%ihexl(jj, l) = ilox_ihexl(jj, l)
              ilocal_recon3(i)%ihexb(jj, l) = ilox_ihexb(jj, l)
            else
              ilocal_recon3(i)%volume(1, 1) = ilox_volume(1, 1)
              ilocal_recon3(i)%ihexgc(jj, l) = ilox_ihexg(jj, l)
              ilocal_recon3(i)%ihexlc(jj, l) = ilox_ihexl(jj, l)
              ilocal_recon3(i)%ihexbc(jj, l) = ilox_ihexb(jj, l)
            end if
          else
            if ((ees .ne. 5) .or. (jj .eq. 1)) then
              ilocal_recon3(i)%ihexg(jj, l) = ilox_ihexg(jj, l)
              ilocal_recon3(i)%ihexl(jj, l) = ilox_ihexl(jj, l)
              ilocal_recon3(i)%ihexb(jj, l) = ilox_ihexb(jj, l)
              ilocal_recon3(i)%ihexn(jj, l) = ilox_ihexn(jj, l)
            else
              ilocal_recon3(i)%ihexgc(jj, l) = ilox_ihexg(jj, l)
              ilocal_recon3(i)%ihexlc(jj, l) = ilox_ihexl(jj, l)
              ilocal_recon3(i)%ihexbc(jj, l) = ilox_ihexb(jj, l)
              ilocal_recon3(i)%ihexnc(jj, l) = ilox_ihexn(jj, l)
            end if
            eltype = ilox_ishape(jj, l)
            elem_listd = 0.0d0; vext = 0.0d0; nodes_list = 0.0d0
            select case (eltype)
            case (5)
              elem_dec = 2; jx = 4
            case (6)
              elem_dec = 1; jx = 3
            end select
            do k = 1, jx
              nodes_list(k, 1) = ilon_x(jj, l, k)
              nodes_list(k, 2) = ilon_y(jj, l, k)
              vext(k, :) = nodes_list(k, :)
            end do
            call decompose2(n, eltype, nodes_list, elem_listd)
            dumv2 = 0.0d0
            do k = 1, elem_dec
              vext(1:3, 1:2) = elem_listd(k, 1:3, 1:2)
              dumv2 = dumv2 + trianglevolume(n, vext)
            end do
            ilox_volume(jj, l) = dumv2
            if ((ees .ne. 5) .or. (jj .eq. 1)) then
              if (idum .eq. 1) then
                ilocal_recon3(i)%volume(1, l) = ilox_volume(1, l)
              else
                ilocal_recon3(i)%volume(1, 1) = ilox_volume(1, 1)
              end if
            else
              ilocal_recon3(i)%volume(1, 1) = ilox_volume(1, 1)
            end if
          end if
        end if
      end do
    end do
    if (ielem(n, i)%interior .eq. 0) then
      call compute_centre2d(i, cords)
      vext(1, 1:dims) = cords(1:dims)
      do k = 1, ielem(n, i)%ifca
        j = ielem(n, i)%ineigh(k)
        call compute_centre2d(j, cords)
        vext(2, 1:dims) = cords(1:dims)
        dist1 = distance2(n, vext)
        if (rungekutta .ge. 2) then
          ielem(n, i)%dih(k) = dist1
        end if
      end do
    else
      call compute_centre2d(i, cords)
      vext(1, 1:dims) = cords(1:dims)
      do k = 1, ielem(n, i)%ifca
        if (ielem(n, i)%ineighg(k) .eq. 0) then        !boundaries except other cpus and periodics
          facexx = k
          ixxfff = 2
          call compute_centre2df(n, iconsi, facexx, ixxfff, cords)
          vext(2, 1:dims) = cords(1:dims)
          dist1 = distance2(n, vext)
          if (rungekutta .ge. 2) then
            ielem(n, i)%dih(k) = dist1*2.0d0
          end if
        end if
        if ((ielem(n, i)%ineighg(k) .gt. 0) .and. (ielem(n, i)%ibounds(k) .eq. 0)) then        !non periodic boundaries
        if (ielem(n, i)%ineighb(k) .eq. n) then                !within my cpu
          j = ielem(n, i)%ineigh(k)
          call compute_centre2d(j, cords)
          vext(2, 1:dims) = cords(1:dims)
          dist1 = distance2(n, vext)
          if (rungekutta .ge. 2) then
            ielem(n, i)%dih(k) = dist1
          end if
        else                                                !from another cpu
          do in1 = 1, ielem(n, i)%inumneighbours
            if (ielem(n, i)%ineighg(k) .eq. ilocal_recon3(i)%ihexg(1, in1)) then
              ielem(n, i)%indexi(k) = in1
              if (rungekutta .ge. 2) then
                vext(2, 1) = ilox_xxc(1, in1); vext(2, 2) = ilox_yyc(1, in1)
                dist1 = distance2(n, vext)
                ielem(n, i)%dih(k) = dist1
              end if
            end if
          end do
        end if
        end if
        if ((ielem(n, i)%ineighg(k) .gt. 0) .and. (ielem(n, i)%ibounds(k) .gt. 0)) then        !periodic boundaries within my cpu
        if (ielem(n, i)%ineighb(k) .eq. n) then
          j = ielem(n, i)%ineigh(k)
          call compute_centre2d(j, cords)
          vext(2, 1:dims) = cords(1:dims)
          if (abs(vext(2, 1) - vext(1, 1)) .gt. xper*oo2) then
            vext(2, 1) = vext(2, 1) + (xper*sign(1.0, vext(1, 1) - xper*oo2))
          end if
          if (abs(vext(2, 2) - vext(1, 2)) .gt. yper*oo2) then
            vext(2, 2) = vext(2, 2) + (yper*sign(1.0, vext(1, 2) - yper*oo2))
          end if
          dist1 = distance2(n, vext)
          if (rungekutta .ge. 2) then
            ielem(n, i)%dih(k) = dist1
          end if
        else        !periodic boundaries from another cpu
          do in1 = 1, ielem(n, i)%inumneighbours
            if (ielem(n, i)%ineighg(k) .eq. ilocal_recon3(i)%ihexg(1, in1)) then
              ielem(n, i)%indexi(k) = in1
              if (rungekutta .ge. 2) then
                vext(2, 1) = ilox_xxc(1, in1); vext(2, 2) = ilox_yyc(1, in1); 
                if (abs(vext(2, 1) - vext(1, 1)) .gt. xper*oo2) then
                  vext(2, 1) = vext(2, 1) + (xper*sign(1.0, vext(1, 1) - xper*oo2))
                end if
                if (abs(vext(2, 2) - vext(1, 2)) .gt. yper*oo2) then
                  vext(2, 2) = vext(2, 2) + (yper*sign(1.0, vext(1, 2) - yper*oo2))
                end if
                dist1 = distance2(n, vext)
                ielem(n, i)%dih(k) = dist1
              end if
            end if
          end do
        end if
        end if
      end do
    end if
    do jj = 1, ielem(n, i)%admis
      if ((ees .ne. 5) .or. (jj .eq. 1)) then
        itarget = ielem(n, i)%inumneighbours
      else
        itarget = numneighbours2
      end if
      do l = 1, itarget
        tempcentres = zero
        tempcentres(1) = ilox_xxc(jj, l)
        tempcentres(2) = ilox_yyc(jj, l)
        tempcentres(:) = matmul(ilocal_recon3(i)%invccjac(:, :), tempcentres(:) - vext(1, :))
        ilox_xxc(jj, l) = tempcentres(1)
        ilox_yyc(jj, l) = tempcentres(2)
      end do
    end do
  end subroutine localise_sten2d
  subroutine direct_side(n)
!> @brief
!> this subroutine establishes the distance betwen cell centres for each face
    implicit none
    integer, intent(in)::n
    integer::i, j, k, kmaxe, facexx, ixxfff
    real, dimension(1:8, 1:dimensiona)::vext
    real, dimension(1:dimensiona)::cords
    real::dist1
    kmaxe = xmpielrank(n)
!$omp do
    do i = 1, kmaxe
      if (ielem(n, i)%interior .eq. 0) then
        call compute_centre3d(i, cords)
        vext(1, 1:dims) = cords(1:dims)
        do k = 1, ielem(n, i)%ifca
          j = ielem(n, i)%ineigh(k)
          call compute_centre3d(j, cords)
          vext(2, 1:dims) = cords(1:dims)
          dist1 = distance3(n, vext)
          ielem(n, i)%dih(k) = dist1
        end do
      else
        call compute_centre3d(i, cords)
        vext(1, 1:dims) = cords(1:dims)
        do k = 1, ielem(n, i)%ifca
          if (ielem(n, i)%ineighg(k) .eq. 0) then        !boundaries except other cpus and periodics
            facexx = k
            select case (ielem(n, i)%types_faces(k))
            case (5)
              ixxfff = 4
            case (6)
              ixxfff = 3
            end select
            call compute_centre3df(n, i, k, ixxfff, cords)
            vext(2, 1:dims) = cords(1:dims)
            dist1 = distance3(n, vext)
            ielem(n, i)%dih(k) = dist1*2.0d0
          end if
          if ((ielem(n, i)%ineighg(k) .gt. 0) .and. (ielem(n, i)%ibounds(k) .eq. 0)) then        !non periodic boundaries
          if (ielem(n, i)%ineighb(k) .eq. n) then                !within my cpu
            j = ielem(n, i)%ineigh(k)
            call compute_centre3d(j, cords)
            vext(2, 1:dims) = cords(1:dims)
            dist1 = distance3(n, vext)
            ielem(n, i)%dih(k) = dist1
          else                                                !from another cpu
            vext(2, 1:dims) = solchanger(ielem(n, i)%ineighn(k))%centres(ielem(n, i)%q_face(k)%q_mapl(1), 1:dims)
            dist1 = distance3(n, vext)
            ielem(n, i)%dih(k) = dist1
          end if
          end if
          if ((ielem(n, i)%ineighg(k) .gt. 0) .and. (ielem(n, i)%ibounds(k) .gt. 0)) then        !periodic boundaries within my cpu
          if (ielem(n, i)%ineighb(k) .eq. n) then
            j = ielem(n, i)%ineigh(k)
            call compute_centre3d(j, cords)
            vext(2, 1:dims) = cords(1:dims)
            if (abs(vext(2, 1) - vext(1, 1)) .gt. xper*oo2) then
              vext(2, 1) = vext(2, 1) + (xper*sign(1.0, vext(1, 1) - xper*oo2))
            end if
            if (abs(vext(2, 2) - vext(1, 2)) .gt. yper*oo2) then
              vext(2, 2) = vext(2, 2) + (yper*sign(1.0, vext(1, 2) - yper*oo2))
            end if
            if (abs(vext(2, 3) - vext(1, 3)) .gt. zper*oo2) then
              vext(2, 3) = vext(2, 3) + (zper*sign(1.0, vext(1, 3) - zper*oo2))
            end if
            dist1 = distance3(n, vext)
            ielem(n, i)%dih(k) = dist1
          else        !periodic boundaries from another cpu
            vext(2, 1:dims) = solchanger(ielem(n, i)%ineighn(k))%centres(ielem(n, i)%q_face(k)%q_mapl(1), 1:dims)
            if (abs(vext(2, 1) - vext(1, 1)) .gt. xper*oo2) then
              vext(2, 1) = vext(2, 1) + (xper*sign(1.0, vext(1, 1) - xper*oo2))
            end if
            if (abs(vext(2, 2) - vext(1, 2)) .gt. yper*oo2) then
              vext(2, 2) = vext(2, 2) + (yper*sign(1.0, vext(1, 2) - yper*oo2))
            end if
            if (abs(vext(2, 3) - vext(1, 3)) .gt. zper*oo2) then
              vext(2, 3) = vext(2, 3) + (zper*sign(1.0, vext(1, 3) - zper*oo2))
            end if
            dist1 = distance3(n, vext)
            ielem(n, i)%dih(k) = dist1
          end if
          end if
        end do
      end if
    end do
  end subroutine direct_side
  subroutine direct_side2d(n)
!> @brief
!> this subroutine establishes the distance betwen cell centres for each edge
    implicit none
    integer, intent(in)::n
    integer::i, j, k, kmaxe, facexx, ixxfff
    real, dimension(1:8, 1:dimensiona)::vext
    real, dimension(1:dimensiona)::cords
    real::dist1
    kmaxe = xmpielrank(n)
    do i = 1, kmaxe
      if (ielem(n, i)%interior .eq. 0) then
        call compute_centre2d(i, cords)
        vext(1, 1:dims) = cords(1:dims)
        do k = 1, ielem(n, i)%ifca
          j = ielem(n, i)%ineigh(k)
          call compute_centre2d(j, cords)
          vext(2, 1:dims) = cords(1:dims)
          dist1 = distance2(n, vext)
          ielem(n, i)%dih(k) = dist1
        end do
      else
        call compute_centre2d(i, cords)
        vext(1, 1:dims) = cords(1:dims)
        do k = 1, ielem(n, i)%ifca
          if (ielem(n, i)%ineighg(k) .eq. 0) then        !boundaries except other cpus and periodics
            facexx = k
            ixxfff = 2
            call compute_centre2df(n, i, facexx, ixxfff, cords)
            vext(2, 1:dims) = cords(1:dims)
            dist1 = distance2(n, vext)
            ielem(n, i)%dih(k) = dist1*2.0d0
          end if
          if ((ielem(n, i)%ineighg(k) .gt. 0) .and. (ielem(n, i)%ibounds(k) .eq. 0)) then        !non periodic boundaries
          if (ielem(n, i)%ineighb(k) .eq. n) then                !within my cpu
            j = ielem(n, i)%ineigh(k)
            call compute_centre2d(j, cords)
            vext(2, 1:dims) = cords(1:dims)
            dist1 = distance2(n, vext)
            ielem(n, i)%dih(k) = dist1
          else                                                !from another cpu
            vext(2, 1:dims) = solchanger(ielem(n, i)%ineighn(k))%centres(ielem(n, i)%q_face(k)%q_mapl(1), 1:dims)
            dist1 = distance2(n, vext)
            ielem(n, i)%dih(k) = dist1
          end if
          end if
          if ((ielem(n, i)%ineighg(k) .gt. 0) .and. (ielem(n, i)%ibounds(k) .gt. 0)) then        !periodic boundaries within my cpu
          if (ielem(n, i)%ineighb(k) .eq. n) then
            j = ielem(n, i)%ineigh(k)
            call compute_centre2d(j, cords)
            vext(2, 1:dims) = cords(1:dims)
            if (abs(vext(2, 1) - vext(1, 1)) .gt. xper/2.d0) then
              vext(2, 1) = vext(2, 1) + (xper*sign(1.0d0, vext(1, 1) - xper/2.d0))
            end if
            if (abs(vext(2, 2) - vext(1, 2)) .gt. yper/2.d0) then
              vext(2, 2) = vext(2, 2) + (yper*sign(1.0d0, vext(1, 2) - yper/2.d0))
            end if
            dist1 = distance2(n, vext)
            ielem(n, i)%dih(k) = dist1
          else        !periodic boundaries from another cpu
            vext(2, 1:dims) = solchanger(ielem(n, i)%ineighn(k))%centres(ielem(n, i)%q_face(k)%q_mapl(1), 1:dims)
            if (abs(vext(2, 1) - vext(1, 1)) .gt. xper*oo2) then
              vext(2, 1) = vext(2, 1) + (xper*sign(1.0d0, vext(1, 1) - xper/2.0d0))
            end if
            if (abs(vext(2, 2) - vext(1, 2)) .gt. yper*oo2) then
              vext(2, 2) = vext(2, 2) + (yper*sign(1.0d0, vext(1, 2) - yper/2.0d0))
            end if
            dist1 = distance2(n, vext)
            ielem(n, i)%dih(k) = dist1
          end if
          end if
        end do
      end if
    end do
  end subroutine direct_side2d
  subroutine grads_assign(n)
    integer, intent(in)::n
    integer::i, j, k, l, jj, kmaxe
    kmaxe = xmpielrank(n)
    if (dimensiona .eq. 3) then
      do i = 1, kmaxe
        call checkgrads(n, i)
      end do
    else
      do i = 1, kmaxe
        call checkgrads2d(n, i)
      end do
    end if
  end subroutine grads_assign
  subroutine checkgrads(n, iconsi)
!> @brief
!> this subroutine assigns the correct viscous gradient approximation flag for each cell based on some additional geometrical characteristics
    implicit none
    integer, intent(in)::n, iconsi
    real::dxx1, dxx2, tempg1, dist1, dist2, oo2, surfmin, surfmax
    integer::i, j, k, l, jj, icount3, nnd, ixf4, idc, idc2
    real, dimension(1:dimensiona)::cords
    real, dimension(1:8, 1:dimensiona)::vext, nodes_list
    i = iconsi
    ielem(n, i)%ggs = greengo
    dxx1 = -tolbig; dxx2 = tolbig
    call compute_centre3d(i, cords)
    vext(1, 1:dims) = cords(1:dims)
    i = iconsi
    if (ielem(n, i)%interior .eq. 1) then
      do k = 1, ielem(n, i)%ifca
        if (ielem(n, i)%types_faces(k) .eq. 5) then
          ixf4 = 4
        else
          ixf4 = 3
        end if
        if (ielem(n, i)%ineighg(k) .gt. 0) then
          if (ielem(n, i)%ineighb(k) .ne. n) then
            if (ielem(n, i)%ineighg(k) .gt. ielem(n, i)%ihexgl) then
           ielem(n, i)%nodes_faces(k, 1:ixf4) = iexboundhirr(ielem(n, i)%ineighn(k))%vertpp(ielem(n, i)%q_face(k)%q_mapl(1), 1:ixf4)
              ielem(n, i)%reorient(k) = 1
            end if
          else
            if (ielem(n, i)%ineighg(k) .gt. ielem(n, i)%ihexgl) then
              ielem(n, i)%nodes_faces(k, 1:ixf4) = ielem(n, ielem(n, i)%ineigh(k))%nodes_faces(ielem(n, i)%ineighn(k), 1:ixf4)
              ielem(n, i)%reorient(k) = 1
            end if
          end if
        end if
      end do
    else
      do k = 1, ielem(n, i)%ifca
        if (ielem(n, i)%types_faces(k) .eq. 5) then
          ixf4 = 4
        else
          ixf4 = 3
        end if
        if (ielem(n, i)%ineighg(k) .gt. 0) then
          if (ielem(n, i)%ineighg(k) .gt. ielem(n, i)%ihexgl) then
            ielem(n, i)%nodes_faces(k, 1:ixf4) = ielem(n, ielem(n, i)%ineigh(k))%nodes_faces(ielem(n, i)%ineighn(k), 1:ixf4)
            ielem(n, i)%reorient(k) = 1
          end if
        end if
      end do
    end if
    surfmin = 1.0e16
    surfmax = 1.0e-16
    do l = 1, ielem(n, i)%ifca
      select case (ielem(n, i)%types_faces(l))
      case (5)
        nnd = 4
      case (6)
        nnd = 3
      end select
      if (ielem(n, i)%interior .eq. 1) then
        if ((ielem(n, i)%ineighg(l) .gt. 0) .and. (ielem(n, i)%ibounds(l) .gt. 0)) then
        end if
      else
        call compute_centre3df(n, i, l, nnd, cords)
        vext(2, 1:dims) = cords(1:dims)
      end if
      surfmin = min(surfmin, ielem(n, i)%surf(l))
      surfmax = max(surfmax, ielem(n, i)%surf(l))
      dist1 = distance3(n, vext)
      if (dist1 .lt. dxx2) then
        dxx2 = dist1
      end if
      if (dist1 .gt. dxx1) then
        dxx1 = dist1
      end if
    end do
    ielem(n, i)%condition = 1.0
    ielem(n, i)%condition = surfmax/surfmin
    if (code_profile .eq. 9) then
      if ((ielem(n, i)%condition .gt. 30) .or. (ielem(n, i)%ishape .eq. 4)) then
        ielem(n, i)%full = 0
      end if
    end if
    if (turbulence .gt. 0) then
      if (ielem(n, i)%ishape .eq. 2) then
        ielem(n, i)%condition = surfmax/surfmin
        if (ielem(n, i)%condition .gt. 30) then
          ielem(n, i)%hybrid = 1
        end if
      end if
      if (ielem(n, i)%ishape .eq. 3) then
        ielem(n, i)%condition = surfmax/surfmin
        if (ielem(n, i)%condition .gt. 10) then
          ielem(n, i)%hybrid = 1
        end if
      end if
    end if
    tempg1 = ielem(n, i)%condition!max((dxx1/dxx2),(dxx2/dxx1))
    ielem(n, i)%erx = tempg1
    if (code_profile .eq. 88) then
      if ((ielem(n, i)%ishape .eq. 3)) then
        ielem(n, i)%full = 0
      end if
    end if
    if (tempg1 .gt. gridar1) then
      ielem(n, i)%ggs = 1
      if ((iadapt .eq. 1) .or. (code_profile .eq. 88) .or. (code_profile .eq. 98)) then
        ielem(n, i)%full = 0
      end if
    end if
    if (fastest .eq. 0) then
      dxx1 = -tolbig; dxx2 = tolbig
      jj = 1
      do l = 1, ielem(n, i)%inumneighbours
        if (ilocal_recon3(i)%volume(jj, l) .lt. dxx2) then
          dxx2 = ilocal_recon3(i)%volume(jj, l)
        end if
        if (ilocal_recon3(i)%volume(jj, l) .gt. dxx1) then
          dxx1 = ilocal_recon3(i)%volume(jj, l)
        end if
      end do
      tempg1 = max((dxx1/dxx2), (dxx2/dxx1))
      !ielem(n,i)%walldist=tempg1
    end if
    idc = 0
    idc2 = 0
    if (ielem(n, i)%interior .eq. 1) then
      do j = 1, ielem(n, i)%ifca
        if (ielem(n, i)%ibounds(j) .gt. 0) then
          if (ibound(n, ielem(n, i)%ibounds(j))%icode .eq. 4) then
            idc = idc + 1
          else
            idc2 = idc2 + 1
          end if
        end if
      end do
    end if
    if ((idc .gt. 1)) then     !until ge is fully adaptive
      ielem(n, i)%ggs = 1
    end if
  end subroutine checkgrads
  subroutine check3(n, iconsi)
!> @brief
!> this subroutine assigns the ordering for the faces
    implicit none
    integer, intent(in)::n
    integer, intent(inout)::iconsi
    real::dxx1, dxx2, tempg1, dist1, dist2, oo2
    integer::i, j, k, l, jj, icount3, nnd, ixf4
    real, dimension(1:dimensiona)::cords
    real, dimension(1:8, 1:dimensiona)::vext
    do iconsi = 1, xmpielrank(n)
      i = iconsi
      ielem(n, i)%ggs = greengo
      dxx1 = -tolbig; dxx2 = tolbig
      call compute_centre3d(i, cords)
      vext(1, 1:dims) = cords(1:dims)
      i = iconsi
      if (ielem(n, i)%interior .eq. 1) then
      do k = 1, ielem(n, i)%ifca
        if (ielem(n, i)%types_faces(k) .eq. 5) then
          ixf4 = 4
        else
          ixf4 = 3
        end if
        if (ielem(n, i)%ineighg(k) .gt. 0) then
          if (ielem(n, i)%ineighb(k) .ne. n) then
            if (ielem(n, i)%ineighg(k) .gt. ielem(n, i)%ihexgl) then
           ielem(n, i)%nodes_faces(k, 1:ixf4) = iexboundhirr(ielem(n, i)%ineighn(k))%vertpp(ielem(n, i)%q_face(k)%q_mapl(1), 1:ixf4)
              ielem(n, i)%reorient(k) = 1
            end if
          else
            if (ielem(n, i)%ineighg(k) .gt. ielem(n, i)%ihexgl) then
              ielem(n, i)%nodes_faces(k, 1:ixf4) = ielem(n, ielem(n, i)%ineigh(k))%nodes_faces(ielem(n, i)%ineighn(k), 1:ixf4)
              ielem(n, i)%reorient(k) = 1
            end if
          end if
        end if
      end do
      else
      do k = 1, ielem(n, i)%ifca
        if (ielem(n, i)%types_faces(k) .eq. 5) then
          ixf4 = 4
        else
          ixf4 = 3
        end if
        if (ielem(n, i)%ineighg(k) .gt. 0) then
          if (ielem(n, i)%ineighg(k) .gt. ielem(n, i)%ihexgl) then
            ielem(n, i)%nodes_faces(k, 1:ixf4) = ielem(n, ielem(n, i)%ineigh(k))%nodes_faces(ielem(n, i)%ineighn(k), 1:ixf4)
            ielem(n, i)%reorient(k) = 1
          end if
        end if
      end do
      end if
    end do
  end subroutine check3
  subroutine checkgrads2d(n, iconsi)
!> @brief
!> this subroutine assigns the correct viscous gradient approximation flag for each cell based on some additional geometrical characteristics
    implicit none
    integer, intent(in)::n, iconsi
    real::dxx1, dxx2, tempg1, dist1, dist2, oo2, surfmin, surfmax
    integer::i, j, k, l, jj, icount3, nnd, ixf4
    real, dimension(1:dimensiona)::cords
    real, dimension(1:8, 1:dimensiona)::vext, nodes_list
    i = iconsi
    tempg1 = 0.0; 
    ielem(n, i)%ggs = greengo
    dxx1 = -tolbig; dxx2 = tolbig
    call compute_centre2d(i, cords)
    vext(1, 1:dims) = cords(1:dims)
    if (ielem(n, i)%interior .eq. 1) then
    do k = 1, ielem(n, i)%ifca
      ixf4 = 2
      if (ielem(n, i)%ineighg(k) .gt. 0) then
        if (ielem(n, i)%ineighb(k) .ne. n) then
          if (ielem(n, i)%ineighg(k) .gt. ielem(n, i)%ihexgl) then
           ielem(n, i)%nodes_faces(k, 1:ixf4) = iexboundhirr(ielem(n, i)%ineighn(k))%vertpp(ielem(n, i)%q_face(k)%q_mapl(1), 1:ixf4)
            ielem(n, i)%reorient(k) = 1
          end if
        else
          if (ielem(n, i)%ineighg(k) .gt. ielem(n, i)%ihexgl) then
            ielem(n, i)%nodes_faces(k, 1:ixf4) = ielem(n, ielem(n, i)%ineigh(k))%nodes_faces(ielem(n, i)%ineighn(k), 1:ixf4)
            ielem(n, i)%reorient(k) = 1
          end if
        end if
      end if
    end do
    else
    do k = 1, ielem(n, i)%ifca
      ixf4 = 2
      if (ielem(n, i)%ineighg(k) .gt. 0) then
        if (ielem(n, i)%ineighg(k) .gt. ielem(n, i)%ihexgl) then
          ielem(n, i)%nodes_faces(k, 1:ixf4) = ielem(n, ielem(n, i)%ineigh(k))%nodes_faces(ielem(n, i)%ineighn(k), 1:ixf4)
          ielem(n, i)%reorient(k) = 1
        end if
      end if
    end do
    end if
    surfmin = 1.0e16
    surfmax = 1.0e-16
    do l = 1, ielem(n, i)%ifca
      nnd = 2
      surfmin = min(surfmin, ielem(n, i)%surf(l))
      surfmax = max(surfmax, ielem(n, i)%surf(l))
      if (ielem(n, i)%interior .eq. 1) then
        if ((ielem(n, i)%ineighg(l) .gt. 0) .and. (ielem(n, i)%ibounds(l) .gt. 0)) then
          if (ielem(n, i)%ineighb(l) .ne. n) then
            do k = 1, nnd
              nodes_list(k, 1:dims) = inoder(ielem(n, iconsi)%nodes_faces(l, k))%cord(1:dims)
            end do
            do k = 1, nnd
            if (abs(nodes_list(k, 1) - vext(1, 1)) .gt. xper*oo2) then
              nodes_list(k, 1) = nodes_list(k, 1) + (xper*sign(1.0d0, vext(1, 1) - xper*oo2))
            end if
            if (abs(nodes_list(k, 2) - vext(1, 2)) .gt. yper*oo2) then
              nodes_list(k, 2) = nodes_list(k, 2) + (yper*sign(1.0d0, vext(1, 2) - yper*oo2))
            end if
            end do
            cords = cordinates2(n, nodes_list, nnd)
            vext(2, 1:dims) = cords(1:dims)
          else
            call compute_centre2df(n, i, l, nnd, cords)
            vext(2, 1:dims) = cords(1:dims)
          end if
        else
          call compute_centre2df(n, i, l, nnd, cords)
          vext(2, 1:dims) = cords(1:dims)
        end if
      else
        call compute_centre2df(n, i, l, nnd, cords)
        vext(2, 1:dims) = cords(1:dims)
      end if
      dist1 = distance2(n, vext)
      if (dist1 .lt. dxx2) then
        dxx2 = dist1
      end if
      if (dist1 .gt. dxx1) then
        dxx1 = dist1
      end if
    end do
    ielem(n, i)%condition = 1.0
    ielem(n, i)%condition = surfmax/surfmin
    tempg1 = ielem(n, i)%condition
    if (tempg1 .gt. gridar1) then
      ielem(n, i)%ggs = 1
      if ((iadapt .eq. 1) .or. (code_profile .eq. 88) .or. (code_profile .eq. 98)) then
        ielem(n, i)%full = 0
      end if
    end if
    if (fastest .eq. 0) then
      dxx1 = -tolbig; dxx2 = tolbig
      jj = 1
      do l = 1, ielem(n, i)%inumneighbours
        if (ilocal_recon3(i)%volume(jj, l) .lt. dxx2) then
          dxx2 = ilocal_recon3(i)%volume(jj, l)
        end if
        if (ilocal_recon3(i)%volume(jj, l) .gt. dxx1) then
          dxx1 = ilocal_recon3(i)%volume(jj, l)
        end if
      end do
if
    end if
  end subroutine checkgrads2d
end module local
