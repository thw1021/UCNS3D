module boundary
  ! @brief
  ! this module includes the subroutines for reading and establishing the boundary conditions
  ! for all the cells that need to be bounded including the periodic ones
  use declaration
  use library
  use transform
  implicit none

contains
  subroutine read_bound(n, imaxb, ibound, xmpielrank)
    ! @brief
    ! this subroutine reads the boundary conditions for all the bounded elements
    implicit none
    integer, intent(in)::n, imaxb
    type(bound_number), allocatable, dimension(:, :), intent(inout)::ibound
    character(len=12)::bndfile, ibx_code
    integer, allocatable, dimension(:), intent(in)::xmpielrank
    integer::i, j, ji, k, lm, iex, kmaxn, kk, kmaxe, kkk, jjj, jfx, jb, kxk, itr1, jj, jj1
    integer::ioy, ibid, ib1, ib2, ib3, ib4, ibx1, itl, ibgw, ibgw2, ibleed, ibdum
    integer, dimension(4)::ib_n
    integer::nb1, nb2, nb3, nb4

    kmaxe = xmpielrank(n)
    kkk = 0
    totiw = 0
    ibgw = 0
    ibgw2 = 0

    if (dimensiona .eq. 3) then
      bndfile = 'grid.bnd'
      if (binio .eq. 0)
        open (10, file=bndfile, form='formatted', status='old', action='read', iostat=ioy)
      if (binio .eq. 1)
        open (10, file=bndfile, form='unformatted', status='old', action='read', iostat=ioy)
      itl = 0
      if (binio .eq. 0) then
        do ji = 1, imaxb
          read (10, *) ibid, ib1, ib2, ib3, ib4
          if ((inoder(ib1)%itor .gt. 0) .and. (inoder(ib2)%itor .gt. 0) .and. (inoder(ib3)%itor .gt. 0) .and. (inoder(ib4)%itor .gt. 0)) then
            itl = itl + 1
          end if
        end do
      else
        do ji = 1, imaxb
          read (10) ibid, ib1, ib2, ib3, ib4
          if ((inoder(ib1)%itor .gt. 0) .and. (inoder(ib2)%itor .gt. 0) .and. (inoder(ib3)%itor .gt. 0) .and. (inoder(ib4)%itor .gt. 0)) then
            itl = itl + 1
          end if
        end do
      end if
      if (itl .gt. 0) then
        allocate(ibound(n:n, itl))
        ibound(n:n, :)%inum = 0
      end if
      close (10)

      itl = 0
      if (binio .eq. 0)
        open (10, file=bndfile, form='formatted', status='old', action='read', iostat=ioy)
      if (binio .eq. 1)
        open (10, file=bndfile, form='unformatted', status='old', action='read', iostat=ioy)

      if (binio .eq. 0) then
        do ji = 1, imaxb
          read (10, *) ibid, ib_n(1), ib_n(2), ib_n(3), ib_n(4), ibx1
          if (ibx1 .eq. 4) then
            ibgw = ibgw + 1
          end if
          if ((inoder(ib_n(1))%itor .gt. 0) .and. (inoder(ib_n(2))%itor .gt. 0) .and. (inoder(ib_n(3))%itor .gt. 0) .and. (inoder(ib_n(4))%itor .gt. 0)) then
            itl = itl + 1
            if (ibx1 .eq. 4) then
              ibound(n, itl)%inum = ibgw
              totiw = totiw + 1
            end if
            ibound(n, itl)%icode = ibx1
            ibound(n, itl)%ibid = ibid
            if (ib_n(3) .eq. ib_n(4)) then
              ibound(n, itl)%ishape = 6
              allocate(ibound(n, itl)%ibl(1:3))
              ibound(n, itl)%ibl(1:3) = ib_n(1:3)
            else
              ibound(n, itl)%ishape = 5
              allocate(ibound(n, itl)%ibl(1:4))
              ibound(n, itl)%ibl(1:4) = ib_n(1:4)
            end if
          end if
        end do
      else
        do ji = 1, imaxb
          read (10) ibid, ib_n(1), ib_n(2), ib_n(3), ib_n(4), ibx1
          if (ibx1 .eq. 4) then
            ibgw = ibgw + 1
          end if
          if ((inoder(ib_n(1))%itor .gt. 0) .and. (inoder(ib_n(2))%itor .gt. 0) .and. (inoder(ib_n(3))%itor .gt. 0) .and. (inoder(ib_n(4))%itor .gt. 0)) then
            itl = itl + 1
            if (ibx1 .eq. 4) then
              ibound(n, itl)%inum = ibgw
              totiw = totiw + 1
            end if
            ibound(n, itl)%icode = ibx1
            ibound(n, itl)%ibid = ibid
            if (ib_n(3) .eq. ib_n(4)) then
              ibound(n, itl)%ishape = 6
              allocate(ibound(n, itl)%ibl(1:3))
              ibound(n, itl)%ibl(1:3) = ib_n(1:3)
            else
              ibound(n, itl)%ishape = 5
              allocate(ibound(n, itl)%ibl(1:4))
              ibound(n, itl)%ibl(1:4) = ib_n(1:4)
            end if
          end if
        end do
      end if
      close (10)
      n_boundaries = itl
      itl = 0
      do ji = 1, n_boundaries
        if ((ibound(n, ji)%icode .eq. 5) .or. (ibound(n, ji)%icode .eq. 50)) then
          allocate(ibound(n, ji)%localn(2))
          allocate(ibound(n, ji)%cpun(2))
          itl = itl + 1
          ibound(n, ji)%localn = 0
          ibound(n, ji)%cpun = 0
        end if
      end do
      do i = 1, kmaxe
        if (ielem(n, i)%interior .eq. 1) then
          allocate(ielem(n, i)%ibounds(ielem(n, i)%ifca))
          ielem(n, i)%ibounds = 0
          ielem(n, i)%nofbc = 0
        end if
      end do

      itl = 0
      jj1 = 0
      do ji = 1, n_boundaries
        do jfx = 1, inoder2(ibound(n, ji)%ibl(1))%numberofneib
          j = inoder2(ibound(n, ji)%ibl(1))%neibids(jfx)
          if (ielem(n, j)%interior .eq. 1) then
            do jj = 1, ielem(n, j)%ifca
              if (ielem(n, j)%ineighg(jj) .eq. 0) then
                if (ielem(n, j)%types_faces(jj) .eq. ibound(n, ji)%ishape) then
                  if (ibound(n, ji)%ishape .eq. 6) then
                    nb1 = ielem(n, j)%nodes_faces(jj, 1)
                    nb2 = ielem(n, j)%nodes_faces(jj, 2)
                    nb3 = ielem(n, j)%nodes_faces(jj, 3)
                    ib1 = ibound(n, ji)%ibl(1)
                    ib2 = ibound(n, ji)%ibl(2)
                    ib3 = ibound(n, ji)%ibl(3)
                    if (((nb1 .eq. ib1) .or. (nb1 .eq. ib2) .or. (nb1 .eq. ib3)) .and. &
                          ((nb2 .eq. ib1) .or. (nb2 .eq. ib2) .or. (nb2 .eq. ib3)) .and. &
                          ((nb3 .eq. ib1) .or. (nb3 .eq. ib2) .or. (nb3 .eq. ib3))) then
                      ielem(n, j)%ibounds(jj) = ji
                      ibound(n, ji)%which = j
                      ibound(n, ji)%face = jj
                      ielem(n, j)%nofbc = ielem(n, j)%nofbc + 1
                      if ((ibound(n, ji)%icode .eq. 5) .or. (ibound(n, ji)%icode .eq. 50)) then
                        ibound(n, ji)%localn(1) = j; ibound(n, ji)%cpun(1) = n
                        jj1 = jj1 + 1
                        go to 51
                      end if
                    end if
                  end if
                  if (ibound(n, ji)%ishape .eq. 5) then
                    nb1 = ielem(n, j)%nodes_faces(jj, 1)
                    nb2 = ielem(n, j)%nodes_faces(jj, 2)
                    nb3 = ielem(n, j)%nodes_faces(jj, 3)
                    nb4 = ielem(n, j)%nodes_faces(jj, 4)
                    ib1 = ibound(n, ji)%ibl(1)
                    ib2 = ibound(n, ji)%ibl(2)
                    ib3 = ibound(n, ji)%ibl(3)
                    ib4 = ibound(n, ji)%ibl(4)
                    if (((nb1 .eq. ib1) .or. (nb1 .eq. ib2) .or. (nb1 .eq. ib3) .or. (nb1 .eq. ib4)) .and. &
                          ((nb2 .eq. ib1) .or. (nb2 .eq. ib2) .or. (nb2 .eq. ib3) .or. (nb2 .eq. ib4)) .and. &
                          ((nb3 .eq. ib1) .or. (nb3 .eq. ib2) .or. (nb3 .eq. ib3) .or. (nb3 .eq. ib4)) .and. &
                          ((nb4 .eq. ib1) .or. (nb4 .eq. ib2) .or. (nb4 .eq. ib3) .or. (nb4 .eq. ib4))) then
                      ielem(n, j)%ibounds(jj) = ji
                      ibound(n, ji)%which = j
                      ibound(n, ji)%face = jj
                      ielem(n, j)%nofbc = ielem(n, j)%nofbc + 1
                      if ((ibound(n, ji)%icode .eq. 5) .or. (ibound(n, ji)%icode .eq. 50)) then
                        ibound(n, ji)%localn(1) = j; ibound(n, ji)%cpun(1) = n
                        jj1 = jj1 + 1
                        go to 51
                      end if
                    end if
                  end if
                end if
              end if
            end do
          end if
        end do
  51    continue
      end do
    end if

    if (dimensiona .eq. 2) then
      bndfile = 'grid.bnd'
      if (binio .eq. 0)
        open (10, file=bndfile, form='formatted', status='old', action='read', iostat=ioy)
      if (binio .eq. 1)
        open (10, file=bndfile, form='unformatted', status='old', action='read', iostat=ioy)
      itl = 0
      if (binio .eq. 0) then
        do ji = 1, imaxb
          read (10, *) ibid, ib1, ib2
          if ((inoder(ib1)%itor .gt. 0) .and. (inoder(ib2)%itor .gt. 0)) then
            itl = itl + 1
          end if
        end do
      else
        do ji = 1, imaxb
          read (10) ibid, ib1, ib2
          if ((inoder(ib1)%itor .gt. 0) .and. (inoder(ib2)%itor .gt. 0)) then
            itl = itl + 1
          end if
        end do
      end if

      if (itl .gt. 0) then
        allocate(ibound(n:n, itl))
        ibound(n:n, :)%inum = 0
      end if

      close (10)

      itl = 0
      if (binio .eq. 0)
        open (10, file=bndfile, form='formatted', status='old', action='read', iostat=ioy)
      if (binio .eq. 1)
        open (10, file=bndfile, form='unformatted', status='old', action='read', iostat=ioy)

      if (binio .eq. 0) then
        do ji = 1, imaxb
          read (10, *) ibid, ib_n(1), ib_n(2), ib_n(3), ib_n(4), ibx1
          if (bleed .eq. 1) then
            ibdum = ibx1
            if (ibx1 .eq. 4) then
              if ((inoder(ib_n(1))%itor .gt. 0) .and. (inoder(ib_n(2))%itor .gt. 0)) then
                !now check if this is in a bleed zone
                do ibleed = 1, bleed_number
                  if (((inoder(ib_n(1))%cord(1).ge.bleed_start(ibleed,1)).and.(inoder(ib_n(1))%cord(1).le.bleed_end(ibleed,1))).and.((inoder(ib_n(2))%cord(1).ge.bleed_start(ibleed,1)).and.(inoder(ib_n(2))%cord(1).le.bleed_end(ibleed,1))).and.((inoder(ib_n(1))%cord(2).ge.bleed_start(ibleed,2)).and.(inoder(ib_n(1))%cord(2).le.bleed_end(ibleed,2))).and.((inoder(ib_n(2))%cord(2).ge.bleed_start(ibleed,2)).and.(inoder(ib_n(2))%cord(2).le.bleed_end(ibleed,2))))then
                    ibdum = 99
                  end if
                end do
              end if
            end if
            ibx1 = ibdum
          end if

          if ((ibx1 .eq. 4) .or. (ibx1 .eq. 99)) then
            ibgw = ibgw + 1
          end if
          if ((inoder(ib_n(1))%itor .gt. 0) .and. (inoder(ib_n(2))%itor .gt. 0)) then
            itl = itl + 1
            if ((ibx1 .eq. 4) .or. (ibx1 .eq. 99)) then
              ibound(n, itl)%inum = ibgw
              totiw = totiw + 1
            end if
            ibound(n, itl)%icode = ibx1
            ibound(n, itl)%ibid = ibid
            ibound(n, itl)%ishape = 7
            allocate(ibound(n, itl)%ibl(1:2))
            ibound(n, itl)%ibl(1:2) = ib_n(1:2)
          end if
        end do
      else
        do ji = 1, imaxb
          read (10) ibid, ib_n(1), ib_n(2), ib_n(3), ib_n(4), ibx1
          if (bleed .eq. 1) then
            ibdum = ibx1
            if (ibx1 .eq. 4) then
              if ((inoder(ib_n(1))%itor .gt. 0) .and. (inoder(ib_n(2))%itor .gt. 0)) then
                ! now check if this is in a bleed zone
                do ibleed = 1, bleed_number
                  if (((inoder(ib_n(1))%cord(1).ge.bleed_start(ibleed,1)).and.(inoder(ib_n(1))%cord(1).le.bleed_end(ibleed,1))).and.((inoder(ib_n(2))%cord(1).ge.bleed_start(ibleed,1)).and.(inoder(ib_n(2))%cord(1).le.bleed_end(ibleed,1))).and.((inoder(ib_n(1))%cord(2).ge.bleed_start(ibleed,2)).and.(inoder(ib_n(1))%cord(2).le.bleed_end(ibleed,2))).and.((inoder(ib_n(2))%cord(2).ge.bleed_start(ibleed,2)).and.(inoder(ib_n(2))%cord(2).le.bleed_end(ibleed,2)))) then
                    ibdum = 99
                  end if
                end do
              end if
            end if
            ibx1 = ibdum
          end if
          if ((ibx1 .eq. 4) .or. (ibx1 .eq. 99)) then
            ibgw = ibgw + 1
          end if
          if ((inoder(ib_n(1))%itor .gt. 0) .and. (inoder(ib_n(2))%itor .gt. 0)) then
            itl = itl + 1
            if ((ibx1 .eq. 4) .or. (ibx1 .eq. 99)) then
              ibound(n, itl)%inum = ibgw; totiw = totiw + 1
            end if
            ibound(n, itl)%icode = ibx1
            ibound(n, itl)%ibid = ibid
            ibound(n, itl)%ishape = 7
            allocate(ibound(n, itl)%ibl(1:2))
            ibound(n, itl)%ibl(1:2) = ib_n(1:2)
          end if
        end do
      end if
      close (10)
      n_boundaries = itl
      do ji = 1, n_boundaries
        if ((ibound(n, ji)%icode .eq. 5) .or. (ibound(n, ji)%icode .eq. 50)) then
          allocate(ibound(n, ji)%localn(2))
          allocate(ibound(n, ji)%cpun(2))
          ibound(n, ji)%localn = 0
          ibound(n, ji)%cpun = 0
        end if
      end do

      do i = 1, kmaxe
        if (ielem(n, i)%interior .eq. 1) then
          allocate(ielem(n, i)%ibounds(ielem(n, i)%ifca))
          if (bleed .eq. 1) then
            allocate(ielem(n, i)%bleedn(ielem(n, i)%ifca))
            ielem(n, i)%bleedn(:) = 0
          end if
          ielem(n, i)%ibounds = 0
          ielem(n, i)%nofbc = 0
        end if
      end do

      itl = 0
      do ji = 1, n_boundaries
        do jfx = 1, inoder2(ibound(n, ji)%ibl(1))%numberofneib
          j = inoder2(ibound(n, ji)%ibl(1))%neibids(jfx)
          if (ielem(n, j)%interior .eq. 1) then
            do jj = 1, ielem(n, j)%ifca
              if (ielem(n, j)%ineighg(jj) .eq. 0) then
                nb1 = ielem(n, j)%nodes_faces(jj, 1)
                nb2 = ielem(n, j)%nodes_faces(jj, 2)
                ib1 = ibound(n, ji)%ibl(1)
                ib2 = ibound(n, ji)%ibl(2)
                if (((nb1 .eq. ib1) .or. (nb1 .eq. ib2)) .and. &
                      ((nb2 .eq. ib1) .or. (nb2 .eq. ib2))) then
                  ielem(n, j)%ibounds(jj) = ji
                  ibound(n, ji)%which = j
                  ibound(n, ji)%face = jj
                  ielem(n, j)%nofbc = ielem(n, j)%nofbc + 1
                  if ((ibound(n, ji)%icode .eq. 5) .or. (ibound(n, ji)%icode .eq. 50)) then
                    ibound(n, ji)%localn(1) = j; ibound(n, ji)%cpun(1) = n
                    itl = itl + 1
                  end if
                end if
              end if
            end do
          end if
        end do
      end do
    end if
    totwalls = ibgw

    if (totiw .gt. 0) then
      allocate(ibound_t(totiw))
      allocate(ibound_t2(totiw))
      totiw = 0
      do i = 1, kmaxe
        if (ielem(n, i)%interior .eq. 1) then
          do j = 1, ielem(n, i)%ifca
            if (ielem(n, i)%ibounds(j) .gt. 0) then
              if ((ibound(n, ielem(n, i)%ibounds(j))%icode .eq. 4) .or. (ibound(n, ielem(n, i)%ibounds(j))%icode .eq. 99)) then
                totiw = totiw + 1
                ibound_t(totiw) = i
                ibound_t2(totiw) = j
              end if
            end if
          end do
        end if
      end do
    end if

    call mpi_barrier(mpi_comm_world, ierror)
    ! deallocate(ibound)
  end subroutine read_bound

  subroutine apply_boundary(n, xper, yper, zper, iperiodicity, xmpielrank)
    ! @brief
    ! this subroutine assigns the correct boundary condition code for all the bounded elements
    implicit none
    real, intent(in)::xper, yper, zper
    integer, intent(in)::iperiodicity, n
    real::small, tolerance, dist, temp_x
    integer::i, k, j, kk, ii, kmaxe, jj1, jj2, ji, l, ibleed
    integer, allocatable, dimension(:), intent(in)::xmpielrank
    integer::dum1, dum2, n_node
    real, dimension(1:8, 1:dimensiona)::vext, nodes_list

    kmaxe = xmpielrank(n)
    tolerance = tolsmall ! the tolerance can have a significant impact on the periodic boundary conditions matching rules
    jj1 = 0
    if (dimensiona .eq. 3) then
      jj2 = 0
      do i = 1, n_boundaries
        if ((ibound(n, i)%icode .eq. 5) .or. (ibound(n, i)%icode .eq. 50)) then
          jj2 = jj2 + 1
        end if
      end do
      do i = 1, kmaxe                        ! for all elements
        if (ielem(n, i)%interior .eq. 1) then                ! that have at least one unknwon neighbour
          if (ielem(n, i)%nofbc .gt. 0) then                ! that have at least established a boundary condition code
            do j = 1, ielem(n, i)%ifca                        ! loop all their faces
              if (ielem(n, i)%ibounds(j) .gt. 0) then
                if ((ibound(n, ielem(n, i)%ibounds(j))%icode .eq. 5) .or. (ibound(n, ielem(n, i)%ibounds(j))%icode .eq. 50)) then        !if any of them has a periodic boundary condition then
                  if (ibound(n, ielem(n, i)%ibounds(j))%ishape .eq. 5) then
                    n_node = 4
                  else
                    n_node = 3
                  end if
                  do kk = 1, n_node
                    nodes_list(kk, 1:3) = inoder(ibound(n, ielem(n, i)%ibounds(j))%ibl(kk))%cord(1:3)
                  end do
                  vext(1, 1:3) = cordinates3(n, nodes_list, n_node)
                  do ii = 1, n_boundaries                                ! loop all the boundaries
                    if ((ii .ne. ielem(n, i)%ibounds(j)) .and. ((ibound(n, ii)%icode .eq. 5) .or. (ibound(n, ii)%icode .eq. 50)) .and. (ibound(n, ielem(n, i)%ibounds(j))%ishape .eq. ibound(n, ii)%ishape)) then
                      if ((ibound(n, ii)%localn(1) .gt. 0)) then        ! excluding itself, and of same shape type
                        if (ielem(n, ibound(n, ii)%localn(1))%ihexgl .ne. ielem(n, i)%ihexgl) then
                          do kk = 1, n_node
                            nodes_list(kk, 1:3) = inoder(ibound(n, ii)%ibl(kk))%cord(1:3)
                          end do
                          vext(2, 1:3) = cordinates3(n, nodes_list, n_node)
                          dist = distance3(n, vext)
                          if (per_rot .eq. 0) then
                            if (((abs(vext(2, 1) - xper) .lt. tolsmall) .or. (abs((abs(vext(2, 1) - xper)) - xper) .lt. tolsmall)) .and. ((abs(vext(1, 1) - xper) .lt. tolsmall) .or. (abs((abs(vext(1, 1) - xper)) - xper) .lt. tolsmall))) then
                              if ((abs(vext(2, 2) - vext(1, 2)) .lt. tolsmall) .and. (abs(vext(2, 3) - vext(1, 3)) .lt. tolsmall)) then
                                ibound(n, ii)%localn(2) = i; ibound(n, ii)%cpun(2) = n
                                ielem(n, i)%ineighg(j) = ielem(n, ibound(n, ii)%localn(1))%ihexgl
                                jj1 = jj1 + 1
                                go to 101
                              end if
                            end if
                            if (((abs(vext(2, 2) - yper) .lt. tolsmall) .or. (abs((abs(vext(2, 2) - yper)) - yper) .lt. tolsmall)) .and. ((abs(vext(1, 2) - yper) .lt. tolsmall) .or. (abs((abs(vext(1, 2) - yper)) - yper) .lt. tolsmall))) then
                              if ((abs(vext(2, 1) - vext(1, 1)) .lt. tolsmall) .and. (abs(vext(2, 3) - vext(1, 3)) .lt. tolsmall)) then
                                ibound(n, ii)%localn(2) = i; ibound(n, ii)%cpun(2) = n
                                ielem(n, i)%ineighg(j) = ielem(n, ibound(n, ii)%localn(1))%ihexgl
                                jj1 = jj1 + 1
                                go to 101
                              end if
                            end if
                          if (((abs(vext(2, 3) - zper) .lt. tolsmall) .or. (abs((abs(vext(2, 3) - zper)) - zper) .lt. tolsmall)) .and. ((abs(vext(1, 3) - zper) .lt. tolsmall) .or. (abs((abs(vext(1, 3) - zper)) - zper) .lt. tolsmall))) then
                            if ((abs(vext(2, 2) - vext(1, 2)) .lt. tolsmall) .and. (abs(vext(2, 1) - vext(1, 1)) .lt. tolsmall)) then
                              ibound(n, ii)%localn(2) = i; ibound(n, ii)%cpun(2) = n
                              ielem(n, i)%ineighg(j) = ielem(n, ibound(n, ii)%localn(1))%ihexgl
                              jj1 = jj1 + 1
                              go to 101
                            end if
                          end if
                        else
                          vext(2, :) = rotate_per(vext(2, :), ibound(n, ii)%icode, angle_per)
                          if ((abs(vext(1, 1) - vext(2, 1)) .lt. tol_per) .and. &
                              (abs(vext(1, 2) - vext(2, 2)) .lt. tol_per) .and. &
                              (abs(vext(1, 3) - vext(2, 3)) .lt. tol_per)) then
                            ibound(n, ii)%localn(2) = i
                            ibound(n, ii)%cpun(2) = n
                            ielem(n, i)%ineighg(j) = ielem(n, ibound(n, ii)%localn(1))%ihexgl
                            jj1 = jj1 + 1
                            go to 101
                          end if
                        end if
                      end if
                    end if
                  end if
                end do
101           continue
              end if
            end if
          end do
        end if
      end if
    end do
    else
      jj2 = 0
      do i = 1, n_boundaries
        if (ibound(n, i)%icode .eq. 5) then
          jj2 = jj2 + 1
        end if
      end do
      jj1 = 0
      do i = 1, kmaxe                        !> all elements
        if (ielem(n, i)%interior .eq. 1) then                ! that have at least one unknwon neighbour
          if (ielem(n, i)%nofbc .gt. 0) then                ! that have at least established a boundary condition code
            do j = 1, ielem(n, i)%ifca                        ! loop all their boundary faces
              if (ielem(n, i)%ibounds(j) .gt. 0) then
                if (bleed .eq. 1) then
                  if (ibound(n, ielem(n, i)%ibounds(j))%icode .eq. 99) then
                    !assign the bleed zone
                    n_node = 2
                    do kk = 1, n_node
                      nodes_list(kk, 1:2) = inoder(ibound(n, ielem(n, i)%ibounds(j))%ibl(kk))%cord(1:2)
                    end do
                    vext(1, :) = cordinates2(n, nodes_list, n_node)
                    do ibleed = 1, bleed_number
                      if (((vext(1,1).ge.bleed_start(ibleed,1)).and.(vext(1,1).le.bleed_end(ibleed,1))).and.((vext(1,2).ge.bleed_start(ibleed,2)).and.(vext(1,2).le.bleed_end(ibleed,2))))then
                        ielem(n, i)%bleedn(j) = ibleed        !assign the bleed number
                      end if
                    end do
                  end if
                end if
                if (ibound(n, ielem(n, i)%ibounds(j))%icode .eq. 5) then        ! if any of them has a periodic boundary condition then
                  n_node = 2
                  do kk = 1, n_node
                    nodes_list(kk, 1:2) = inoder(ibound(n, ielem(n, i)%ibounds(j))%ibl(kk))%cord(1:2)
                  end do
                  vext(1, :) = cordinates2(n, nodes_list, n_node)
                  do ii = 1, n_boundaries                                ! loop all the boundaries
                    if (((ii .ne. ielem(n, i)%ibounds(j)) .and. (ibound(n, ii)%icode .eq. 5))) then
                      if (ibound(n, ii)%localn(1) .gt. 0) then        ! excluding itself, and of same shape type
                        if (ielem(n, ibound(n, ii)%localn(1))%ihexgl .ne. ielem(n, i)%ihexgl) then
                          do kk = 1, n_node
                            nodes_list(kk, 1:2) = inoder(ibound(n, ii)%ibl(kk))%cord(1:2)
                          end do
                          vext(2, :) = cordinates2(n, nodes_list, n_node)
                          dist = distance2(n, vext)
                          if (((abs(vext(2, 1) - xper) .lt. tolsmall) .or. (abs((abs(vext(2, 1) - xper)) - xper) .lt. tolsmall)) .and. ((abs(vext(1, 1) - xper) .lt. tolsmall) .or. (abs((abs(vext(1, 1) - xper)) - xper) .lt. tolsmall))) then
                            if ((abs(vext(2, 2) - vext(1, 2)) .lt. tolsmall)) then
                              ibound(n, ii)%localn(2) = i; ibound(n, ii)%cpun(2) = n
                              ielem(n, i)%ineighg(j) = ielem(n, ibound(n, ii)%localn(1))%ihexgl
                              jj1 = jj1 + 1
                              go to 201
                            end if
                          end if
                          if (((abs(vext(2, 2) - yper) .lt. tolsmall) .or. (abs((abs(vext(2, 2) - yper)) - yper) .lt. tolsmall)) .and. ((abs(vext(1, 2) - yper) .lt. tolsmall) .or. (abs((abs(vext(1, 2) - yper)) - yper) .lt. tolsmall))) then
                            if ((abs(vext(2, 1) - vext(1, 1)) .lt. tolsmall)) then
                              ibound(n, ii)%localn(2) = i; ibound(n, ii)%cpun(2) = n
                              ielem(n, i)%ineighg(j) = ielem(n, ibound(n, ii)%localn(1))%ihexgl
                              jj1 = jj1 + 1
                              go to 201
                            end if
                          end if
                        end if
                      end if
                    end if
                  end do
201               continue
                end if
              end if
            end do
          end if
        end if
      end do
    end if
  end subroutine apply_boundary
end module boundary
