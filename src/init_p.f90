module initialisation
  use declaration
  use library
  use transform
  use profile
  implicit none

  contains

  subroutine initialise(n)
    implicit none
    integer, intent(in)::n
    real, allocatable, dimension(:)::rg, arg
    character(len=20)::proc, restfile, proc3
    integer::prev_turbequation, initial, iii, i, k, jx, qqp, inc, kmaxe, jkn, ki, iterr, jx2, ind1, kx, icompwrt, iconsidered
    integer::eltype, elem_dec, count_1
    real::voltemp
    real, dimension(1:nof_variables + turbulenceequations + passivescalar)::veccos
    real, dimension(1:dimensiona)::pox, poy, poz
    real, dimension(1:dimensiona)::cords
    real, dimension(1:8, 1:dimensiona)::vext, nodes_list
    real, dimension(1:6, 1:4, 1:dimensiona)::elem_listd
    real, dimension(1:dimensiona, 1:numberofpoints)::qpoints
    real, dimension(1:numberofpoints)::wequa3d

    if (rungekutta .eq. 4) then
      ind1 = 7
    else
      ind1 = 5
    end if

    if (lamps .eq. 1) then
      iii = 1
    else
      iii = 0
    end if

    prev_turbequation = 0
    if (prev_turbmodel .eq. 1) then
      prev_turbequation = 1
    end if
    if (prev_turbmodel .eq. 2) then
      prev_turbequation = 2
    end if
    if (lamps .eq. 1) then
      allocate(rg(nof_variables + prev_turbequation + 1))
    else
      allocate(rg(nof_variables + prev_turbequation))
    end if

    veccos(1:nof_variables) = 0.0
    kmaxe = xmpielrank(n)

    if (restart .eq. 0) then

      if (ischeme .lt. 2) then
        if (itestcase .eq. 0) then
          do initial = 1, kmaxe
            u_c(initial)%val(1, 1) = 1.0d0
            u_e(initial)%val(1, 1) = u_c(initial)%val(1, 1)
          end do
        end if
        if (itestcase .eq. 1) then
          do initial = 1, kmaxe
            pox(1) = ielem(n, initial)%xxc
            poy(1) = ielem(n, initial)%yyc
            if (dimensiona .eq. 3) then
              poz(1) = ielem(n, initial)%zzc
              u_c(initial)%val(1, 1) = linear_init3d(n, pox, poy, poz)
            else
              u_c(initial)%val(1, 1) = linear_init2d(n, pox, poy, poz)
            end if

            u_e(initial)%val(1, 1) = u_c(initial)%val(1, 1)
          end do
        end if
        if (itestcase .eq. 2) then
          do initial = 1, kmaxe
            pox(1) = ielem(n, initial)%xxc
            poy(1) = ielem(n, initial)%yyc
            if (dimensiona .eq. 3) then
              poz(1) = ielem(n, initial)%zzc
              u_c(initial)%val(1, 1) = linear_init3d(n, pox, poy, poz)
            else
              u_c(initial)%val(1, 1) = linear_init2d(n, pox, poy, poz)
            end if

            u_e(initial)%val(1, 1) = u_c(initial)%val(1, 1)
          end do
        end if
        if (itestcase .ge. 3) then
          do initial = 1, kmaxe
            veccos(:) = zero
            pox(1) = ielem(n, initial)%xxc
            poy(1) = ielem(n, initial)%yyc
            if (dimensiona .eq. 3) then
              poz(1) = ielem(n, initial)%zzc
              call initialise_euler3d(n, veccos, pox, poy, poz)
            else
              call initialise_euler2d(n, veccos, pox, poy, poz)
            end if

            if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
              u_c(initial)%val(1, 1:nof_variables) = veccos(1:nof_variables)
              u_ct(initial)%val(1,1:0+turbulenceequations+passivescalar)=veccos(nof_variables+1:nof_variables+turbulenceequations+passivescalar)
            else
              u_c(initial)%val(1, :) = veccos(:)
              if (itestcase .ge. 3) u_e(initial)%val(1, :) = u_c(initial)%val(1, :)
            end if
          end do
        end if
      else

        do i = 1, kmaxe
          iconsidered = i
          vext = zero
          nodes_list = zero
          eltype = ielem(n, i)%ishape
          elem_dec = ielem(n, i)%vdec
          elem_listd = zero
          voltemp = zero
          jx = ielem(n, i)%nonodes
          do k = 1, jx
            jx2 = ielem(n, i)%nodes(k)
            nodes_list(k, :) = inoder(jx2)%cord(:)
            vext(k, :) = nodes_list(k, :)
          end do
          if (dimensiona .eq. 3) then
            call decompose3(n, eltype, nodes_list, elem_listd)
          else
            call decompose2(n, eltype, nodes_list, elem_listd)
          end if
          select case (ielem(n, i)%ishape)
          case (1)
            count_1 = 0
            if (ielem(n, i)%mode .eq. 0) then
              call quadraturehexa(n, igqrules, vext, qpoints, wequa3d)
              voltemp = hexavolume(n, vext, qpoints, wequa3d)/ielem(n, i)%totvolume
              qqp = qp_hexa
              do inc = 1, qqp
                pox(1) = qpoints(1, inc)
                poy(1) = qpoints(2, inc)
                poz(1) = qpoints(3, inc)
                call assign_initial(i, inc, voltemp, wequa3d, pox, poy, poz, count_1)
              end do

            else
              count_1 = 0
              do k = 1, elem_dec
                vext(1:4, 1:3) = elem_listd(k, 1:4, 1:3)

                call quadraturetetra(n, igqrules, vext, qpoints, wequa3d)
                if (dg .eq. 1) then
                  voltemp = tetravolume(n, vext)
                else
                  voltemp = tetravolume(n, vext)/ielem(n, i)%totvolume
                end if
                qqp = qp_tetra
                do inc = 1, qqp
                  count_1 = count_1 + 1
                  pox(1) = qpoints(1, inc)
                  poy(1) = qpoints(2, inc)
                  poz(1) = qpoints(3, inc)

                  call assign_initial(i, inc, voltemp, wequa3d, pox, poy, poz, count_1)
                end do
              end do
            end if

          case (2)
            vext(1:4, 1:3) = elem_listd(1, 1:4, 1:3)
            call quadraturetetra(n, igqrules, vext, qpoints, wequa3d)
            if (dg .eq. 1) then
              voltemp = tetravolume(n, vext)
              count_1 = 0
            else
              voltemp = tetravolume(n, vext)/ielem(n, i)%totvolume
            end if

            qqp = qp_tetra
            do inc = 1, qqp
              count_1 = count_1 + 1
              pox(1) = qpoints(1, inc)
              poy(1) = qpoints(2, inc)
              poz(1) = qpoints(3, inc)
              call assign_initial(i, inc, voltemp, wequa3d, pox, poy, poz, count_1)
            end do
          case (3)

            count_1 = 0
            do k = 1, elem_dec
              vext(1:4, 1:3) = elem_listd(k, 1:4, 1:3)
              call quadraturetetra(n, igqrules, vext, qpoints, wequa3d)
              if (dg .eq. 1) then
                voltemp = tetravolume(n, vext)
              else
                voltemp = tetravolume(n, vext)/ielem(n, i)%totvolume
              end if
              qqp = qp_tetra
              do inc = 1, qqp
                count_1 = count_1 + 1
                pox(1) = qpoints(1, inc)
                poy(1) = qpoints(2, inc)
                poz(1) = qpoints(3, inc)
                call assign_initial(i, inc, voltemp, wequa3d, pox, poy, poz, count_1)
              end do
            end do

          case (4)
            if (ielem(n, i)%mode .eq. 0) then
              call quadratureprism(n, igqrules, vext, qpoints, wequa3d)
              voltemp = prismvolume(n, vext, qpoints, wequa3d)/ielem(n, i)%totvolume
              qqp = qp_prism
              do inc = 1, qqp
                pox(1) = qpoints(1, inc)
                poy(1) = qpoints(2, inc)
                poz(1) = qpoints(3, inc)
                count_1 = inc
                call assign_initial(i, inc, voltemp, wequa3d, pox, poy, poz, count_1)
              end do
            else
              count_1 = 0
              do k = 1, elem_dec
                vext(1:4, 1:3) = elem_listd(k, 1:4, 1:3)
                call quadraturetetra(n, igqrules, vext, qpoints, wequa3d)
                if (dg .eq. 1) then
                  voltemp = tetravolume(n, vext)
                else
                  voltemp = tetravolume(n, vext)/ielem(n, i)%totvolume
                end if
                qqp = qp_tetra
                do inc = 1, qqp
                  count_1 = count_1 + 1
                  pox(1) = qpoints(1, inc)
                  poy(1) = qpoints(2, inc)
                  poz(1) = qpoints(3, inc)
                  call assign_initial(i, inc, voltemp, wequa3d, pox, poy, poz, count_1)
                end do
              end do
            end if

          case (5)
            if (ielem(n, i)%mode .eq. 0) then
              call quadraturequad(n, igqrules, vext, qpoints, wequa3d)
              voltemp = 1.0d0
              qqp = qp_quad
              do inc = 1, qqp
                count_1 = qqp
                pox(1) = qpoints(1, inc) !pox,poy required for linear_init2d
                poy(1) = qpoints(2, inc)
                call assign_initial(i, inc, voltemp, wequa3d, pox, poy, poz, count_1)
              end do
            else
              count_1 = 0
              voltemp = 0.0d0
              wequa3d = 0.0d0
              qpoints = 0.0d0
              do k = 1, elem_dec
                vext(1:3, 1:2) = elem_listd(k, 1:3, 1:2)
                call quadraturetriangle(n, igqrules, vext, qpoints, wequa3d)
                if (dg .eq. 1) then
                  voltemp = trianglevolume(n, vext)
                else
                  voltemp = trianglevolume(n, vext)/ielem(n, i)%totvolume
                end if

                qqp = qp_triangle
                do inc = 1, qqp
                  count_1 = count_1 + 1
                  pox(1) = qpoints(1, inc) !pox,poy required for linear_init2d
                  poy(1) = qpoints(2, inc)
                  call assign_initial(i, inc, voltemp, wequa3d, pox, poy, poz, count_1)
                end do
              end do
            end if

          case (6)
            count_1 = 0
            call quadraturetriangle(n, igqrules, vext, qpoints, wequa3d)
            if (dg .eq. 1) then
              voltemp = trianglevolume(n, vext)
            else
              voltemp = trianglevolume(n, vext)/ielem(n, i)%totvolume
            end if
            qqp = qp_triangle
            do inc = 1, qqp
              count_1 = count_1 + 1
              pox(1) = qpoints(1, inc) !pox,poy required for linear_init2d
              poy(1) = qpoints(2, inc)
              call assign_initial(i, inc, voltemp, wequa3d, pox, poy, poz, count_1)
            end do
          end select
        end do
      end if

      res_time = zero
    end if
    initialres = 0.0d0

    deallocate(rg)
  end subroutine initialise

  subroutine assign_initial(iconsidered, inc, voltemp, wequa3d, pox, poy, poz, count_1)
    implicit none
    integer, intent(in)::iconsidered, inc, count_1
    real, dimension(1:dimensiona), intent(in)::pox, poy, poz
    real, intent(in)::voltemp
    real, dimension(1:nof_variables + turbulenceequations + passivescalar)::veccos
    real, dimension(1:idegfree + 1)::basis_vector
    real, dimension(1, 1:idegfree + 1)::tempsol
    real, dimension(1:numberofpoints), intent(in)::wequa3d
    integer::i, kx, icompwrt

    i = iconsidered

    if (itestcase .le. 2) then
      if (dg .eq. 1) then
        basis_vector(1) = 1.0d0
        icompwrt = -2

        if (dimensiona .eq. 2) then
        basis_vector(2:idegfree + 1) = basis_rec2d(n, qp_array(i)%x(count_1), qp_array(i)%y(count_1), iorder, i, idegfree, icompwrt)
        else
            basis_vector(2:idegfree+1) = basis_rec(n,qp_array(i)%x(count_1),qp_array(i)%y(count_1),qp_array(i)%z(count_1),iorder,i,idegfree,icompwrt)
        end if
        if (dimensiona .eq. 2) then
          tempsol(1, :) = linear_init2d(n, pox, poy, poz)*wequa3d(inc)*(voltemp)*basis_vector(1:idegfree + 1)
        else
          tempsol(1, :) = linear_init3d(n, pox, poy, poz)*wequa3d(inc)*(voltemp)*basis_vector(1:idegfree + 1)
        end if
        u_c(i)%valdg(1, 1, :) = u_c(i)%valdg(1, 1, :) + matmul(m_1(i)%val(:, :), tempsol(1, :))
      else
        if (dimensiona .eq. 2) then
          u_c(i)%val(1, 1) = u_c(i)%val(1, 1) + linear_init2d(n, pox, poy, poz)*wequa3d(inc)*(voltemp)
        else
          u_c(i)%val(1, 1) = u_c(i)%val(1, 1) + linear_init3d(n, pox, poy, poz)*wequa3d(inc)*(voltemp)
        end if
      end if
      u_e(i)%val(1, 1) = u_c(i)%val(1, 1)
      if (dg .eq. 1) u_e(i)%val(1, 1) = u_c(i)%valdg(1, 1, 1)
    else
      if (dimensiona .eq. 2) then
        call initialise_euler2d(n, veccos, pox, poy, poz)
      else
        call initialise_euler3d(n, veccos, pox, poy, poz)
      end if
      if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
        u_c(i)%val(1, 1:nof_variables) = u_c(i)%val(1, 1:nof_variables) + veccos(1:nof_variables)*wequa3d(inc)*(voltemp)
        u_ct(i)%val(1, 1:0 + turbulenceequations + passivescalar) = u_ct(i)%val(1, 1:0 + turbulenceequations + passivescalar) + &
                                veccos(nof_variables + 1:nof_variables + turbulenceequations + passivescalar)*wequa3d(inc)*(voltemp)
      else
        if (dg .eq. 1) then
          basis_vector(1) = 1.0d0
          icompwrt = -2
          if (dimensiona .eq. 2) then
            basis_vector(2:idegfree + 1) = basis_rec2d(n, qp_array(i)%x(count_1), qp_array(i)%y(count_1), iorder, i, idegfree, icompwrt)

          else
            basis_vector(2:idegfree+1) = basis_rec(n,qp_array(i)%x(count_1),qp_array(i)%y(count_1),qp_array(i)%z(count_1),iorder,i,idegfree,icompwrt)
          end if
          do kx = 1, nof_variables
            tempsol(1, :) = veccos(kx)*wequa3d(inc)*(voltemp)*basis_vector(1:idegfree + 1)
            u_c(i)%valdg(1, kx, :) = u_c(i)%valdg(1, kx, :) + matmul(m_1(i)%val(:, :), tempsol(1, :))
          end do
        else
          u_c(i)%val(1, 1:nof_variables) = u_c(i)%val(1, 1:nof_variables) + (veccos(1:nof_variables)*wequa3d(inc)*(voltemp))
          if (itestcase .ge. 3) u_e(i)%val(1, :) = u_c(i)%val(1, :)
        end if
      end if
    end if
  end subroutine assign_initial
end module initialisation
