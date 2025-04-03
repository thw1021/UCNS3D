module fluxes
  use library
  use transform
  use local
  use riemann
  use flow_operations
  use declaration
  use basis
  use dg_functions
  implicit none

contains
  subroutine solution_integ(i, solution_integ2)
    implicit none
    real, dimension(1:nof_variables), intent(inout)::solution_integ2
    integer, intent(in)::i
    solution_integ2 = dg_vol_integral2(n, i)
  end subroutine solution_integ

  subroutine solution_integ_s(i, solution_integ_strong)
    implicit none
    real, dimension(1:nof_variables), intent(inout)::solution_integ_strong
    integer, intent(in)::i
    solution_integ_strong = dg_vol_integral_strong(n, i)
  end subroutine solution_integ_s

  subroutine solution_integ_w(i, solution_integ_weak)
    implicit none
    real, dimension(1:nof_variables), intent(inout)::solution_integ_weak
    integer, intent(in)::i
    solution_integ_weak = dg_vol_integral_weak(n, i)
  end subroutine solution_integ_w

  subroutine calculate_fluxeshi(n)
!> @brief
!> this subroutine computes the fluxes for linear-advection equation
    implicit none
    integer, intent(in)::n
    real::godflux2, sum_detect
    integer::i, l, ngp, kmaxe, iqp
    real, dimension(1:numberofpoints2)::weights_q, weights_t, weights_temp, weights_dg
    real, dimension(1:8, 1:dimensiona)::vext
    real, dimension(1:dimensiona, 1:numberofpoints2)::qpoints2d
    real, dimension(1:numberofpoints2)::wequa2d
    real::angle1, angle2, nx, ny, nz, normalvect
    integer::facex, pointx, iconsidered
    real, dimension(1:nof_variables)::cleft, cright, hllcflux, rhllcflux
    real, allocatable, dimension(:, :)::dg_rhs, dg_rhs_vol_integ, dg_rhs_surf_integ
    kmaxe = xmpielrank(n)
    call quadraturequad3d(n, igqrules, vext, qpoints2d, wequa2d)
    weights_q(1:qp_quad) = wequa2d(1:qp_quad)
    call quadraturetriang(n, igqrules, vext, qpoints2d, wequa2d)
    weights_t(1:qp_triangle) = wequa2d(1:qp_triangle)
    if (dg .eq. 1) then
      allocate(dg_rhs(1:num_dg_dofs,1:nof_variables), dg_rhs_vol_integ(1:num_dg_dofs,1:nof_variables), dg_rhs_surf_integ(1:num_dg_dofs,1:nof_variables))
    end if
    do i = 1, kmaxe
      iconsidered = i
      if (dg .eq. 1) then
        rhs(i)%valdg = zero
        dg_rhs = zero
        dg_rhs_surf_integ = zero
        dg_rhs_vol_integ = zero
      end if

      if (dg .eq. 1) then
        dg_rhs_vol_integ = dg_vol_integral(n, iconsidered)
      end if

      if (ielem(n, i)%interior .eq. 0) then
        rhs(i)%val(1) = zero
        do l = 1, ielem(n, i)%ifca
          godflux2 = zero
          angle1 = ielem(n, i)%faceanglex(l)
          angle2 = ielem(n, i)%faceangley(l)
          facex = l
          normalvect = ((cos(angle1)*sin(angle2))*lamx) + ((sin(angle1)*sin(angle2))*lamy) + ((cos(angle2))*lamz)
          if (ielem(n, i)%types_faces(l) .eq. 5) then
            iqp = qp_quad
            weights_temp(1:iqp) = weights_q(1:iqp)
          else
            iqp = qp_triangle
            weights_temp(1:iqp) = weights_t(1:iqp)
          end if

          if (dg .eq. 1) weights_dg(1:iqp) = weights_temp(1:iqp)
          do ngp = 1, iqp
            pointx = ngp
            if (dg .eq. 1) then
              cleft = ilocal_recon3(i)%uleft_dg(1:nof_variables, l, ngp)
              cright = ilocal_recon3(ielem(n, i)%ineigh(l))%uleft_dg(1:nof_variables, ielem(n, i)%ineighn(l), ngp)
            else
              cleft(1) = ilocal_recon3(i)%uleft(1, l, ngp)
              cright(1) = ilocal_recon3(ielem(n, i)%ineigh(l))%uleft(1, ielem(n, i)%ineighn(l), ngp)
            end if
            call exact_riemann_solver(n, cleft, cright, normalvect, hllcflux)
            if (dg .eq. 1) then
              rhllcflux(1) = hllcflux(1)
              dg_rhs_surf_integ = dg_rhs_surf_integ + dg_surf_flux(n, iconsidered, facex, pointx, weights_temp, rhllcflux)
            else
              godflux2 = godflux2 + (hllcflux(1)*(weights_temp(ngp)*ielem(n, i)%surf(l)))
            end if

          end do
            rhs(i)%val(1) = rhs(i)%val(1) + godflux2
        end do
      end if

      if (ielem(n, i)%interior .eq. 1) then
        rhs(i)%val(1) = zero
        do l = 1, ielem(n, i)%ifca
          facex = l
          angle1 = ielem(n, i)%faceanglex(l)
          angle2 = ielem(n, i)%faceangley(l)
          normalvect = ((cos(angle1)*sin(angle2))*lamx) + ((sin(angle1)*sin(angle2))*lamy) + ((cos(angle2))*lamz)

          if (ielem(n, i)%types_faces(l) .eq. 5) then
            iqp = qp_quad
            weights_temp(1:iqp) = weights_q(1:iqp)
          else
            iqp = qp_triangle
            weights_temp(1:iqp) = weights_t(1:iqp)
          end if

          if (dg .eq. 1) weights_dg(1:iqp) = weights_temp(1:iqp)
            godflux2 = zero
          do ngp = 1, iqp
            pointx = ngp
            if (dg .eq. 1) then
              cleft(1:nof_variables) = ilocal_recon3(i)%uleft_dg(1:nof_variables, l, ngp)
            else
              cleft(1:nof_variables) = ilocal_recon3(i)%uleft(1:nof_variables, l, ngp)
            end if
            if (ielem(n, i)%ineighb(l) .eq. n) then        !my cpu only
              if (ielem(n, i)%ibounds(l) .gt. 0) then        !check for boundaries
                if (ibound(n, ielem(n, i)%ibounds(l))%icode .eq. 5) then        !periodic in my cpu
                  if (dg .eq. 1) then
                    cright = ilocal_recon3(ielem(n, i)%ineigh(l))%uleft_dg(1:nof_variables, ielem(n, i)%ineighn(l), ngp)
                  else
                    cright(1:nof_variables) = ilocal_recon3(ielem(n, i)%ineigh(l))%uleft(1:nof_variables, ielem(n, i)%ineighn(l), ngp)
                  end if
                else
                  !not periodic ones in my cpu
                  cright(1:nof_variables) = cleft(1:nof_variables)
                end if
              else
                if (dg .eq. 1) then
                  cright(1:nof_variables) = ilocal_recon3(ielem(n, i)%ineigh(l))%uleft_dg(1:nof_variables, ielem(n, i)%ineighn(l), ngp)
                else
                  cright(1:nof_variables) = ilocal_recon3(ielem(n, i)%ineigh(l))%uleft(1:nof_variables, ielem(n, i)%ineighn(l), ngp)
                end if
              end if
            else        !in other cpus they can only be periodic or mpi neighbours
              if (ielem(n, i)%ibounds(l) .gt. 0) then        !check for boundaries
                if (ibound(n, ielem(n, i)%ibounds(l))%icode .eq. 5) then        !periodic in other cpu
                  if (dg .eq. 1) then
                    cright(1:nof_variables) = iexboundhir(ielem(n, i)%ineighn(l))%facesol_dg(ielem(n, i)%q_face(l)%q_mapl(ngp), 1:nof_variables)
                  else
                    cright(1:nof_variables) = iexboundhir(ielem(n, i)%ineighn(l))%facesol(ielem(n, i)%q_face(l)%q_mapl(ngp), 1:nof_variables)
                  end if
                end if
              else
                if (dg .eq. 1) then
                    cright(1:nof_variables) = iexboundhir(ielem(n, i)%ineighn(l))%facesol_dg(ielem(n, i)%q_face(l)%q_mapl(ngp), 1:nof_variables)
                else
                  cright(1:nof_variables) = iexboundhir(ielem(n, i)%ineighn(l))%facesol(ielem(n, i)%q_face(l)%q_mapl(ngp), 1:nof_variables)
                end if
              end if
            end if
            call exact_riemann_solver(n, cleft, cright, normalvect, hllcflux)

            if (dg .eq. 1) then
              rhllcflux(1) = hllcflux(1)
              dg_rhs_surf_integ = dg_rhs_surf_integ + dg_surf_flux(n, iconsidered, facex, pointx, weights_temp, rhllcflux)
            else
              godflux2 = godflux2 + (hllcflux(1)*(weights_temp(ngp)*ielem(n, i)%surf(l)))
            end if
          end do
            rhs(i)%val(1) = rhs(i)%val(1) + godflux2
        end do
      end if

      if (dg .eq. 1) then
        dg_rhs = dg_rhs_surf_integ - dg_rhs_vol_integ
        rhs(i)%valdg = rhs(i)%valdg + dg_rhs
      end if

    end do
    !$omp end do
    if (dg .eq. 1) then
      deallocate(dg_rhs, dg_rhs_vol_integ, dg_rhs_surf_integ)
    end if

  end subroutine calculate_fluxeshi

  subroutine calculate_fluxeshi2d(n)
!> @brief
!> this subroutine computes the fluxes for linear-advection equation in 2d
    implicit none
    integer, intent(in)::n
    real::godflux2, sum_detect, lamxl, lamyl
    integer::i, l, k, ngp, kmaxe, iqp, neighbor_index, neighbor_face_index
    real, dimension(1:numberofpoints2)::weights_temp, weights_dg !quadrature weights for interfaces
    real, dimension(1:8, 1:dimensiona)::vext
    real, dimension(1:dimensiona, 1:numberofpoints2)::qpoints2d
    real, dimension(1:numberofpoints2)::wequa2d
    real::angle1, angle2, nx, ny, nz, normalvect
    integer::facex, pointx, iconsidered
    real, dimension(1:nof_variables)::cleft, cright, hllcflux, rhllcflux
    real, allocatable, dimension(:, :)::dg_rhs, dg_rhs_vol_integ, dg_rhs_surf_integ

    if (dg .eq. 1) then
        allocate(dg_rhs(1:num_dg_dofs,1:nof_variables), dg_rhs_vol_integ(1:num_dg_dofs,1:nof_variables), dg_rhs_surf_integ(1:num_dg_dofs,1:nof_variables))
    end if
    kmaxe = xmpielrank(n)
    call quadratureline(n, igqrules, vext, qpoints2d, wequa2d)
    weights_temp(1:qp_line_n) = wequa2d(1:qp_line_n)
    !$omp do
    do i = 1, kmaxe
      if (initcond .eq. 3) then
        lamxl = -ielem(n, i)%yyc + 0.5d0
        lamyl = ielem(n, i)%xxc - 0.5
      else
        lamxl = lamx
        lamyl = lamy
      end if

      if (dg .eq. 1) then
        rhs(i)%valdg = zero
        dg_rhs = zero
        dg_rhs_surf_integ = zero
        dg_rhs_vol_integ = zero
      else
        rhs(i)%val = zero
      end if

      iconsidered = i

      if (dg .eq. 1) then
        dg_rhs_vol_integ = dg_vol_integral(n, iconsidered)
      end if

      if (ielem(n, i)%interior .eq. 0) then ! element is interior
        do l = 1, ielem(n, i)%ifca
          godflux2 = zero
          nx = ielem(n, i)%faceanglex(l)
          ny = ielem(n, i)%faceangley(l)
          facex = l
          normalvect = (nx*lamxl) + (ny*lamyl)
          iqp = qp_line_n
          neighbor_index = ielem(n, i)%ineigh(l)
          neighbor_face_index = ielem(n, i)%ineighn(l)

          do ngp = 1, iqp
            pointx = ngp

            if (dg .eq. 1) then
              cleft = ilocal_recon3(i)%uleft_dg(1:nof_variables, l, ngp)
              cright = ilocal_recon3(ielem(n, i)%ineigh(l))%uleft_dg(1:nof_variables, ielem(n, i)%ineighn(l), ngp)
            else !fv
              cleft(1) = ilocal_recon3(i)%uleft(1, l, ngp)
              cright(1) = ilocal_recon3(ielem(n, i)%ineigh(l))%uleft(1, ielem(n, i)%ineighn(l), ngp)
            end if
            call exact_riemann_solver(n, cleft, cright, normalvect, hllcflux)
            if (dg .eq. 1) then
              ! riemann flux at interface quadrature points times basis
              rhllcflux(1) = hllcflux(1)
              dg_rhs_surf_integ = dg_rhs_surf_integ + dg_surf_flux(n, iconsidered, facex, pointx, weights_temp, rhllcflux)
            else !fv
              godflux2 = godflux2 + (hllcflux(1)*(weights_temp(ngp)*ielem(n, i)%surf(l)))
            end if
          end do
          if (dg .ne. 1) rhs(i)%val(1) = rhs(i)%val(1) + godflux2
        end do

      else if (ielem(n, i)%interior .eq. 1) then

        do l = 1, ielem(n, i)%ifca
          facex = l
          nx = ielem(n, i)%faceanglex(l)
          ny = ielem(n, i)%faceangley(l)
          normalvect = (nx*lamxl) + (ny*lamyl)
          iqp = qp_line_n
          godflux2 = zero
          do ngp = 1, iqp
            pointx = ngp
            if (dg .eq. 1) then
              cleft = ilocal_recon3(i)%uleft_dg(1:nof_variables, l, ngp)
            else
              cleft(1) = ilocal_recon3(i)%uleft(1, l, ngp)
            end if

            if (ielem(n, i)%ineighb(l) .eq. n) then        !my cpu only
              if (ielem(n, i)%ibounds(l) .gt. 0) then        !check for boundaries
                if (ibound(n, ielem(n, i)%ibounds(l))%icode .eq. 5) then        !periodic in my cpu
                  if (dg .eq. 1) then
                    cright = ilocal_recon3(ielem(n, i)%ineigh(l))%uleft_dg(1:nof_variables, ielem(n, i)%ineighn(l), ngp)
                  else !fv
                    cright(1) = ilocal_recon3(ielem(n, i)%ineigh(l))%uleft(1, ielem(n, i)%ineighn(l), ngp)
                  end if
                else !not periodic ones in my cpu
                  cright(1:nof_variables) = cleft(1:nof_variables)
                end if
              else
                if (dg .eq. 1) then
                  cright = ilocal_recon3(ielem(n, i)%ineigh(l))%uleft_dg(1:nof_variables, ielem(n, i)%ineighn(l), ngp)
                else !fv
                  cright(1) = ilocal_recon3(ielem(n, i)%ineigh(l))%uleft(1, ielem(n, i)%ineighn(l), ngp)
                end if
              end if
            else !in other cpus they can only be periodic or mpi neighbours
              if (ielem(n, i)%ibounds(l) .gt. 0) then        !check for boundaries
                if (ibound(n, ielem(n, i)%ibounds(l))%icode .eq. 5) then        !periodic in other cpu
                  if (dg .eq. 1) then
                    cright(1:nof_variables) = iexboundhir(ielem(n, i)%ineighn(l))%facesol_dg(ielem(n, i)%q_face(l)%q_mapl(ngp), 1:nof_variables)
                  else
                    cright(1:nof_variables) = iexboundhir(ielem(n, i)%ineighn(l))%facesol(ielem(n, i)%q_face(l)%q_mapl(ngp), 1:nof_variables)
                  end if
                end if
              else
                if (dg .eq. 1) then
                   cright(1:nof_variables) = iexboundhir(ielem(n, i)%ineighn(l))%facesol_dg(ielem(n, i)%q_face(l)%q_mapl(ngp), 1:nof_variables)
                else
                   cright(1:nof_variables) = iexboundhir(ielem(n, i)%ineighn(l))%facesol(ielem(n, i)%q_face(l)%q_mapl(ngp), 1:nof_variables)
                end if
              end if
            end if
            call exact_riemann_solver(n, cleft, cright, normalvect, hllcflux)
            if (dg .eq. 1) then
              !riemann flux at interface quadrature points times basis
              rhllcflux(1) = hllcflux(1)
              dg_rhs_surf_integ = dg_rhs_surf_integ + dg_surf_flux(n, iconsidered, facex, pointx, weights_temp, rhllcflux)
            else !fv
              godflux2 = godflux2 + (hllcflux(1)*(weights_temp(ngp)*ielem(n, i)%surf(l)))
            end if
          end do
          if (dg .ne. 1) rhs(i)%val(1) = rhs(i)%val(1) + godflux2
        end do
      end if
      if (dg .eq. 1) dg_rhs = dg_rhs_surf_integ - dg_rhs_vol_integ
      if (dg .eq. 1) rhs(i)%valdg = rhs(i)%valdg + dg_rhs

    end do
    !$omp end do

    if (dg .eq. 1) then
      deallocate(dg_rhs, dg_rhs_vol_integ, dg_rhs_surf_integ)
    end if

  end subroutine calculate_fluxeshi2d

  subroutine calculate_fluxeshi_convective(n)
!> @brief
!> this subroutine computes the convective fluxes for hyperbolic conservation laws
    implicit none
    integer, intent(in)::n
    real, dimension(1:nof_variables + turbulenceequations + passivescalar)::godflux2, dg_vol_rec, rhllcflux, hllcflux
    integer::i, l, ngp, kmaxe, iqp, ii, ikas, igoflux, icaseb, jx, jx2, b_code, srf
    real::sum_detect, norms, tempxx
    real, dimension(1:numberofpoints2)::weights_q, weights_t, weights_temp, weights_dg
    integer::iconsidered, facex, pointx
    real::angle1, angle2, nx, ny, nz, mp_source1, mp_source2, mp_source3
    real, dimension(1:nof_variables + turbulenceequations + passivescalar)::cleft, cright, cleft_rot, cright_rot
    real, dimension(1:nof_variables)::leftv, rightv, srf_speedrot
    real, dimension(1:turbulenceequations + passivescalar)::cturbl, cturbr
    real, dimension(1:dimensiona)::pox, poy, poz
    real, dimension(1:nof_variables)::srf_speed
    real, dimension(1:8, 1:dimensiona)::vext
    real, dimension(1:dimensiona, 1:numberofpoints2)::qpoints2d
    real, dimension(1:numberofpoints2)::wequa2d
    real::x1, y1, z1
    real, allocatable, dimension(:, :)::dg_rhs, dg_rhs_vol_integ, dg_rhs_surf_integ

    if (dg .eq. 1) then
        allocate(dg_rhs(1:num_dg_dofs,1:nof_variables), dg_rhs_vol_integ(1:num_dg_dofs,1:nof_variables), dg_rhs_surf_integ(1:num_dg_dofs,1:nof_variables))
    end if

    kmaxe = xmpielrank(n)

    call quadraturequad3d(n, igqrules, vext, qpoints2d, wequa2d)

    weights_q(1:qp_quad) = wequa2d(1:qp_quad)

    call quadraturetriang(n, igqrules, vext, qpoints2d, wequa2d)
    weights_t(1:qp_triangle) = wequa2d(1:qp_triangle)

    do i = 1, xmpielrank(n)
      rhs(i)%val(:) = zero; if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) rhst(i)%val(:) = zero

    end do
    !$omp barrier
    !$omp do
    do ii = 1, nof_interior        !for all the interior elements
      i = el_int(ii)
      iconsidered = i
      mp_source3 = zero

      if (mrf .eq. 1) then
        srf = ilocal_recon3(i)%mrf
      end if
      if (dg .eq. 1) then
        rhs(i)%valdg = zero
        dg_rhs = zero
        dg_rhs_surf_integ = zero
        dg_rhs_vol_integ = zero

        dg_rhs_vol_integ = dg_vol_integral(n, iconsidered)

      end if

      do l = 1, ielem(n, i)%ifca !for all their faces
        b_code = 0
        godflux2 = zero
        mp_source2 = zero
        angle1 = ielem(n, i)%faceanglex(l)
        angle2 = ielem(n, i)%faceangley(l)
        nx = (cos(angle1)*sin(angle2))
        ny = (sin(angle1)*sin(angle2))
        nz = (cos(angle2))
        if (ielem(n, i)%types_faces(l) .eq. 5) then
          iqp = qp_quad_n
          weights_temp(1:iqp) = weights_q(1:iqp)
        else
          iqp = qp_triangle_n
          weights_temp(1:iqp) = weights_t(1:iqp)
        end if

        facex = l
        if (dg .eq. 1) weights_dg(1:iqp) = weights_temp(1:iqp)

        do ngp = 1, iqp        !for all the gaussian quadrature points
          pointx = ngp
          call get_states_interior(n,b_code,iconsidered,facex,pointx,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speedrot,cleft,cright)

          call rotatef(n, cright_rot, cright, angle1, angle2)        !rotate wrt to normalvector of face and solve 1d riemann problem
          call rotatef(n, cleft_rot, cleft, angle1, angle2)        !rotate wrt to normalvector of face and solve 1d riemann problem

          if ((lmach .eq. 1)) then    !application of the low mach number correction
            leftv(1:nof_variables) = cleft_rot(1:nof_variables); rightv(1:nof_variables) = cright_rot(1:nof_variables)
            call lmacht(n, leftv, rightv)
            cleft_rot(1:nof_variables) = leftv(1:nof_variables); cright_rot(1:nof_variables) = rightv(1:nof_variables);
          end if

          select case (iriemann)
          case (1)                        !hllc
            call hllc_riemann_solver(n, iconsidered, facex, cleft_rot, cright_rot, hllcflux, mp_source1, srf_speedrot)
            call rotateb(n, rhllcflux, hllcflux, angle1, angle2)
            if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
               rhllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=hllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)
            end if

          case (9)                        !hll
            call hll_riemann_solver(n, cleft_rot, cright_rot, hllcflux, mp_source1, srf_speedrot)
            call rotateb(n, rhllcflux, hllcflux, angle1, angle2)
            if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
               rhllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=hllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)
            end if

          case (2)                        !rusanov
            call rusanov_riemann_solver(n, iconsidered, facex, cleft_rot, cright_rot, hllcflux, mp_source1, srf_speedrot)
            call rotateb(n, rhllcflux, hllcflux, angle1, angle2)
            if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
               rhllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=hllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)
            end if

          case (3)                        !roe
            call rotateb(n, cleft, cleft_rot, angle1, angle2)
            call rotateb(n, cright, cright_rot, angle1, angle2)
            call roe_riemann_solver(n, iconsidered, facex, cleft, cright, hllcflux, mp_source1, srf_speedrot, nx, ny, nz)
            rhllcflux = hllcflux

          case (4)                        !roe
            call rotateb(n, cleft, cleft_rot, angle1, angle2)
            call rotateb(n, cright, cright_rot, angle1, angle2)
            call rroe_riemann_solver(n, iconsidered, facex, cleft, cright, hllcflux, mp_source1, srf_speedrot, nx, ny, nz)
            rhllcflux = hllcflux

          case (5)                        !roe
            call rotateb(n, cleft, cleft_rot, angle1, angle2)
            call rotateb(n, cright, cright_rot, angle1, angle2)
            call troe_riemann_solver(n, cleft, cright, hllcflux, mp_source1, srf_speedrot, nx, ny, nz)
            rhllcflux = hllcflux

          end select

          if (dg .eq. 1) then
            dg_rhs_surf_integ = dg_rhs_surf_integ + dg_surf_flux(n, iconsidered, facex, pointx, weights_temp, rhllcflux)
          else
            godflux2(1:nof_variables) = godflux2(1:nof_variables) + (rhllcflux(1:nof_variables)*(weights_temp(ngp)*ielem(n, i)%surf(l)))
          end if

          if (multispecies .eq. 1) then
            mp_source2 = mp_source2 + mp_source1*(weights_temp(ngp)*ielem(n, i)%surf(l))
          end if

          if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
            if (icoupleturb .eq. 0) then        !first order upwind flux
              norms = 0.5*(cleft_rot(2) + cright_rot(2))
              if (ilocal_recon3(iconsidered)%mrf .eq. 1) then
                norms = norms - srf_speedrot(2)
              end if
                rhllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=0.5*((norms*(cturbl(:)+cturbr(:)))+(abs(norms)*(cturbl(:)-(cturbr(:)))))
            end if

           godflux2(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=godflux2(nof_variables+1:nof_variables+turbulenceequations+passivescalar)+&
          (rhllcflux(nof_variables + 1:nof_variables + turbulenceequations + passivescalar)*(weights_temp(ngp)*ielem(n, i)%surf(l)))
          end if
        end do

        rhs(i)%val(1:nof_variables) = rhs(i)%val(1:nof_variables) + godflux2(1:nof_variables)
        if (multispecies .eq. 1) then
          mp_source3 = mp_source3 + mp_source2
        end if
          if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
            rhst(i)%val(1:turbulenceequations + passivescalar) = rhst(i)%val(1:turbulenceequations + passivescalar) + &
            godflux2(nof_variables + 1:nof_variables + turbulenceequations + passivescalar)
        end if
      end do

      if (multispecies .eq. 1) then
        rhs(i)%val(8) = rhs(i)%val(8) - (u_c(i)%val(1, 8)*mp_source3)
      end if

      if (dg .eq. 1) then
        dg_rhs = dg_rhs_surf_integ - dg_rhs_vol_integ
        rhs(i)%valdg = rhs(i)%valdg + dg_rhs

        if (multispecies .eq. 1) then
          rhs(i)%valdg(1, 8) = rhs(i)%valdg(1, 8) - (u_c(i)%val(1, 8)*mp_source3)
        end if
      end if

    end do
    !$omp end do
    !$omp barrier
    !$omp do
    do ii = 1, nof_bounded
      i = el_bnd(ii)
      iconsidered = i
      mp_source3 = zero

      if (mrf .eq. 1) then
        srf = ilocal_recon3(i)%mrf
      end if
      if (dg .eq. 1) then
        rhs(i)%valdg = zero
        dg_rhs = zero
        dg_rhs_surf_integ = zero
        dg_rhs_vol_integ = zero
        dg_rhs_vol_integ = dg_vol_integral(n, iconsidered)
      end if

      do l = 1, ielem(n, i)%ifca
        facex = l
        angle1 = ielem(n, i)%faceanglex(l)
        angle2 = ielem(n, i)%faceangley(l)
        nx = (cos(angle1)*sin(angle2))
        ny = (sin(angle1)*sin(angle2))
        nz = (cos(angle2))

        if (ielem(n, i)%types_faces(l) .eq. 5) then
          iqp = qp_quad_n
          weights_temp(1:iqp) = weights_q(1:iqp)
        else
          iqp = qp_triangle_n
          weights_temp(1:iqp) = weights_t(1:iqp)
        end if

        if (dg .eq. 1) weights_dg(1:iqp) = weights_temp(1:iqp)
        godflux2 = zero
        mp_source2 = zero
        do ngp = 1, iqp
          pointx = ngp
          b_code = 0

          call get_states_bounds(n,b_code,iconsidered,facex,pointx,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speedrot,cleft,cright)
          call rotatef(n, cright_rot, cright, angle1, angle2)        !rotate wrt to normalvector of face and solve 1d riemann problem
          call rotatef(n, cleft_rot, cleft, angle1, angle2)        !rotate wrt to normalvector of face and solve 1d riemann problem

          if ((lmach .eq. 1)) then    !application of the low mach number correction
            leftv(1:nof_variables) = cleft_rot(1:nof_variables); rightv(1:nof_variables) = cright_rot(1:nof_variables)
            call lmacht(n, leftv, rightv)
            cleft_rot(1:nof_variables) = leftv(1:nof_variables); cright_rot(1:nof_variables) = rightv(1:nof_variables);
          end if

          if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
            cleft_rot(nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = cturbl(1:turbulenceequations + passivescalar)
            cright_rot(nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = cturbr(1:turbulenceequations + passivescalar)
          end if

          select case (iriemann)

          case (1)                        !hllc
            call hllc_riemann_solver(n, iconsidered, facex, cleft_rot, cright_rot, hllcflux, mp_source1, srf_speedrot)
            call rotateb(n, rhllcflux, hllcflux, angle1, angle2)
            if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
              rhllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=hllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)
            end if

          case (9)                        !hllc
            call hll_riemann_solver(n, cleft_rot, cright_rot, hllcflux, mp_source1, srf_speedrot)
            call rotateb(n, rhllcflux, hllcflux, angle1, angle2)
            if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
              rhllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=hllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)
            end if

          case (2)                        !rusanov
            call rusanov_riemann_solver(n, iconsidered, facex, cleft_rot, cright_rot, hllcflux, mp_source1, srf_speedrot)
            call rotateb(n, rhllcflux, hllcflux, angle1, angle2)
            if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
              rhllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=hllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)
            end if

          case (3)                        !roe
            call rotateb(n, cleft, cleft_rot, angle1, angle2)
            call rotateb(n, cright, cright_rot, angle1, angle2)
            call roe_riemann_solver(n, iconsidered, facex, cleft, cright, hllcflux, mp_source1, srf_speedrot, nx, ny, nz)
            rhllcflux = hllcflux

          case (4)                !roe
            if (b_code .gt. 0) then
              call rusanov_riemann_solver(n, iconsidered, facex, cleft_rot, cright_rot, hllcflux, mp_source1, srf_speedrot)
              call rotateb(n, rhllcflux, hllcflux, angle1, angle2)
            else
              call rotateb(n, cleft, cleft_rot, angle1, angle2)
              call rotateb(n, cright, cright_rot, angle1, angle2)
              call rroe_riemann_solver(n, iconsidered, facex, cleft, cright, hllcflux, mp_source1, srf_speedrot, nx, ny, nz)
              rhllcflux = hllcflux
            end if

          case (5)                        !roe
            call rotateb(n, cleft, cleft_rot, angle1, angle2)
            call rotateb(n, cright, cright_rot, angle1, angle2)
            if (b_code .le. 0) then
              call troe_riemann_solver(n, cleft, cright, hllcflux, mp_source1, srf_speedrot, nx, ny, nz)
            else
              call roe_riemann_solver(n, iconsidered, facex, cleft, cright, hllcflux, mp_source1, srf_speedrot, nx, ny, nz)
            end if

            rhllcflux = hllcflux

          end select
          if (dg .eq. 1) then
            dg_rhs_surf_integ = dg_rhs_surf_integ + dg_surf_flux(n, iconsidered, facex, pointx, weights_temp, rhllcflux)
          else
            godflux2(1:nof_variables) = godflux2(1:nof_variables) + (rhllcflux(1:nof_variables)*(weights_temp(ngp)*ielem(n, i)%surf(l)))
          end if
          if (multispecies .eq. 1) then
            mp_source2 = mp_source2 + mp_source1*(weights_temp(ngp)*ielem(n, i)%surf(l))
          end if
            if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
              if (icoupleturb .eq. 0) then        !first order upwind flux
                norms = 0.5*(cleft_rot(2) + cright_rot(2))
              if (ilocal_recon3(iconsidered)%mrf .eq. 1) then
                norms = norms - srf_speedrot(2)
              end if
                 rhllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=0.5*((norms*(cturbl(:)+cturbr(:)))+(abs(norms)*(cturbl(:)-(cturbr(:)))))
            end if

            if ((b_code .eq. 3) .or. (b_code .eq. 4)) then
              rhllcflux(nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = zero
            end if
              godflux2(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=godflux2(nof_variables+1:nof_variables+turbulenceequations+passivescalar)+&
             (rhllcflux(nof_variables + 1:nof_variables + turbulenceequations + passivescalar)*(weights_temp(ngp)*ielem(n, i)%surf(l)))
          end if
        end do

        rhs(i)%val(1:nof_variables) = rhs(i)%val(1:nof_variables) + godflux2(1:nof_variables)
        if (multispecies .eq. 1) then
          mp_source3 = mp_source3 + mp_source2
        end if

        if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
          rhst(i)%val(1:turbulenceequations + passivescalar) = rhst(i)%val(1:turbulenceequations + passivescalar) + &
          godflux2(nof_variables + 1:nof_variables + turbulenceequations + passivescalar)
        end if
      end do
      if (multispecies .eq. 1) then
        rhs(i)%val(8) = rhs(i)%val(8) - (u_c(i)%val(1, 8)*mp_source3)
      end if

      if (dg .eq. 1) then
        dg_rhs = dg_rhs_surf_integ - dg_rhs_vol_integ
        rhs(i)%valdg = rhs(i)%valdg + dg_rhs
        if (multispecies .eq. 1) then
          rhs(i)%valdg(1, 8) = rhs(i)%valdg(1, 8) - (u_c(i)%val(1, 8)*mp_source3)
        end if
      end if
    end do

    if (dg .eq. 1) then
      deallocate(dg_rhs, dg_rhs_vol_integ, dg_rhs_surf_integ)
    end if

  end subroutine calculate_fluxeshi_convective

  subroutine calculate_fluxeshi_convective2d(n)
!> @brief
!> this subroutine computes the convective fluxes for hyperbolic conservation laws in 2d
    implicit none
    integer, intent(in)::n
    real, dimension(1:nof_variables + turbulenceequations + passivescalar)::godflux2, dg_vol_rec, rhllcflux, hllcflux
    integer::i, l, ngp, kmaxe, iqp, ii, ikas, igoflux, icaseb, kxk, b_code
    real::sum_detect, norms
    real, dimension(1:numberofpoints2)::weights_temp, weights_dg
    integer::iconsidered, facex, pointx
    real::angle1, angle2, nx, ny, nz, mp_source1, mp_source2, mp_source3
    real, dimension(1:nof_variables + turbulenceequations + passivescalar)::cleft, cright, cleft_rot, cright_rot
    real, dimension(1:nof_variables)::leftv, rightv, srf_speedrot
    real, dimension(1:turbulenceequations + passivescalar)::cturbl, cturbr
    real, dimension(1:dimensiona)::pox, poy, poz
    real::x1, y1, z1
    real, dimension(1:nof_variables)::srf_speed
    real, dimension(1:8, 1:dimensiona)::vext
    real, dimension(1:dimensiona, 1:numberofpoints2)::qpoints2d
    real, dimension(1:numberofpoints2)::wequa2d
    real, allocatable, dimension(:, :)::dg_rhs, dg_rhs_vol_integ, dg_rhs_surf_integ

    if (dg .eq. 1) then
        allocate(dg_rhs(1:num_dg_dofs,1:nof_variables), dg_rhs_vol_integ(1:num_dg_dofs,1:nof_variables), dg_rhs_surf_integ(1:num_dg_dofs,1:nof_variables))
    end if
    kmaxe = xmpielrank(n)

    call quadratureline(n, igqrules, vext, qpoints2d, wequa2d)
    weights_temp = wequa2d(1:qp_line_n)

    if (reduce_comp .eq. 1) then
      weights_temp = 1.0d0
    end if

    do i = 1, kmaxe
      rhs(i)%val(:) = zero; if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) rhst(i)%val(:) = zero
      iconsidered = i
    end do

    !$omp barrier
    !$omp do
    do ii = 1, nof_interior        !for all the interior elements
      i = el_int(ii)
      iconsidered = i
      if (dg .eq. 1) then
        rhs(i)%valdg = zero
        dg_rhs = zero
        dg_rhs_surf_integ = zero
        dg_rhs_vol_integ = zero
        dg_rhs_vol_integ = dg_vol_integral(n, iconsidered)
      end if

      mp_source3 = zero
      do l = 1, ielem(n, i)%ifca !for all their faces

        godflux2 = zero
        mp_source2 = zero
        nx = ielem(n, i)%faceanglex(l)
        ny = ielem(n, i)%faceangley(l)
        angle1 = nx
        angle2 = ny
        b_code = 0
        iqp = qp_line_n
        do ngp = 1, iqp        !for all the gaussian quadrature points
          pointx = ngp
          facex = l
          pointx = ngp

          call get_states_interior2d(n,b_code,iconsidered,facex,pointx,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speedrot,cleft,cright)
          call rotatef2d(n, cright_rot, cright, angle1, angle2)        !rotate wrt to normalvector of face and solve 1d riemann problem
          call rotatef2d(n, cleft_rot, cleft, angle1, angle2)        !rotate wrt to normalvector of face and solve 1d riemann problem

          if ((lmach .eq. 1)) then    !application of the low mach number correction
            leftv(1:nof_variables) = cleft_rot(1:nof_variables); rightv(1:nof_variables) = cright_rot(1:nof_variables)
            call lmacht2d(n, leftv, rightv)
            cleft_rot(1:nof_variables) = leftv(1:nof_variables); cright_rot(1:nof_variables) = rightv(1:nof_variables);
          end if

          select case (iriemann)
          case (1)                        !hllc
            call hllc_riemann_solver2d(n, cleft_rot, cright_rot, hllcflux, mp_source1, srf_speedrot)
            call rotateb2d(n, rhllcflux, hllcflux, angle1, angle2)
            if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
            rhllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=hllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)
            end if

          case (9)                        !hll
            call hll_riemann_solver2d(n, cleft_rot, cright_rot, hllcflux, mp_source1, srf_speedrot)
            call rotateb2d(n, rhllcflux, hllcflux, angle1, angle2)
            if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
              rhllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=hllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)
            end if

          case (2)                        !rusanov
            call rusanov_riemann_solver2d(n, cleft_rot, cright_rot, hllcflux, mp_source1, srf_speedrot)
            call rotateb2d(n, rhllcflux, hllcflux, angle1, angle2)
            if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
              rhllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=hllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)
            end if

          case (3)                        !roe
            call rotateb2d(n, cleft, cleft_rot, angle1, angle2)
            call rotateb2d(n, cright, cright_rot, angle1, angle2)
            call roe_riemann_solver2d(n, cleft, cright, hllcflux, mp_source1, srf_speedrot, nx, ny)
            rhllcflux(1:nof_variables) = hllcflux(1:nof_variables)

          case (4)                        !roe
            call rotateb2d(n, cleft, cleft_rot, angle1, angle2)
            call rotateb2d(n, cright, cright_rot, angle1, angle2)
            call rroe_riemann_solver2d(n, cleft, cright, hllcflux, mp_source1, srf_speedrot, nx, ny, b_code)
            rhllcflux = hllcflux

          end select

          if (dg .eq. 1) then
            ! riemann flux at interface quadrature points times basis
            dg_rhs_surf_integ = dg_rhs_surf_integ + dg_surf_flux(n, iconsidered, facex, pointx, weights_temp, rhllcflux)
          else !fv
            godflux2(1:nof_variables) = godflux2(1:nof_variables) + (rhllcflux(1:nof_variables)*(weights_temp(ngp)*ielem(n, i)%surf(l)))
          end if

          if (multispecies .eq. 1) then
            mp_source2 = mp_source2 + mp_source1*(weights_temp(ngp)*ielem(n, i)%surf(l))
          end if
          if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
            if (icoupleturb .eq. 0) then        !first order upwind flux
              norms = 0.5*(cleft_rot(2) + cright_rot(2))
                rhllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=0.5*((norms*(cturbl(:)+cturbr(:)))+(abs(norms)*(cturbl(:)-(cturbr(:)))))
            end if
              godflux2(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=godflux2(nof_variables+1:nof_variables+turbulenceequations+passivescalar)+&
             (rhllcflux(nof_variables + 1:nof_variables + turbulenceequations + passivescalar)*(weights_temp(ngp)*ielem(n, i)%surf(l)))
          end if

        end do
        rhs(i)%val(1:nof_variables) = rhs(i)%val(1:nof_variables) + godflux2(1:nof_variables)
        if (multispecies .eq. 1) then
          mp_source3 = mp_source3 + mp_source2
        end if
        if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
          rhst(i)%val(1:turbulenceequations + passivescalar) = rhst(i)%val(1:turbulenceequations + passivescalar) + &
          godflux2(nof_variables + 1:nof_variables + turbulenceequations + passivescalar)
        end if
      end do
      if (multispecies .eq. 1) then
        rhs(i)%val(7) = rhs(i)%val(7) - (u_c(i)%val(1, 7)*mp_source3)!*ielem(n,i)%totvolume)
      end if

      if (dg .eq. 1) then
        dg_rhs = dg_rhs_surf_integ - dg_rhs_vol_integ
        rhs(i)%valdg = rhs(i)%valdg + dg_rhs
        if (multispecies .eq. 1) then
          rhs(i)%valdg(1, 7) = rhs(i)%valdg(1, 7) - (u_c(i)%val(1, 7)*mp_source3)
        end if
      end if

    end do
    !$omp end do
    !$omp do
    do ii = 1, nof_bounded
      i = el_bnd(ii)
      iconsidered = i
      mp_source3 = zero

      if (dg .eq. 1) then
        rhs(i)%valdg = zero
        dg_rhs = zero
        dg_rhs_surf_integ = zero
        dg_rhs_vol_integ = zero
        dg_rhs_vol_integ = dg_vol_integral(n, iconsidered)
      end if

      do l = 1, ielem(n, i)%ifca
        facex = l
        igoflux = 0
        b_code = 0
        nx = ielem(n, i)%faceanglex(l)
        ny = ielem(n, i)%faceangley(l)
        angle1 = nx
        angle2 = ny
        iqp = qp_line_n
        godflux2 = zero
        mp_source2 = zero
        do ngp = 1, iqp
          pointx = ngp
          facex = l

          call get_states_bounds2d(n,b_code,iconsidered,facex,pointx,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speedrot,cleft,cright)
          call rotatef2d(n, cright_rot, cright, angle1, angle2)        !rotate wrt to normalvector of face and solve 1d riemann problem
          call rotatef2d(n, cleft_rot, cleft, angle1, angle2)        !rotate wrt to normalvector of face and solve 1d riemann problem

          if ((lmach .eq. 1)) then    !application of the low mach number correction
            leftv(1:nof_variables) = cleft_rot(1:nof_variables); rightv(1:nof_variables) = cright_rot(1:nof_variables)
            call lmacht2d(n, leftv, rightv)
            cleft_rot(1:nof_variables) = leftv(1:nof_variables); cright_rot(1:nof_variables) = rightv(1:nof_variables);
          end if

          if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
            cleft_rot(nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = cturbl(1:turbulenceequations + passivescalar)
            cright_rot(nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = cturbr(1:turbulenceequations + passivescalar)
          end if

          select case (iriemann)
          case (1)                        !hllc
            call hllc_riemann_solver2d(n, cleft_rot, cright_rot, hllcflux, mp_source1, srf_speedrot)
            call rotateb2d(n, rhllcflux, hllcflux, angle1, angle2)
            if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
              rhllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=hllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)
            end if

          case (9)                        !hll
            call hll_riemann_solver2d(n, cleft_rot, cright_rot, hllcflux, mp_source1, srf_speedrot)
            call rotateb2d(n, rhllcflux, hllcflux, angle1, angle2)
            if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
              rhllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=hllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)
            end if

          case (2)                        !rusanov
            call rusanov_riemann_solver2d(n, cleft_rot, cright_rot, hllcflux, mp_source1, srf_speedrot)
            call rotateb2d(n, rhllcflux, hllcflux, angle1, angle2)
            if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
              rhllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=hllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)
            end if

          case (3)                        !roe
            call rotateb2d(n, cleft, cleft_rot, angle1, angle2)
            call rotateb2d(n, cright, cright_rot, angle1, angle2)
            call roe_riemann_solver2d(n, cleft, cright, hllcflux, mp_source1, srf_speedrot, nx, ny)
            rhllcflux = hllcflux

          case (4)                        !roe
            if ((b_code .le. 0)) then
              call rotateb2d(n, cleft, cleft_rot, angle1, angle2)
              call rotateb2d(n, cright, cright_rot, angle1, angle2)
              call rroe_riemann_solver2d(n, cleft, cright, hllcflux, mp_source1, srf_speedrot, nx, ny, b_code)
              rhllcflux = hllcflux
            else
              call rusanov_riemann_solver2d(n, cleft_rot, cright_rot, hllcflux, mp_source1, srf_speedrot)
              call rotateb2d(n, rhllcflux, hllcflux, angle1, angle2)
            end if

          end select

          if (dg .eq. 1) then
            dg_rhs_surf_integ = dg_rhs_surf_integ + dg_surf_flux(n, iconsidered, facex, pointx, weights_temp, rhllcflux)
          else

        godflux2(1:nof_variables) = godflux2(1:nof_variables) + (rhllcflux(1:nof_variables)*(weights_temp(ngp)*ielem(n, i)%surf(l)))
          end if

          if (multispecies .eq. 1) then
            mp_source2 = mp_source2 + mp_source1*(weights_temp(ngp)*ielem(n, i)%surf(l))
          end if

          if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
            if (icoupleturb .eq. 0) then        !first order upwind flux

              norms = 0.5*(cleft_rot(2) + cright_rot(2))
                rhllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=0.5*((norms*(cturbl(:)+cturbr(:)))+(abs(norms)*(cturbl(:)-(cturbr(:)))))
            end if

            if ((b_code .eq. 4) .or. (b_code .eq. 3)) then
              rhllcflux(nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = zero
            end if

            godflux2(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=godflux2(nof_variables+1:nof_variables+turbulenceequations+passivescalar)+&
            (rhllcflux(nof_variables + 1:nof_variables + turbulenceequations + passivescalar)*(weights_temp(ngp)*ielem(n, i)%surf(l)))

          end if
        end do

        rhs(i)%val(1:nof_variables) = rhs(i)%val(1:nof_variables) + godflux2(1:nof_variables)
        if (multispecies .eq. 1) then
          mp_source3 = mp_source3 + mp_source2
        end if

        if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
          rhst(i)%val(1:turbulenceequations + passivescalar) = rhst(i)%val(1:turbulenceequations + passivescalar) + &
          godflux2(nof_variables + 1:nof_variables + turbulenceequations + passivescalar)
        end if
      end do
      if (multispecies .eq. 1) then
        rhs(i)%val(7) = rhs(i)%val(7) - (u_c(i)%val(1, 7)*mp_source3)!*ielem(n,i)%totvolume)
      end if

      if (dg .eq. 1) then
        dg_rhs = dg_rhs_surf_integ - dg_rhs_vol_integ
        rhs(i)%valdg = rhs(i)%valdg + dg_rhs
        if (multispecies .eq. 1) then
          rhs(i)%valdg(1, 7) = rhs(i)%valdg(1, 7) - (u_c(i)%val(1, 7)*mp_source3)
        end if
      end if

    end do
    !$omp end do
    if (dg .eq. 1) then
      deallocate(dg_rhs, dg_rhs_vol_integ, dg_rhs_surf_integ)
    end if

  end subroutine calculate_fluxeshi_convective2d

  subroutine calculate_fluxeshi_diffusive(n)
!> @brief
!> this subroutine computes the diffusive fluxes for euler-navier-stokes equations
    implicit none
    integer, intent(in)::n
    real, dimension(1:nof_variables + turbulenceequations + passivescalar)::godflux2, dg_vol_rec, rhllcflux, hllcflux
    integer::i, l, ngp, kmaxe, iqp, ii, nvar, kc, iex, ittt, ikas, igoflux, icaseb, kk, b_code, srf
    real::sum_detect, norms
    integer::iconsidered, facex, pointx
    real::angle1, angle2, nx, ny, nz, mp_source1, mp_source2, mp_source3
    real, dimension(1:nof_variables + turbulenceequations + passivescalar)::cleft, cright, cleft_rot, cright_rot
    real, dimension(1:nof_variables)::leftv, rightv, srf_speedrot
    real, dimension(1:turbulenceequations + passivescalar)::cturbl, cturbr
    real, dimension(1:dimensiona)::pox, poy, poz
    real, dimension(1:nof_variables)::srf_speed
    real, dimension(1:numberofpoints2)::weights_q, weights_t, weights_temp, weights_dg
    real, dimension(1:8, 1:dimensiona)::vext
    real, dimension(1:dimensiona, 1:numberofpoints2)::qpoints2d
    real, dimension(1:numberofpoints2)::wequa2d
    real::mp_pinfl, mp_pinfr, gammal, gammar
    real, dimension(1:4)::viscl, laml
    real, dimension(1:2)::turbmv
    real, dimension(1)::etvm
    real, dimension(1:20)::eddyfl, eddyfr
    real, dimension(1:nof_variables - 1, 1:dims)::lcvgrad, rcvgrad
    real, dimension(turbulenceequations + passivescalar, 1:dims)::lcvgrad_t, rcvgrad_t
    real, dimension(5)::fxv, fyv, fzv, tem_pn, rtem_pn
    real, dimension(3, 3)::taul, taur, tau
    real, dimension(3)::q, nnn, nall
    real::ux, uy, uz, vx, vy, vz, wx, wy, wz, rho12, u12, v12, w12, damp, vdamp, tempxx
    real, allocatable, dimension(:, :)::dg_rhs, dg_rhs_vol_integ, dg_rhs_surf_integ

    if (dg .eq. 1) then
        allocate(dg_rhs(1:num_dg_dofs,1:nof_variables), dg_rhs_vol_integ(1:num_dg_dofs,1:nof_variables), dg_rhs_surf_integ(1:num_dg_dofs,1:nof_variables))
    end if

    kmaxe = xmpielrank(n)

    call quadraturequad3d(n, igqrules, vext, qpoints2d, wequa2d)

    weights_q(1:qp_quad) = wequa2d(1:qp_quad)

    call quadraturetriang(n, igqrules, vext, qpoints2d, wequa2d)
    weights_t(1:qp_triangle) = wequa2d(1:qp_triangle)

    if (reduce_comp .eq. 1) then
      weights_t = 1.0d0; weights_q = 1.0d0
    end if

    !$omp do
    do ii = 1, nof_interior        !for all the interior elements
      i = el_int(ii)
      iconsidered = i

      if (dg .eq. 1) then
        dg_rhs = zero
        dg_rhs_surf_integ = zero
        dg_rhs_vol_integ = zero
      end if

      damp = lamx
      if (br2_yn .eq. 2) damp = 0.0d0
      if (mrf .eq. 1) then
        srf = ilocal_recon3(i)%mrf
      end if
      do l = 1, ielem(n, i)%ifca !for all their faces

        godflux2 = zero
        angle1 = ielem(n, i)%faceanglex(l)
        angle2 = ielem(n, i)%faceangley(l)
        nx = (cos(angle1)*sin(angle2))
        ny = (sin(angle1)*sin(angle2))
        nz = (cos(angle2))
        nnn(1) = nx; nnn(2) = ny; nnn(3) = nz
        if (ielem(n, i)%types_faces(l) .eq. 5) then
          iqp = qp_quad_n
          weights_temp(1:iqp) = weights_q(1:iqp)
        else
          iqp = qp_triangle_n
          weights_temp(1:iqp) = weights_t(1:iqp)
        end if

        if (dg .eq. 1) weights_dg = weights_temp

        do ngp = 1, iqp        !for all the gaussian quadrature points

          facex = l
          pointx = ngp

        call calculate_interior_viscous(n,b_code,iconsidered,facex,pointx,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speedrot,cleft,cright,lcvgrad,rcvgrad,lcvgrad_t,rcvgrad_t)

          if ((lmach .eq. 1)) then    !application of the low mach number correction
            call rotatef(n, cright_rot, cright, angle1, angle2)        !rotate wrt to normalvector of face and solve 1d riemann problem
            call rotatef(n, cleft_rot, cleft, angle1, angle2)        !rotate wrt to normalvector of face and solve 1d riemann problem
            leftv(1:nof_variables) = cleft_rot(1:nof_variables); rightv(1:nof_variables) = cright_rot(1:nof_variables)
            call lmacht(n, leftv, rightv)
            cleft_rot(1:nof_variables) = leftv(1:nof_variables); cright_rot(1:nof_variables) = rightv(1:nof_variables);
            call rotateb(n, cleft, cleft_rot, angle1, angle2)
            call rotateb(n, cright, cright_rot, angle1, angle2)
          end if

          leftv(1:nof_variables) = cleft(1:nof_variables); rightv(1:nof_variables) = cright(1:nof_variables)
          call cons2prim2(n, leftv, rightv, mp_pinfl, mp_pinfr, gammal, gammar)
          call sutherland(n, leftv, rightv, viscl, laml)

          if (turbulence .eq. 1) then
            if (turbulencemodel .eq. 1) then
              turbmv(1) = cturbl(1); turbmv(2) = cturbr(1); eddyfl(2) = turbmv(1); eddyfr(2) = turbmv(2)
              call eddyvisco(n, viscl, laml, turbmv, etvm, eddyfl, eddyfr, leftv, rightv)
            end if
            if (turbulencemodel .eq. 2) then
              eddyfl(1) = ielem(n, i)%walldist; eddyfl(2) = cturbl(1); eddyfl(3) = cturbl(2)
              eddyfl(4:6) = lcvgrad(1, 1:3); eddyfl(7:9) = lcvgrad(2, 1:3)
              eddyfl(10:12) = lcvgrad(3, 1:3); eddyfl(13:15) = lcvgrad_t(1, 1:3)
              eddyfl(16:18) = lcvgrad_t(2, 1:3)

              eddyfr(1) = ielem(n, i)%walldist; eddyfr(2) = cturbr(1); eddyfr(3) = cturbr(2)
              eddyfr(4:6) = rcvgrad(1, 1:3); eddyfr(7:9) = rcvgrad(2, 1:3); eddyfr(10:12) = rcvgrad(3, 1:3)
              eddyfr(13:15) = rcvgrad_t(1, 1:3); eddyfl(16:18) = rcvgrad_t(2, 1:3)
              call eddyvisco(n, viscl, laml, turbmv, etvm, eddyfl, eddyfr, leftv, rightv)
            end if
          end if

          taul = zero; tau = zero; taur = zero; q = zero; ux = zero; uy = zero; uz = zero; vx = zero; vy = zero; vz = zero; wx = zero; wy = zero; wz = zero;
          fxv = zero; fyv = zero; fzv = zero; rho12 = zero;
          u12 = zero; v12 = zero; w12 = zero

          vdamp = (4.0/3.0)!*(( (viscl(1))+(viscl(2)))))
          nall(1) = nx; nall(2) = ny; nall(3) = nz
          lcvgrad(1,1:3)=((lcvgrad(1,1:3)+rcvgrad(1,1:3))/(2.0d0))+damp*((vdamp/abs(ielem(n,i)%dih(l)))*nall(1:3)*(rightv(2)-leftv(2)))
          lcvgrad(2,1:3)=((lcvgrad(2,1:3)+rcvgrad(2,1:3))/(2.0d0))+damp*((vdamp/abs(ielem(n,i)%dih(l)))*nall(1:3)*(rightv(3)-leftv(3)))
          lcvgrad(3,1:3)=((lcvgrad(3,1:3)+rcvgrad(3,1:3))/(2.0d0))+damp*((vdamp/abs(ielem(n,i)%dih(l)))*nall(1:3)*(rightv(4)-leftv(4)))
          lcvgrad(4,1:3)=((lcvgrad(4,1:3)+rcvgrad(4,1:3))/(2.0d0))+damp*((vdamp/abs(ielem(n,i)%dih(l)))*nall(1:3)*((rightv(5)/rightv(1))-(leftv(5)/leftv(1))))

          if (turbulence .eq. 1) then
            q(1:3) = -oo2*((laml(3) + (laml(4)))*lcvgrad(4, 1:3))
          else
            q(1:3) = -oo2*((laml(1) + (laml(2)))*lcvgrad(4, 1:3))
          end if

          fxv(5) = fxv(5) - q(1); fyv(5) = fyv(5) - q(2); fzv(5) = fzv(5) - q(3)

          !left state derivatives
          ! determine taul!!
          ux = lcvgrad(1, 1); uy = lcvgrad(1, 2); uz = lcvgrad(1, 3);
          vx = lcvgrad(2, 1); vy = lcvgrad(2, 2); vz = lcvgrad(2, 3);
          wx = lcvgrad(3, 1); wy = lcvgrad(3, 2); wz = lcvgrad(3, 3);
          ! tau_xx
          taul(1, 1) = (4.0d0/3.0d0)*ux - (2.0d0/3.0d0)*vy - (2.0d0/3.0d0)*wz
          ! tau_yy
          taul(2, 2) = (4.0d0/3.0d0)*vy - (2.0d0/3.0d0)*ux - (2.0d0/3.0d0)*wz
          ! tau_zz
          taul(3, 3) = (4.0d0/3.0d0)*wz - (2.0d0/3.0d0)*ux - (2.0d0/3.0d0)*vy

          ! tau_xy
          taul(1, 2) = (uy + vx); taul(2, 1) = taul(1, 2)

          ! tau_xz
          taul(1, 3) = (wx + uz); taul(3, 1) = taul(1, 3)

          ! tau_yz
          taul(2, 3) = (vz + wy); taul(3, 2) = taul(2, 3)
          !end determine taul

          ! average and multiplay by viscosity
          if (turbulence .eq. 1) then
            tau = oo2*(((viscl(1) + viscl(3)) + (viscl(2) + viscl(4))))*taul
          else
            tau = oo2*(((viscl(1)) + (viscl(2))))*taul
          end if

          ! now addition into momentum fluxes
          do kc = 2, 4
            fxv(kc) = fxv(kc) + tau(1, kc - 1)
            fyv(kc) = fyv(kc) + tau(2, kc - 1)
            fzv(kc) = fzv(kc) + tau(3, kc - 1)
          end do

          ! compute interface velocities
          rho12 = oo2*(cleft(1) + cright(1))
          u12 = oo2*(cleft(2) + cright(2))/rho12
          v12 = oo2*(cleft(3) + cright(3))/rho12
          w12 = oo2*(cleft(4) + cright(4))/rho12

          fxv(5) = fxv(5) + u12*tau(1, 1) + v12*tau(1, 2) + w12*tau(1, 3)
          fyv(5) = fyv(5) + u12*tau(2, 1) + v12*tau(2, 2) + w12*tau(2, 3)
          fzv(5) = fzv(5) + u12*tau(3, 1) + v12*tau(3, 2) + w12*tau(3, 3)

          hllcflux(1:nof_variables) = (nx*fxv + ny*fyv + nz*fzv)

          if (dg .eq. 1) then
            rhllcflux(1:nof_variables) = hllcflux(1:nof_variables)
            dg_rhs_surf_integ = dg_rhs_surf_integ + dg_surf_flux(n, iconsidered, facex, pointx, weights_temp, rhllcflux)
          else

         godflux2(1:nof_variables) = godflux2(1:nof_variables) + (hllcflux(1:nof_variables)*(weights_temp(ngp)*ielem(n, i)%surf(l)))

          end if

          if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
            if (turbulence .eq. 1) then
              hllcflux(nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = &
              ((oo2*(viscl(1) + viscl(2))) + (oo2*(viscl(3) + viscl(4))))* &
              (((lcvgrad_t(1:turbulenceequations + passivescalar, 1) + rcvgrad_t(1:turbulenceequations + passivescalar, 1))*oo2*nx) + &
              ((lcvgrad_t(1:turbulenceequations + passivescalar, 2) + rcvgrad_t(1:turbulenceequations + passivescalar, 2))*oo2*ny) + &
              ((lcvgrad_t(1:turbulenceequations + passivescalar, 3) + rcvgrad_t(1:turbulenceequations + passivescalar, 3))*oo2*nz))
            else
              hllcflux(nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = &
             ((oo2*(viscl(1) + viscl(2))))* &
             (((lcvgrad_t(1:turbulenceequations + passivescalar, 1) + rcvgrad_t(1:turbulenceequations + passivescalar, 1))*oo2*nx) + &
             ((lcvgrad_t(1:turbulenceequations + passivescalar, 2) + rcvgrad_t(1:turbulenceequations + passivescalar, 2))*oo2*ny) + &
             ((lcvgrad_t(1:turbulenceequations + passivescalar, 3) + rcvgrad_t(1:turbulenceequations + passivescalar, 3))*oo2*nz))

            end if
            if (turbulencemodel .eq. 1) then
              hllcflux(6) = hllcflux(6)/sigma
            end if
              godflux2(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=godflux2(nof_variables+1:nof_variables+turbulenceequations+passivescalar)+&
             (hllcflux(nof_variables + 1:nof_variables + turbulenceequations + passivescalar)*(weights_temp(ngp)*ielem(n, i)%surf(l)))
          end if
        end do

        rhs(i)%val(1:nof_variables) = rhs(i)%val(1:nof_variables) - godflux2(1:nof_variables)

        if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
          rhst(i)%val(1:turbulenceequations + passivescalar) = rhst(i)%val(1:turbulenceequations + passivescalar) - &
          godflux2(nof_variables + 1:nof_variables + turbulenceequations + passivescalar)
        end if
      end do

      if (dg .eq. 1) then
        dg_rhs = dg_rhs_surf_integ
        rhs(i)%valdg = rhs(i)%valdg - dg_rhs
      end if

    end do
    !$omp end do

    !$omp do
    do ii = 1, nof_bounded
      i = el_bnd(ii)
      iconsidered = i
      if (dg .eq. 1) then
        dg_rhs = zero
        dg_rhs_surf_integ = zero
        dg_rhs_vol_integ = zero
      end if

      if (mrf .eq. 1) then
        srf = ilocal_recon3(i)%mrf
      end if

      do l = 1, ielem(n, i)%ifca

        igoflux = 0

        damp = lamx
        if (br2_yn .eq. 2) damp = 0.0d0
        angle1 = ielem(n, i)%faceanglex(l)
        angle2 = ielem(n, i)%faceangley(l)
        nx = (cos(angle1)*sin(angle2))
        ny = (sin(angle1)*sin(angle2))
        nz = (cos(angle2))
        nnn(1) = nx; nnn(2) = ny; nnn(3) = nz
        if (ielem(n, i)%types_faces(l) .eq. 5) then
          iqp = qp_quad_n
          weights_temp(1:iqp) = weights_q(1:iqp)
        else
          iqp = qp_triangle_n
          weights_temp(1:iqp) = weights_t(1:iqp)
        end if

        if (dg .eq. 1) weights_dg = weights_temp

        godflux2 = zero

        b_code = 0

        do ngp = 1, iqp
          facex = l
          pointx = ngp

        call calculate_bounded_viscous(n,b_code,iconsidered,facex,pointx,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speedrot,cleft,cright,lcvgrad,rcvgrad,lcvgrad_t,rcvgrad_t)

          if ((lmach .eq. 1)) then    !application of the low mach number correction
            call rotatef(n, cright_rot, cright, angle1, angle2)        !rotate wrt to normalvector of face and solve 1d riemann problem
            call rotatef(n, cleft_rot, cleft, angle1, angle2)        !rotate wrt to normalvector of face and solve 1d riemann problem
            leftv(1:nof_variables) = cleft_rot(1:nof_variables); rightv(1:nof_variables) = cright_rot(1:nof_variables)
            call lmacht(n, leftv, rightv)
            cleft_rot(1:nof_variables) = leftv(1:nof_variables); cright_rot(1:nof_variables) = rightv(1:nof_variables);
            call rotateb(n, cleft, cleft_rot, angle1, angle2)
            call rotateb(n, cright, cright_rot, angle1, angle2)
          end if

          leftv(1:nof_variables) = cleft(1:nof_variables); rightv(1:nof_variables) = cright(1:nof_variables)
          call cons2prim2(n, leftv, rightv, mp_pinfl, mp_pinfr, gammal, gammar)
          call sutherland(n, leftv, rightv, viscl, laml)

          if (turbulence .eq. 1) then
            if (turbulencemodel .eq. 1) then
              turbmv(1) = cturbl(1); turbmv(2) = cturbr(1); eddyfl(2) = turbmv(1); eddyfr(2) = turbmv(2)
              call eddyvisco(n, viscl, laml, turbmv, etvm, eddyfl, eddyfr, leftv, rightv)
            end if
            if (turbulencemodel .eq. 2) then
              eddyfl(1) = ielem(n, i)%walldist; eddyfl(2) = cturbl(1); eddyfl(3) = cturbl(2)
              eddyfl(4:6) = lcvgrad(1, 1:3); eddyfl(7:9) = lcvgrad(2, 1:3)
              eddyfl(10:12) = lcvgrad(3, 1:3); eddyfl(13:15) = lcvgrad_t(1, 1:3)
              eddyfl(16:18) = lcvgrad_t(2, 1:3)

              eddyfr(1) = ielem(n, i)%walldist; eddyfr(2) = cturbr(1); eddyfr(3) = cturbr(2)
              eddyfr(4:6) = rcvgrad(1, 1:3); eddyfr(7:9) = rcvgrad(2, 1:3); eddyfr(10:12) = rcvgrad(3, 1:3)
              eddyfr(13:15) = rcvgrad_t(1, 1:3); eddyfl(16:18) = rcvgrad_t(2, 1:3)
              call eddyvisco(n, viscl, laml, turbmv, etvm, eddyfl, eddyfr, leftv, rightv)
            end if
          end if

          taul = zero; tau = zero; taur = zero; q = zero; ux = zero; uy = zero; uz = zero; vx = zero; vy = zero; vz = zero; wx = zero; wy = zero; wz = zero;
          fxv = zero; fyv = zero; fzv = zero; rho12 = zero;
          u12 = zero; v12 = zero; w12 = zero

          if ((b_code .lt. 5) .and. (b_code .gt. 0)) then
            damp = zero
          end if

          vdamp = (4.0/3.0)!*(( (viscl(1))+(viscl(2)))))
          nall(1) = nx; nall(2) = ny; nall(3) = nz
          lcvgrad(1,1:3)=((lcvgrad(1,1:3)+rcvgrad(1,1:3))/(2.0d0))+damp*((vdamp/abs(ielem(n,i)%dih(l)))*nall(1:3)*(rightv(2)-leftv(2)))
          lcvgrad(2,1:3)=((lcvgrad(2,1:3)+rcvgrad(2,1:3))/(2.0d0))+damp*((vdamp/abs(ielem(n,i)%dih(l)))*nall(1:3)*(rightv(3)-leftv(3)))
          lcvgrad(3,1:3)=((lcvgrad(3,1:3)+rcvgrad(3,1:3))/(2.0d0))+damp*((vdamp/abs(ielem(n,i)%dih(l)))*nall(1:3)*(rightv(4)-leftv(4)))
          lcvgrad(4,1:3)=((lcvgrad(4,1:3)+rcvgrad(4,1:3))/(2.0d0))+damp*((vdamp/abs(ielem(n,i)%dih(l)))*nall(1:3)*((rightv(5)/rightv(1))-(leftv(5)/leftv(1))))

          if (turbulence .eq. 1) then
            q(1:3) = -oo2*((laml(3) + (laml(4)))*lcvgrad(4, 1:3))
          else
            q(1:3) = -oo2*((laml(1) + (laml(2)))*lcvgrad(4, 1:3))
          end if

          fxv(5) = fxv(5) - q(1); fyv(5) = fyv(5) - q(2); fzv(5) = fzv(5) - q(3)

          !left state derivatives
          ! determine taul!!
          ux = lcvgrad(1, 1); uy = lcvgrad(1, 2); uz = lcvgrad(1, 3);
          vx = lcvgrad(2, 1); vy = lcvgrad(2, 2); vz = lcvgrad(2, 3);
          wx = lcvgrad(3, 1); wy = lcvgrad(3, 2); wz = lcvgrad(3, 3);
          ! tau_xx
          taul(1, 1) = (4.0d0/3.0d0)*ux - (2.0d0/3.0d0)*vy - (2.0d0/3.0d0)*wz
          ! tau_yy
          taul(2, 2) = (4.0d0/3.0d0)*vy - (2.0d0/3.0d0)*ux - (2.0d0/3.0d0)*wz
          ! tau_zz
          taul(3, 3) = (4.0d0/3.0d0)*wz - (2.0d0/3.0d0)*ux - (2.0d0/3.0d0)*vy

          ! tau_xy
          taul(1, 2) = (uy + vx); taul(2, 1) = taul(1, 2)

          ! tau_xz
          taul(1, 3) = (wx + uz); taul(3, 1) = taul(1, 3)

          ! tau_yz
          taul(2, 3) = (vz + wy); taul(3, 2) = taul(2, 3)
          !end determine taul

          ! average and multiplay by viscosity
          if (turbulence .eq. 1) then
            tau = oo2*(((viscl(1) + viscl(3)) + (viscl(2) + viscl(4))))*taul
          else
            tau = oo2*(((viscl(1)) + (viscl(2))))*taul
          end if

          ! now addition into momentum fluxes
          do kc = 2, 4
            fxv(kc) = fxv(kc) + tau(1, kc - 1)
            fyv(kc) = fyv(kc) + tau(2, kc - 1)
            fzv(kc) = fzv(kc) + tau(3, kc - 1)
          end do

          ! compute interface velocities
          rho12 = oo2*(cleft(1) + cright(1))
          u12 = oo2*(cleft(2) + cright(2))/rho12
          v12 = oo2*(cleft(3) + cright(3))/rho12
          w12 = oo2*(cleft(4) + cright(4))/rho12

          fxv(5) = fxv(5) + u12*tau(1, 1) + v12*tau(1, 2) + w12*tau(1, 3)
          fyv(5) = fyv(5) + u12*tau(2, 1) + v12*tau(2, 2) + w12*tau(2, 3)
          fzv(5) = fzv(5) + u12*tau(3, 1) + v12*tau(3, 2) + w12*tau(3, 3)

          hllcflux(1:nof_variables) = (nx*fxv + ny*fyv + nz*fzv)

          if (dg .eq. 1) then
            rhllcflux(1:nof_variables) = hllcflux(1:nof_variables)
            dg_rhs_surf_integ = dg_rhs_surf_integ + dg_surf_flux(n, iconsidered, facex, pointx, weights_temp, rhllcflux)
          else
            godflux2(1:nof_variables) = godflux2(1:nof_variables) + (hllcflux(1:nof_variables)*(weights_temp(ngp)*ielem(n, i)%surf(l)))
          end if

          if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
            if (turbulence .eq. 1) then
              hllcflux(nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = &
              ((oo2*(viscl(1) + viscl(2))) + (oo2*(viscl(3) + viscl(4))))* &
              (((lcvgrad_t(1:turbulenceequations + passivescalar, 1) + rcvgrad_t(1:turbulenceequations + passivescalar, 1))*oo2*nx) + &
              ((lcvgrad_t(1:turbulenceequations + passivescalar, 2) + rcvgrad_t(1:turbulenceequations + passivescalar, 2))*oo2*ny) + &
              ((lcvgrad_t(1:turbulenceequations + passivescalar, 3) + rcvgrad_t(1:turbulenceequations + passivescalar, 3))*oo2*nz))
            else
              hllcflux(nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = &
              ((oo2*(viscl(1) + viscl(2))))* &
              (((lcvgrad_t(1:turbulenceequations + passivescalar, 1) + rcvgrad_t(1:turbulenceequations + passivescalar, 1))*oo2*nx) + &
              ((lcvgrad_t(1:turbulenceequations + passivescalar, 2) + rcvgrad_t(1:turbulenceequations + passivescalar, 2))*oo2*ny) + &
              ((lcvgrad_t(1:turbulenceequations + passivescalar, 3) + rcvgrad_t(1:turbulenceequations + passivescalar, 3))*oo2*nz))

            end if
            if (turbulencemodel .eq. 1) then
              hllcflux(6) = hllcflux(6)/sigma
            end if

            if ((b_code .eq. 3)) then
              hllcflux(nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = zero
            end if
              godflux2(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=godflux2(nof_variables+1:nof_variables+turbulenceequations+passivescalar)+&
              (hllcflux(nof_variables + 1:nof_variables + turbulenceequations + passivescalar)*(weights_temp(ngp)*ielem(n, i)%surf(l)))
          end if
        end do

        rhs(i)%val(1:nof_variables) = rhs(i)%val(1:nof_variables) - godflux2(1:nof_variables)

        if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
          rhst(i)%val(1:turbulenceequations + passivescalar) = rhst(i)%val(1:turbulenceequations + passivescalar) - &
          godflux2(nof_variables + 1:nof_variables + turbulenceequations + passivescalar)
        end if
      end do

      if (dg .eq. 1) then
        dg_rhs = dg_rhs_surf_integ
        rhs(i)%valdg = rhs(i)%valdg - dg_rhs
      end if

    end do
    !$omp end do

    if (dg .eq. 1) then
      deallocate(dg_rhs, dg_rhs_vol_integ, dg_rhs_surf_integ)
    end if

  end subroutine calculate_fluxeshi_diffusive

  subroutine calculate_fluxeshi_diffusive2d(n)
!> @brief
!> this subroutine computes the diffusive fluxes for euler-navier-stokes equations in 2d
    implicit none
    integer, intent(in)::n
    real, dimension(1:nof_variables + turbulenceequations + passivescalar)::godflux2, dg_vol_rec, rhllcflux, hllcflux
    integer::i, l, ngp, kmaxe, iqp, ii, nvar, kc, iex, ittt, ikas, igoflux, icaseb, kk, b_code
    real::sum_detect, norms
    integer::iconsidered, facex, pointx
    real::angle1, angle2, nx, ny, nz, mp_source1, mp_source2, mp_source3
    real, dimension(1:nof_variables + turbulenceequations + passivescalar)::cleft, cright, cleft_rot, cright_rot
    real, dimension(1:nof_variables)::leftv, rightv, srf_speedrot
    real, dimension(1:turbulenceequations + passivescalar)::cturbl, cturbr
    real, dimension(1:dimensiona)::pox, poy, poz
    real, dimension(1:nof_variables)::srf_speed
    real, dimension(1:numberofpoints2)::weights_q, weights_t, weights_temp, weights_dg
    real, dimension(1:8, 1:dimensiona)::vext
    real, dimension(1:dimensiona, 1:numberofpoints2)::qpoints2d
    real, dimension(1:numberofpoints2)::wequa2d
    real::mp_pinfl, mp_pinfr, gammal, gammar
    real, dimension(1:4)::viscl, laml
    real, dimension(1:2)::turbmv
    real, dimension(1)::etvm
    real, dimension(1:20)::eddyfl, eddyfr
    real, dimension(1:nof_variables - 1, 1:dims)::lcvgrad, rcvgrad
    real, dimension(turbulenceequations + passivescalar, 1:dims)::lcvgrad_t, rcvgrad_t
    real, dimension(4)::fxv, fyv, fzv, tem_pn, rtem_pn
    real, dimension(2, 2)::taul, taur, tau
    real, dimension(2)::q, nall
    real::ux, uy, uz, vx, vy, vz, wx, wy, wz, rho12, u12, v12, w12, damp, vdamp
    real, allocatable, dimension(:, :)::dg_rhs, dg_rhs_vol_integ, dg_rhs_surf_integ

    if (dg .eq. 1) then
        allocate(dg_rhs(1:num_dg_dofs,1:nof_variables), dg_rhs_vol_integ(1:num_dg_dofs,1:nof_variables), dg_rhs_surf_integ(1:num_dg_dofs,1:nof_variables))
    end if

    kmaxe = xmpielrank(n)

    call quadratureline(n, igqrules, vext, qpoints2d, wequa2d)

    weights_temp(1:qp_line_n) = wequa2d(1:qp_line_n)
    if (reduce_comp .eq. 1) then
      weights_temp = 1.0d0
    end if

    !$omp do
    do ii = 1, nof_interior        !for all the interior elements
      i = el_int(ii)
      iconsidered = i

      if (dg .eq. 1) then
        dg_rhs = zero
        dg_rhs_surf_integ = zero
        dg_rhs_vol_integ = zero

      end if

      do l = 1, ielem(n, i)%ifca !for all their faces
!                 if (ielem(n,i)%reorient(l).eq.0)then
        damp = lamx
        if (br2_yn .eq. 2) damp = 0.0d0

        godflux2 = zero
        angle1 = ielem(n, i)%faceanglex(l)
        angle2 = ielem(n, i)%faceangley(l)
        nx = angle1
        ny = angle2

        iqp = qp_line_n

        do ngp = 1, iqp        !for all the gaussian quadrature points
          facex = l
          pointx = ngp

          call calculate_interior_viscous2d(n,b_code,iconsidered,facex,pointx,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speedrot,cleft,cright,lcvgrad,rcvgrad,lcvgrad_t,rcvgrad_t)

          if ((lmach .eq. 1)) then    !application of the low mach number correction
            call rotatef2d(n, cright_rot, cright, angle1, angle2)        !rotate wrt to normalvector of face and solve 1d riemann problem
            call rotatef2d(n, cleft_rot, cleft, angle1, angle2)        !rotate wrt to normalvector of face and solve 1d riemann problem
            leftv(1:nof_variables) = cleft_rot(1:nof_variables); rightv(1:nof_variables) = cright_rot(1:nof_variables)
            call lmacht2d(n, leftv, rightv)
            cleft_rot(1:nof_variables) = leftv(1:nof_variables); cright_rot(1:nof_variables) = rightv(1:nof_variables);
            call rotateb2d(n, cleft, cleft_rot, angle1, angle2)
            call rotateb2d(n, cright, cright_rot, angle1, angle2)
          end if

          leftv(1:nof_variables) = cleft(1:nof_variables); rightv(1:nof_variables) = cright(1:nof_variables)
          call cons2prim2(n, leftv, rightv, mp_pinfl, mp_pinfr, gammal, gammar)
          call sutherland2d(n, leftv, rightv, viscl, laml)

          if (turbulence .eq. 1) then
            if (turbulencemodel .eq. 1) then
              turbmv(1) = cturbl(1); turbmv(2) = cturbr(1); eddyfl(2) = turbmv(1); eddyfr(2) = turbmv(2)
              call eddyvisco2d(n, viscl, laml, turbmv, etvm, eddyfl, eddyfr, leftv, rightv)
            end if
            if (turbulencemodel .eq. 2) then
              eddyfl(1) = ielem(n, i)%walldist; eddyfl(2) = cturbl(1); eddyfl(3) = cturbl(2)
              eddyfl(4:5) = lcvgrad(1, 1:2); eddyfl(6:7) = lcvgrad(2, 1:2)
              eddyfl(8:9) = lcvgrad_t(1, 1:2)
              eddyfl(10:11) = lcvgrad_t(2, 1:2)
              eddyfr(1) = ielem(n, i)%walldist; eddyfr(2) = cturbr(1); eddyfr(3) = cturbr(2)
              eddyfr(4:5) = rcvgrad(1, 1:2); eddyfr(6:7) = rcvgrad(2, 1:2)
              eddyfr(8:9) = rcvgrad_t(1, 1:2); eddyfl(10:11) = rcvgrad_t(2, 1:2)
              call eddyvisco2d(n, viscl, laml, turbmv, etvm, eddyfl, eddyfr, leftv, rightv)
            end if
          end if

          taul = zero; tau = zero; taur = zero; q = zero; ux = zero; uy = zero; uz = zero; vx = zero; vy = zero; vz = zero; wx = zero; wy = zero; wz = zero;
          fxv = zero; fyv = zero; fzv = zero; rho12 = zero;
          u12 = zero; v12 = zero; w12 = zero

          vdamp = (4.0/3.0)!*(( (viscl(1))+(viscl(2)))))
          nall(1) = nx; nall(2) = ny

          lcvgrad(1,1:2)=((lcvgrad(1,1:2)+rcvgrad(1,1:2))/(2.0d0))+damp*((vdamp/abs(ielem(n,i)%dih(l)))*nall(1:2)*(rightv(2)-leftv(2)))
          lcvgrad(2,1:2)=((lcvgrad(2,1:2)+rcvgrad(2,1:2))/(2.0d0))+damp*((vdamp/abs(ielem(n,i)%dih(l)))*nall(1:2)*(rightv(3)-leftv(3)))
          lcvgrad(3,1:2)=((lcvgrad(3,1:2)+rcvgrad(3,1:2))/(2.0d0))+damp*((vdamp/abs(ielem(n,i)%dih(l)))*nall(1:2)*((rightv(4)/rightv(1))-(leftv(4)/leftv(1))))

          if (turbulence .eq. 1) then
            q(1:2) = -oo2*((laml(3) + (laml(4)))*lcvgrad(3, 1:2))
          else
            q(1:2) = -oo2*((laml(1) + (laml(2)))*lcvgrad(3, 1:2))
          end if

          fxv(4) = fxv(4) - q(1); fyv(4) = fyv(4) - q(2)

          !left state derivatives
          ux = lcvgrad(1, 1); uy = lcvgrad(1, 2)
          vx = lcvgrad(2, 1); vy = lcvgrad(2, 2)
          ! determine taul!!

          ! tau_xx
          taul(1, 1) = (4.0d0/3.0d0)*ux - (2.0d0/3.0d0)*vy
          ! tau_yy
          taul(2, 2) = (4.0d0/3.0d0)*vy - (2.0d0/3.0d0)*ux
          ! tau_zz
          ! tau_xy
          taul(1, 2) = (uy + vx); taul(2, 1) = taul(1, 2)

          ! average and multiplay by viscosity
          if (turbulence .eq. 1) then
            tau = oo2*(((viscl(1) + viscl(3))) + ((viscl(2) + viscl(4))))*taul
          else
            tau = oo2*((viscl(1)) + (viscl(2)))*taul
          end if

          ! now addition into momentum fluxes
          do kc = 2, 3
            fxv(kc) = fxv(kc) + tau(1, kc - 1)
            fyv(kc) = fyv(kc) + tau(2, kc - 1)

          end do

          ! compute interface velocities
          rho12 = oo2*(cleft(1) + cright(1))
          u12 = oo2*(cleft(2) + cright(2))/rho12
          v12 = oo2*(cleft(3) + cright(3))/rho12
          fxv(4) = fxv(4) + u12*tau(1, 1) + v12*tau(1, 2)
          fyv(4) = fyv(4) + u12*tau(2, 1) + v12*tau(2, 2)
          hllcflux(1:nof_variables) = (nx*fxv + ny*fyv)

          if (dg .eq. 1) then

            rhllcflux(1:nof_variables) = hllcflux(1:nof_variables)
            dg_rhs_surf_integ = dg_rhs_surf_integ + dg_surf_flux(n, iconsidered, facex, pointx, weights_temp, rhllcflux)
          else

         godflux2(1:nof_variables) = godflux2(1:nof_variables) + (hllcflux(1:nof_variables)*(weights_temp(ngp)*ielem(n, i)%surf(l)))
          end if

          if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
            if (turbulence .eq. 1) then
              hllcflux(nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = &
              ((oo2*(viscl(1) + viscl(2))) + (oo2*(viscl(3) + viscl(4))))* &
              (((lcvgrad_t(1:turbulenceequations + passivescalar, 1) + rcvgrad_t(1:turbulenceequations + passivescalar, 1))*oo2*nx) + &
              ((lcvgrad_t(1:turbulenceequations + passivescalar, 2) + rcvgrad_t(1:turbulenceequations + passivescalar, 2))*oo2*ny))
            else
              hllcflux(nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = &
              ((oo2*(viscl(1) + viscl(2))))* &
              (((lcvgrad_t(1:turbulenceequations + passivescalar, 1) + rcvgrad_t(1:turbulenceequations + passivescalar, 1))*oo2*nx) + &
              ((lcvgrad_t(1:turbulenceequations + passivescalar, 2) + rcvgrad_t(1:turbulenceequations + passivescalar, 2))*oo2*ny))

            end if
            if (turbulencemodel .eq. 1) then
              hllcflux(5) = hllcflux(5)/sigma
            end if
              godflux2(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=godflux2(nof_variables+1:nof_variables+turbulenceequations+passivescalar)+&
             (hllcflux(nof_variables + 1:nof_variables + turbulenceequations + passivescalar)*(weights_temp(ngp)*ielem(n, i)%surf(l)))
          end if
        end do

        rhs(i)%val(1:nof_variables) = rhs(i)%val(1:nof_variables) - godflux2(1:nof_variables)

        if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
          rhst(i)%val(1:turbulenceequations + passivescalar) = rhst(i)%val(1:turbulenceequations + passivescalar) - &
          godflux2(nof_variables + 1:nof_variables + turbulenceequations + passivescalar)
        end if
      end do

      if (dg .eq. 1) then
        dg_rhs = dg_rhs_surf_integ
        rhs(i)%valdg = rhs(i)%valdg - dg_rhs
      end if

    end do
    !$omp end do
    !$omp do
    do ii = 1, nof_bounded
      i = el_bnd(ii)
      iconsidered = i
      if (dg .eq. 1) then
        dg_rhs = zero
        dg_rhs_surf_integ = zero
        dg_rhs_vol_integ = zero
      end if

      do l = 1, ielem(n, i)%ifca
        damp = lamx
        if (br2_yn .eq. 2) damp = 0.0d0
        b_code = 0
        godflux2 = zero
        angle1 = ielem(n, i)%faceanglex(l)
        angle2 = ielem(n, i)%faceangley(l)
        nx = angle1
        ny = angle2
        iqp = qp_line_n
        do ngp = 1, iqp
          facex = l
          pointx = ngp
        call calculate_bounded_viscous2d(n,b_code,iconsidered,facex,pointx,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speedrot,cleft,cright,lcvgrad,rcvgrad,lcvgrad_t,rcvgrad_t)
          if ((lmach .eq. 1)) then    !application of the low mach number correction
            call rotatef2d(n, cright_rot, cright, angle1, angle2)        !rotate wrt to normalvector of face and solve 1d riemann problem
            call rotatef2d(n, cleft_rot, cleft, angle1, angle2)        !rotate wrt to normalvector of face and solve 1d riemann problem
            leftv(1:nof_variables) = cleft_rot(1:nof_variables); rightv(1:nof_variables) = cright_rot(1:nof_variables)
            call lmacht2d(n, leftv, rightv)
            cleft_rot(1:nof_variables) = leftv(1:nof_variables); cright_rot(1:nof_variables) = rightv(1:nof_variables);
            call rotateb2d(n, cleft, cleft_rot, angle1, angle2)
            call rotateb2d(n, cright, cright_rot, angle1, angle2)
          end if

          leftv(1:nof_variables) = cleft(1:nof_variables); rightv(1:nof_variables) = cright(1:nof_variables)
          call cons2prim2(n, leftv, rightv, mp_pinfl, mp_pinfr, gammal, gammar)
          call sutherland2d(n, leftv, rightv, viscl, laml)

          if (turbulence .eq. 1) then
            if (turbulencemodel .eq. 1) then
              turbmv(1) = cturbl(1); turbmv(2) = cturbr(1); eddyfl(2) = turbmv(1); eddyfr(2) = turbmv(2)
              call eddyvisco2d(n, viscl, laml, turbmv, etvm, eddyfl, eddyfr, leftv, rightv)
              if (b_code .eq. 4) then
                viscl(3:4) = zero
                laml(3:4) = zero
              end if
            end if
            if (turbulencemodel .eq. 2) then
              eddyfl(1) = ielem(n, i)%walldist; eddyfl(2) = cturbl(1); eddyfl(3) = cturbl(2)
              eddyfl(4:5) = lcvgrad(1, 1:2); eddyfl(6:7) = lcvgrad(2, 1:2)
              eddyfl(8:9) = lcvgrad_t(1, 1:2)
              eddyfl(10:11) = lcvgrad_t(2, 1:2)

              eddyfr(1) = ielem(n, i)%walldist; eddyfr(2) = cturbr(1); eddyfr(3) = cturbr(2)
              eddyfr(4:5) = rcvgrad(1, 1:2); eddyfr(6:7) = rcvgrad(2, 1:2)
              eddyfr(8:9) = rcvgrad_t(1, 1:2); eddyfl(10:11) = rcvgrad_t(2, 1:2)
              call eddyvisco2d(n, viscl, laml, turbmv, etvm, eddyfl, eddyfr, leftv, rightv)
            end if
          end if

          taul = zero; tau = zero; taur = zero; q = zero; ux = zero; uy = zero; uz = zero; vx = zero; vy = zero; vz = zero; wx = zero; wy = zero; wz = zero;
          fxv = zero; fyv = zero; fzv = zero; rho12 = zero;
          u12 = zero; v12 = zero; w12 = zero

          if ((b_code .lt. 5) .and. (b_code .gt. 0)) then
            damp = zero
          end if

          vdamp = (4.0/3.0)!*(( (viscl(1))+(viscl(2)))))
          nall(1) = nx; nall(2) = ny
          lcvgrad(1,1:2)=((lcvgrad(1,1:2)+rcvgrad(1,1:2))/(2.0d0))+damp*((vdamp/abs(ielem(n,i)%dih(l)))*nall(1:2)*(rightv(2)-leftv(2)))
          lcvgrad(2,1:2)=((lcvgrad(2,1:2)+rcvgrad(2,1:2))/(2.0d0))+damp*((vdamp/abs(ielem(n,i)%dih(l)))*nall(1:2)*(rightv(3)-leftv(3)))
          lcvgrad(3,1:2)=((lcvgrad(3,1:2)+rcvgrad(3,1:2))/(2.0d0))+damp*((vdamp/abs(ielem(n,i)%dih(l)))*nall(1:2)*((rightv(4)/rightv(1))-(leftv(4)/leftv(1))))

          if (turbulence .eq. 1) then
            q(1:2) = -oo2*((laml(3) + (laml(4)))*lcvgrad(3, 1:2))
          else
            q(1:2) = -oo2*((laml(1) + (laml(2)))*lcvgrad(3, 1:2))
          end if

          fxv(4) = fxv(4) - q(1); fyv(4) = fyv(4) - q(2)

          !left state derivatives
          ux = lcvgrad(1, 1); uy = lcvgrad(1, 2)
          vx = lcvgrad(2, 1); vy = lcvgrad(2, 2)
          ! determine taul!!

          ! tau_xx
          taul(1, 1) = (4.0d0/3.0d0)*ux - (2.0d0/3.0d0)*vy
          ! tau_yy
          taul(2, 2) = (4.0d0/3.0d0)*vy - (2.0d0/3.0d0)*ux
          ! tau_zz

          ! tau_xy
          taul(1, 2) = (uy + vx); taul(2, 1) = taul(1, 2)

          ! average and multiplay by viscosity
          if (turbulence .eq. 1) then
            tau = oo2*(((viscl(1) + viscl(3))) + ((viscl(2) + viscl(4))))*taul
          else
            tau = oo2*((viscl(1)) + (viscl(2)))*taul
          end if

          ! now addition into momentum fluxes
          do kc = 2, 3
            fxv(kc) = fxv(kc) + tau(1, kc - 1)
            fyv(kc) = fyv(kc) + tau(2, kc - 1)

          end do

          ! compute interface velocities
          rho12 = oo2*(cleft(1) + cright(1))
          u12 = oo2*(cleft(2) + cright(2))/rho12
          v12 = oo2*(cleft(3) + cright(3))/rho12

          fxv(4) = fxv(4) + u12*tau(1, 1) + v12*tau(1, 2)
          fyv(4) = fyv(4) + u12*tau(2, 1) + v12*tau(2, 2)

          hllcflux(1:nof_variables) = (nx*fxv + ny*fyv)

          if (dg .eq. 1) then
            rhllcflux(1:nof_variables) = hllcflux(1:nof_variables)
            dg_rhs_surf_integ = dg_rhs_surf_integ + dg_surf_flux(n, iconsidered, facex, pointx, weights_temp, rhllcflux)
          else
          godflux2(1:nof_variables) = godflux2(1:nof_variables) + (hllcflux(1:nof_variables)*(weights_temp(ngp)*ielem(n, i)%surf(l)))
          end if

          if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
            if (turbulence .eq. 1) then
              hllcflux(nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = &
               ((oo2*(viscl(1) + viscl(2))) + (oo2*(viscl(3) + viscl(4))))* &
               (((lcvgrad_t(1:turbulenceequations + passivescalar, 1) + rcvgrad_t(1:turbulenceequations + passivescalar, 1))*oo2*nx) + &
               ((lcvgrad_t(1:turbulenceequations + passivescalar, 2) + rcvgrad_t(1:turbulenceequations + passivescalar, 2))*oo2*ny))
            else
              hllcflux(nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = &
              ((oo2*(viscl(1) + viscl(2))))* &
              (((lcvgrad_t(1:turbulenceequations + passivescalar, 1) + rcvgrad_t(1:turbulenceequations + passivescalar, 1))*oo2*nx) + &
               ((lcvgrad_t(1:turbulenceequations + passivescalar, 2) + rcvgrad_t(1:turbulenceequations + passivescalar, 2))*oo2*ny))
            end if
            if (turbulencemodel .eq. 1) then
              hllcflux(5) = hllcflux(5)/sigma
            end if

            if ((b_code .eq. 3)) then
              hllcflux(nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = zero
            end if
             godflux2(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=godflux2(nof_variables+1:nof_variables+turbulenceequations+passivescalar)+&
             (hllcflux(nof_variables + 1:nof_variables + turbulenceequations + passivescalar)*(weights_temp(ngp)*ielem(n, i)%surf(l)))
          end if
        end do

        rhs(i)%val(1:nof_variables) = rhs(i)%val(1:nof_variables) - godflux2(1:nof_variables)

        if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
          rhst(i)%val(1:turbulenceequations + passivescalar) = rhst(i)%val(1:turbulenceequations + passivescalar) - &
          godflux2(nof_variables + 1:nof_variables + turbulenceequations + passivescalar)
        end if
      end do

      if (dg .eq. 1) then
        dg_rhs = dg_rhs_surf_integ
        rhs(i)%valdg = rhs(i)%valdg - dg_rhs
      end if

    end do
    !$omp end do
    if (dg .eq. 1) then
      deallocate(dg_rhs, dg_rhs_vol_integ, dg_rhs_surf_integ)
    end if

  end subroutine calculate_fluxeshi_diffusive2d

  subroutine calculate_fluxeshi_convective_mood(n)
!> @brief
!> this subroutine computes the convective fluxes for hyperbolic conservation laws
    implicit none
    integer, intent(in)::n
    real, dimension(1:nof_variables + turbulenceequations + passivescalar)::godflux2, dg_vol_rec, rhllcflux, hllcflux
    integer::i, l, ngp, kmaxe, iqp, ii, ikas, igoflux, icaseb, jx, jx2, b_code
    real::sum_detect, norms, tempxx
    real, dimension(1:numberofpoints2)::weights_q, weights_t, weights_temp, weights_dg
    integer::iconsidered, facex, pointx, n_node
    real::angle1, angle2, nx, ny, nz, mp_source1, mp_source2, mp_source3
    real, dimension(1:nof_variables + turbulenceequations + passivescalar)::cleft, cright, cleft_rot, cright_rot
    real, dimension(1:nof_variables)::leftv, rightv, srf_speedrot
    real, dimension(1:turbulenceequations + passivescalar)::cturbl, cturbr
    real, dimension(1:dimensiona)::pox, poy, poz
    real, dimension(1:nof_variables)::srf_speed
    real, dimension(1:8, 1:dimensiona)::vext
    real, dimension(1:dimensiona, 1:numberofpoints2)::qpoints2d
    real, dimension(1:numberofpoints2)::wequa2d
    real, dimension(1:8, 1:dimensiona)::nodes_list
    real, dimension(1:dimensiona)::cords
    integer::ibfc
    kmaxe = xmpielrank(n)

    kmaxe = xmpielrank(n)

    call quadraturequad3d(n, igqrules, vext, qpoints2d, wequa2d)

    weights_q(1:qp_quad) = wequa2d(1:qp_quad)

    call quadraturetriang(n, igqrules, vext, qpoints2d, wequa2d)
    weights_t(1:qp_triangle) = wequa2d(1:qp_triangle)

    if (reduce_comp .eq. 1) then
      weights_t = 1.0d0; weights_q = 1.0d0
    end if
    do i = 1, xmpielrank(n)
      if (ielem(n, i)%recalc .eq. 1) then
        rhs(i)%val(:) = zero; if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) rhst(i)%val(:) = zero
      end if
    end do
    !$omp barrier
    !$omp do
    do ii = 1, nof_interior        !for all the interior elements
      i = el_int(ii)
      iconsidered = i
      mp_source3 = zero
      b_code = 0
      if (ielem(n, i)%recalc .eq. 1) then
        do l = 1, ielem(n, i)%ifca !for all their faces
!                                      if (ielem(n,i)%reorient(l).eq.0)then
          godflux2 = zero
          mp_source2 = zero
          angle1 = ielem(n, i)%faceanglex(l)
          angle2 = ielem(n, i)%faceangley(l)
          nx = (cos(angle1)*sin(angle2))
          ny = (sin(angle1)*sin(angle2))
          nz = (cos(angle2))
          if (ielem(n, i)%types_faces(l) .eq. 5) then
            iqp = qp_quad_n
            weights_temp(1:iqp) = weights_q(1:iqp)
          else
            iqp = qp_triangle_n
            weights_temp(1:iqp) = weights_t(1:iqp)
          end if
          do ngp = 1, iqp        !for all the gaussian quadrature points
            cleft(1:nof_variables) = ilocal_recon3(i)%uleft(1:nof_variables, l, ngp)        !left mean flow state
            if ((cascade .eq. 2) .and. (ielem(n, i)%mood .eq. 1)) then
              cleft(1:nof_variables) = u_c(i)%val(3, 1:nof_variables)
            end if
            cright(1:nof_variables) = ilocal_recon3(ielem(n, i)%ineigh(l))%uleft(1:nof_variables, ielem(n, i)%ineighn(l), ngp)
            if ((cascade .eq. 2) .and. (ielem(n, (ielem(n, i)%ineigh(l)))%mood .eq. 1)) then
              cright(1:nof_variables) = u_c(ielem(n, i)%ineigh(l))%val(3, 1:nof_variables)
            end if
            !right mean flow state
            if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
              if (icoupleturb .eq. 1) then

cturbl(1:turbulenceequations + passivescalar) = ilocal_recon3(i)%uleftturb(1:turbulenceequations + passivescalar, l, ngp) !left additional equations flow state
            cturbr(1:turbulenceequations+passivescalar)=ilocal_recon3(ielem(n,i)%ineigh(l))%uleftturb(1:turbulenceequations+passivescalar,ielem(n,i)%ineighn(l),ngp)!right additional equations flow state
            else
              cturbl(1:turbulenceequations + passivescalar) = u_ct(i)%val(1, 1:turbulenceequations + passivescalar)
              cturbr(1:turbulenceequations + passivescalar) = u_ct(ielem(n, i)%ineigh(l))%val(1, 1:turbulenceequations + passivescalar)
            end if

            cleft_rot(nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = cturbl(1:turbulenceequations + passivescalar)
            cright_rot(nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = cturbr(1:turbulenceequations + passivescalar)
            end if

            call rotatef(n, cright_rot, cright, angle1, angle2)        !rotate wrt to normalvector of face and solve 1d riemann problem
            call rotatef(n, cleft_rot, cleft, angle1, angle2)        !rotate wrt to normalvector of face and solve 1d riemann problem

            if ((lmach .eq. 1)) then    !application of the low mach number correction
              leftv(1:nof_variables) = cleft_rot(1:nof_variables); rightv(1:nof_variables) = cright_rot(1:nof_variables)
              call lmacht(n, leftv, rightv)
              cleft_rot(1:nof_variables) = leftv(1:nof_variables); cright_rot(1:nof_variables) = rightv(1:nof_variables);
            end if

            select case (iriemann)

            case (1)                        !hllc
              call hllc_riemann_solver(n, iconsidered, facex, cleft_rot, cright_rot, hllcflux, mp_source1, srf_speedrot)
              call rotateb(n, rhllcflux, hllcflux, angle1, angle2)
              if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
                rhllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=hllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)
              end if

            case (2)                        !rusanov
              call rusanov_riemann_solver(n, iconsidered, facex, cleft_rot, cright_rot, hllcflux, mp_source1, srf_speedrot)
              call rotateb(n, rhllcflux, hllcflux, angle1, angle2)
              if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
                rhllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=hllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)
              end if

            case (9)                        !hll
              call hll_riemann_solver(n, cleft_rot, cright_rot, hllcflux, mp_source1, srf_speedrot)
              call rotateb(n, rhllcflux, hllcflux, angle1, angle2)
              if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
                rhllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=hllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)
              end if

            case (3)                        !roe
              call rotateb(n, cleft, cleft_rot, angle1, angle2)
              call rotateb(n, cright, cright_rot, angle1, angle2)
              call roe_riemann_solver(n, iconsidered, facex, cleft, cright, hllcflux, mp_source1, srf_speedrot, nx, ny, nz)
              rhllcflux = hllcflux

            case (4)                        !roe
              call rotateb(n, cleft, cleft_rot, angle1, angle2)
              call rotateb(n, cright, cright_rot, angle1, angle2)
              call rroe_riemann_solver(n, iconsidered, facex, cleft, cright, hllcflux, mp_source1, srf_speedrot, nx, ny, nz)
              rhllcflux = hllcflux

            case (5)                        !roe
              call rotateb(n, cleft, cleft_rot, angle1, angle2)
              call rotateb(n, cright, cright_rot, angle1, angle2)
              call troe_riemann_solver(n, cleft, cright, hllcflux, mp_source1, srf_speedrot, nx, ny, nz)
              rhllcflux = hllcflux

            end select

        godflux2(1:nof_variables) = godflux2(1:nof_variables) + (rhllcflux(1:nof_variables)*(weights_temp(ngp)*ielem(n, i)%surf(l)))
            if (multispecies .eq. 1) then
              mp_source2 = mp_source2 + mp_source1*(weights_temp(ngp)*ielem(n, i)%surf(l))
            end if
            if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
              if (icoupleturb .eq. 0) then        !first order upwind flux
                norms = (nx*(u_c(i)%val(1, 2)/u_c(i)%val(1, 1))) &
                        + (ny*(u_c(i)%val(1, 3)/u_c(i)%val(1, 1))) &
                        + (nz*(u_c(i)%val(1, 4)/u_c(i)%val(1, 1)))
                if (norms .ge. zero) then
                  rhllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=(norms)*cturbl(1:turbulenceequations+passivescalar)
                end if
                if (norms .lt. zero) then
                  rhllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=(norms)*cturbr(1:turbulenceequations+passivescalar)
                end if
              end if

              godflux2(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=godflux2(nof_variables+1:nof_variables+turbulenceequations+passivescalar)+&
             (rhllcflux(nof_variables + 1:nof_variables + turbulenceequations + passivescalar)*(weights_temp(ngp)*ielem(n, i)%surf(l)))
            end if
          end do

          rhs(i)%val(1:nof_variables) = rhs(i)%val(1:nof_variables) + godflux2(1:nof_variables)
          if (multispecies .eq. 1) then
            mp_source3 = mp_source3 + mp_source2
          end if

!         rhs(ielem(n,i)%ineigh(l))%val(1:nof_variables)=rhs(ielem(n,i)%ineigh(l))%val(1:nof_variables)-godflux2(1:nof_variables)
          if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
            rhst(i)%val(1:turbulenceequations + passivescalar) = rhst(i)%val(1:turbulenceequations + passivescalar) + &
            godflux2(nof_variables + 1:nof_variables + turbulenceequations + passivescalar)
!           rhst(ielem(n,i)%ineigh(l))%val(1:turbulenceequations+passivescalar)=rhst(ielem(n,i)%ineigh(l))%val(1:turbulenceequations+passivescalar)-&
!           godflux2(nof_variables+1:nof_variables+turbulenceequations+passivescalar)
          end if
!                                      end if
        end do
        if (multispecies .eq. 1) then
          rhs(i)%val(8) = rhs(i)%val(8) - (u_c(i)%val(1, 8)*mp_source3)
        end if
      end if
    end do
    !$omp end do

    !$omp barrier
    !$omp do
    do ii = 1, nof_bounded
      i = el_bnd(ii)
      iconsidered = i
      mp_source3 = zero
      if (ielem(n, i)%recalc .eq. 1) then

        do l = 1, ielem(n, i)%ifca
          igoflux = 0
          if (ielem(n, i)%ineighb(l) .eq. n) then        !my cpu only
            if (ielem(n, i)%ibounds(l) .gt. 0) then        !check for boundaries
              if (ibound(n, ielem(n, i)%ibounds(l))%icode .eq. 5) then        !periodic in my cpu
                icaseb = 1        !periodic mine
              else
                icaseb = 3        !physical
              end if
            else
              icaseb = 2!no boundaries interior
            end if
          else
            if (ielem(n, i)%ibounds(l) .gt. 0) then        !check for boundaries
              if (ibound(n, ielem(n, i)%ibounds(l))%icode .eq. 5) then        !periodic in other cpu
                icaseb = 4
              end if
            else
              icaseb = 5
            end if
          end if

          angle1 = ielem(n, i)%faceanglex(l)
          angle2 = ielem(n, i)%faceangley(l)
          nx = (cos(angle1)*sin(angle2))
          ny = (sin(angle1)*sin(angle2))
          nz = (cos(angle2))

          if (ielem(n, i)%types_faces(l) .eq. 5) then
            iqp = qp_quad_n
            weights_temp(1:iqp) = weights_q(1:iqp)
            n_node = 4
          else
            iqp = qp_triangle_n
            weights_temp(1:iqp) = weights_t(1:iqp)
            n_node = 3
          end if
          godflux2 = zero

          mp_source2 = zero

          do ngp = 1, iqp
            b_code = 0
            cleft(1:nof_variables) = ilocal_recon3(i)%uleft(1:nof_variables, l, ngp)
            if ((cascade .eq. 2) .and. (ielem(n, i)%mood .eq. 1)) then
              cleft(1:nof_variables) = u_c(i)%val(3, 1:nof_variables)
            end if

            if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
              if (icoupleturb .eq. 1) then
           cturbl(1:turbulenceequations + passivescalar) = ilocal_recon3(i)%uleftturb(1:turbulenceequations + passivescalar, l, ngp)
              else
                cturbl(1:turbulenceequations + passivescalar) = u_ct(i)%val(1, 1:turbulenceequations + passivescalar)
              end if
            end if

            if (ielem(n, i)%ineighb(l) .eq. n) then        !my cpu only
              if (ielem(n, i)%ibounds(l) .gt. 0) then        !check for boundaries
                if (ibound(n, ielem(n, i)%ibounds(l))%icode .eq. 5) then        !periodic in my cpu
                  cright(1:nof_variables) = ilocal_recon3(ielem(n, i)%ineigh(l))%uleft(1:nof_variables, ielem(n, i)%ineighn(l), ngp)
                  if ((cascade .eq. 2) .and. (ielem(n, (ielem(n, i)%ineigh(l)))%mood .eq. 1)) then
                    cright(1:nof_variables) = u_c(ielem(n, i)%ineigh(l))%val(3, 1:nof_variables)
                  end if
                  if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
                    if (icoupleturb .eq. 1) then
                    if (cascade .eq. 1) then
                      cturbr(1:turbulenceequations + passivescalar) = ilocal_recon3(ielem(n, i)%ineigh(l))%uleftturb &
                          (1:turbulenceequations + passivescalar, ielem(n, i)%ineighn(l), ngp)!right additional equations flow state
                    else
                      cturbr(1:turbulenceequations + passivescalar) = u_ct(ielem(n, i)%ineigh(l))%val(1, 1:turbulenceequations + passivescalar)
                    end if
                    else
                      cturbr(1:turbulenceequations + passivescalar) = u_ct(ielem(n, i)%ineigh(l))%val(1, 1:turbulenceequations + passivescalar)
                    end if
                  end if
                else
                  !not periodic ones in my cpu
                  facex = l; iconsidered = i
                  call coordinates_face_innerx(n, iconsidered, facex, vext, nodes_list)
                  if (ielem(n, iconsidered)%types_faces(facex) .eq. 5) then
                    n_node = 4
                  else
                    n_node = 3
                  end if

                  cords(1:3) = zero
                  cords(1:3) = cordinates3(n, nodes_list, n_node)

                  poy(1) = cords(2)
                  pox(1) = cords(1)
                  poz(1) = cords(3)

                  leftv(1:nof_variables) = cleft(1:nof_variables)
                  b_code = ibound(n, ielem(n, i)%ibounds(l))%icode

                  call boundarys(n,b_code,iconsidered,facex,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speed,srf_speedrot,ibfc)
                  cright(1:nof_variables) = rightv(1:nof_variables)

                end if
              else
                cright(1:nof_variables) = ilocal_recon3(ielem(n, i)%ineigh(l))%uleft(1:nof_variables, ielem(n, i)%ineighn(l), ngp)
                if ((cascade .eq. 2) .and. (ielem(n, (ielem(n, i)%ineigh(l)))%mood .eq. 1)) then
                  cright(1:nof_variables) = u_c(ielem(n, i)%ineigh(l))%val(3, 1:nof_variables)
                end if

                if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
                  if (icoupleturb .eq. 1) then
                    cturbr(1:turbulenceequations + passivescalar) = ilocal_recon3(ielem(n, i)%ineigh(l))%uleftturb &
                          (1:turbulenceequations + passivescalar, ielem(n, i)%ineighn(l), ngp)!right additional equations flow state
                  else
                    cturbr(1:turbulenceequations + passivescalar) = u_ct(ielem(n, i)%ineigh(l))%val(1, 1:turbulenceequations + passivescalar)
                  end if
                end if
              end if
            else        !in other cpus they can only be periodic or mpi neighbours
              if (ielem(n, i)%ibounds(l) .gt. 0) then        !check for boundaries
                if (ibound(n, ielem(n, i)%ibounds(l))%icode .eq. 5) then        !periodic in other cpu
                  cright(1:nof_variables) = iexboundhir(ielem(n, i)%ineighn(l))%facesol(ielem(n, i)%q_face(l)%q_mapl(ngp), 1:nof_variables)
                  if ((cascade .eq. 2) .and. (iexboundhir(ielem(n, i)%ineighn(l))%facesol_m(ielem(n, i)%q_face(l)%q_mapl(ngp), 1) .gt. 0.5)) then
                    cright(1:nof_variables) = iexsolhir(ilocal_recon3(i)%ihexn(1, ielem(n, i)%indexi(l)))%sol &
                      (ilocal_recon3(i)%ihexl(1, ielem(n, i)%indexi(l)), 1:nof_variables)
                  end if
                  if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
                    if (icoupleturb .eq. 1) then
                      cturbr(1:turbulenceequations + passivescalar) = iexboundhir(ielem(n, i)%ineighn(l))%facesol &
                      (ielem(n, i)%q_face(l)%q_mapl(ngp), nof_variables + 1:nof_variables + turbulenceequations + passivescalar)!right additional equations flow state
                    else
                      cturbr(1:turbulenceequations + passivescalar) = iexboundhir(ielem(n, i)%ineighn(l))%facesol &
                      (ielem(n, i)%q_face(l)%q_mapl(ngp), nof_variables + 1:nof_variables + turbulenceequations + passivescalar)!right additional equations flow state
                    end if
                  end if
                end if
              else

           cright(1:nof_variables) = iexboundhir(ielem(n, i)%ineighn(l))%facesol(ielem(n, i)%q_face(l)%q_mapl(ngp), 1:nof_variables)
            if ((cascade .eq. 2) .and. (iexboundhir(ielem(n, i)%ineighn(l))%facesol_m(ielem(n, i)%q_face(l)%q_mapl(ngp), 1) .gt. 0.5)) then
              cright(1:nof_variables) = iexsolhir(ilocal_recon3(i)%ihexn(1, ielem(n, i)%indexi(l)))%sol &
              (ilocal_recon3(i)%ihexl(1, ielem(n, i)%indexi(l)), 1:nof_variables)
                end if

                if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
                  if (icoupleturb .eq. 1) then
                    cturbr(1:turbulenceequations + passivescalar) = iexboundhir(ielem(n, i)%ineighn(l))%facesol &
                    (ielem(n, i)%q_face(l)%q_mapl(ngp), nof_variables + 1:nof_variables + turbulenceequations + passivescalar)!right additional equations flow state
                  else
                    cturbr(1:turbulenceequations + passivescalar) = iexboundhir(ielem(n, i)%ineighn(l))%facesol &
                    (ielem(n, i)%q_face(l)%q_mapl(ngp), nof_variables + 1:nof_variables + turbulenceequations + passivescalar)!right additional equations flow state
                  end if
                end if
              end if
            end if
            call rotatef(n, cright_rot, cright, angle1, angle2)        !rotate wrt to normalvector of face and solve 1d riemann problem
            call rotatef(n, cleft_rot, cleft, angle1, angle2)        !rotate wrt to normalvector of face and solve 1d riemann problem

            if ((lmach .eq. 1)) then    !application of the low mach number correction
              leftv(1:nof_variables) = cleft_rot(1:nof_variables); rightv(1:nof_variables) = cright_rot(1:nof_variables)
              call lmacht(n, leftv, rightv)
              cleft_rot(1:nof_variables) = leftv(1:nof_variables); cright_rot(1:nof_variables) = rightv(1:nof_variables);
            end if
            if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
              cleft_rot(nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = cturbl(1:turbulenceequations + passivescalar)
              cright_rot(nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = cturbr(1:turbulenceequations + passivescalar)
            end if

            select case (iriemann)

            case (1)                        !hllc
              call hllc_riemann_solver(n, iconsidered, facex, cleft_rot, cright_rot, hllcflux, mp_source1, srf_speedrot)
              call rotateb(n, rhllcflux, hllcflux, angle1, angle2)
              if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
              rhllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=hllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)
              end if

            case (2)                        !rusanov
              call rusanov_riemann_solver(n, iconsidered, facex, cleft_rot, cright_rot, hllcflux, mp_source1, srf_speedrot)
              call rotateb(n, rhllcflux, hllcflux, angle1, angle2)
              if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
              rhllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=hllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)
              end if

            case (3)                        !roe
              call rotateb(n, cleft, cleft_rot, angle1, angle2)
              call rotateb(n, cright, cright_rot, angle1, angle2)
              call roe_riemann_solver(n, iconsidered, facex, cleft, cright, hllcflux, mp_source1, srf_speedrot, nx, ny, nz)
              rhllcflux = hllcflux

            case (9)                        !hll
              call hll_riemann_solver(n, cleft_rot, cright_rot, hllcflux, mp_source1, srf_speedrot)
              call rotateb(n, rhllcflux, hllcflux, angle1, angle2)
              if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
              rhllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=hllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)
              end if

            case (4)                        !roe
              call rotateb(n, cleft, cleft_rot, angle1, angle2)
              call rotateb(n, cright, cright_rot, angle1, angle2)
              if (b_code .le. 0) then
                call rroe_riemann_solver(n, iconsidered, facex, cleft, cright, hllcflux, mp_source1, srf_speedrot, nx, ny, nz)
              else
                call roe_riemann_solver(n, iconsidered, facex, cleft, cright, hllcflux, mp_source1, srf_speedrot, nx, ny, nz)
              end if
              rhllcflux = hllcflux

            case (5)                        !roe
              call rotateb(n, cleft, cleft_rot, angle1, angle2)
              call rotateb(n, cright, cright_rot, angle1, angle2)
              if (b_code .le. 0) then
                call troe_riemann_solver(n, cleft, cright, hllcflux, mp_source1, srf_speedrot, nx, ny, nz)
              else
                call roe_riemann_solver(n, iconsidered, facex, cleft, cright, hllcflux, mp_source1, srf_speedrot, nx, ny, nz)
              end if
              rhllcflux = hllcflux

            end select

        godflux2(1:nof_variables) = godflux2(1:nof_variables) + (rhllcflux(1:nof_variables)*(weights_temp(ngp)*ielem(n, i)%surf(l)))
            if (multispecies .eq. 1) then
              mp_source2 = mp_source2 + mp_source1*(weights_temp(ngp)*ielem(n, i)%surf(l))
            end if

            if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
              if (icoupleturb .eq. 0) then        !first order upwind flux
                norms = (nx*(u_c(i)%val(1, 2)/u_c(i)%val(1, 1))) &
                        + (ny*(u_c(i)%val(1, 3)/u_c(i)%val(1, 1))) &
                        + (nz*(u_c(i)%val(1, 4)/u_c(i)%val(1, 1)))
                if (norms .ge. zero) then
                  rhllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=(norms)*cturbl(1:turbulenceequations+passivescalar)
                end if
                if (norms .lt. zero) then
                  rhllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=(norms)*cturbr(1:turbulenceequations+passivescalar)
                end if
              end if
               godflux2(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=godflux2(nof_variables+1:nof_variables+turbulenceequations+passivescalar)+&
              (rhllcflux(nof_variables + 1:nof_variables + turbulenceequations + passivescalar)*(weights_temp(ngp)*ielem(n, i)%surf(l)))
            end if
          end do

          rhs(i)%val(1:nof_variables) = rhs(i)%val(1:nof_variables) + godflux2(1:nof_variables)
          if (multispecies .eq. 1) then
            mp_source3 = mp_source3 + mp_source2
          end if
!         if ((igoflux.eq.1))then
!         rhs(ielem(n,i)%ineigh(l))%val(1:nof_variables)=rhs(ielem(n,i)%ineigh(l))%val(1:nof_variables)-godflux2(1:nof_variables)
!         end if

          if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
            rhst(i)%val(1:turbulenceequations + passivescalar) = rhst(i)%val(1:turbulenceequations + passivescalar) + &
            godflux2(nof_variables + 1:nof_variables + turbulenceequations + passivescalar)
!          if ((igoflux.eq.1))then
!          rhst(ielem(n,i)%ineigh(l))%val(1:turbulenceequations+passivescalar)=rhst(ielem(n,i)%ineigh(l))%val(1:turbulenceequations+passivescalar)-&
!          godflux2(nof_variables+1:nof_variables+turbulenceequations+passivescalar)
!          end if
          end if
!                                     end if
        end do
        if (multispecies .eq. 1) then
          rhs(i)%val(8) = rhs(i)%val(8) - (u_c(i)%val(1, 8)*mp_source3)
        end if
      end if
    end do
    !$omp end do

  end subroutine calculate_fluxeshi_convective_mood

  subroutine calculate_fluxeshi_convective2d_mood(n)
!> @brief
!> this subroutine computes the convective fluxes for hyperbolic conservation laws in 2d
    implicit none
    integer, intent(in)::n
    real, dimension(1:nof_variables + turbulenceequations + passivescalar)::godflux2, dg_vol_rec, rhllcflux, hllcflux
    integer::i, l, ngp, kmaxe, iqp, ii, ikas, igoflux, icaseb, kxk, b_code
    real::sum_detect, norms
    real, dimension(1:numberofpoints2)::weights_temp, weights_dg
    integer::iconsidered, facex, pointx, n_node
    real::angle1, angle2, nx, ny, nz, mp_source1, mp_source2, mp_source3
    real, dimension(1:nof_variables + turbulenceequations + passivescalar)::cleft, cright, cleft_rot, cright_rot
    real, dimension(1:nof_variables)::leftv, rightv, srf_speedrot
    real, dimension(1:turbulenceequations + passivescalar)::cturbl, cturbr
    real, dimension(1:dimensiona)::pox, poy, poz
    real, dimension(1:nof_variables)::srf_speed
    real, dimension(1:8, 1:dimensiona)::vext, nodes_list
    real, dimension(1:dimensiona, 1:numberofpoints2)::qpoints2d
    real, dimension(1:numberofpoints2)::wequa2d
    real, dimension(1:dimensiona)::cords
    integer::ibfc
    kmaxe = xmpielrank(n)

    call quadratureline(n, igqrules, vext, qpoints2d, wequa2d)
    weights_temp = wequa2d(1:qp_line_n)

    if (reduce_comp .eq. 1) then
      weights_temp = 1.0d0
    end if

    do i = 1, kmaxe
    if (ielem(n, i)%recalc .eq. 1) then
      rhs(i)%val(:) = zero; if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) rhst(i)%val(:) = zero
    end if
    end do

    !$omp barrier
    !$omp do
    do ii = 1, nof_interior        !for all the interior elements
      i = el_int(ii)
      iconsidered = i
!                     rhs(i)%val(:)=zero;if ((turbulence.eq.1).or.(passivescalar.gt.0)) rhst(i)%val(:)=zero
      if (ielem(n, i)%recalc .eq. 1) then
        mp_source3 = zero
        do l = 1, ielem(n, i)%ifca !for all their faces

          godflux2 = zero
          mp_source2 = zero
          nx = ielem(n, i)%faceanglex(l)
          ny = ielem(n, i)%faceangley(l)
          angle1 = nx
          angle2 = ny
          b_code = 0
          iqp = qp_line_n

          do ngp = 1, iqp        !for all the gaussian quadrature points
            cleft(1:nof_variables) = ilocal_recon3(i)%uleft(1:nof_variables, l, ngp)        !left mean flow state
            if ((cascade .eq. 2) .and. (ielem(n, i)%mood .eq. 1)) then
              cleft(1:nof_variables) = u_c(i)%val(3, 1:nof_variables)
            end if
            cright(1:nof_variables) = ilocal_recon3(ielem(n, i)%ineigh(l))%uleft(1:nof_variables, ielem(n, i)%ineighn(l), ngp) !right mean flow state
            if ((cascade .eq. 2) .and. (ielem(n, (ielem(n, i)%ineigh(l)))%mood .eq. 1)) then
              cright(1:nof_variables) = u_c(ielem(n, i)%ineigh(l))%val(3, 1:nof_variables)
            end if
            if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
              if (icoupleturb .eq. 1) then
cturbl(1:turbulenceequations + passivescalar) = ilocal_recon3(i)%uleftturb(1:turbulenceequations + passivescalar, l, ngp) !left additional equations flow state
                turbr(1:turbulenceequations+passivescalar)=ilocal_recon3(ielem(n,i)%ineigh(l))%uleftturb(1:turbulenceequations+passivescalar,ielem(n,i)%ineighn(l),ngp)!right additional equations flow state
              else
                cturbl(1:turbulenceequations + passivescalar) = u_ct(i)%val(1, 1:turbulenceequations + passivescalar)
                cturbr(1:turbulenceequations + passivescalar) = u_ct(ielem(n, i)%ineigh(l))%val(1, 1:turbulenceequations + passivescalar)
              end if

              cleft_rot(nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = cturbl(1:turbulenceequations + passivescalar)
              cright_rot(nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = cturbr(1:turbulenceequations + passivescalar)
            end if

            call rotatef2d(n, cright_rot, cright, angle1, angle2)        !rotate wrt to normalvector of face and solve 1d riemann problem
            call rotatef2d(n, cleft_rot, cleft, angle1, angle2)        !rotate wrt to normalvector of face and solve 1d riemann problem

            if ((lmach .eq. 1)) then    !application of the low mach number correction
              leftv(1:nof_variables) = cleft_rot(1:nof_variables); rightv(1:nof_variables) = cright_rot(1:nof_variables)
              call lmacht2d(n, leftv, rightv)
              cleft_rot(1:nof_variables) = leftv(1:nof_variables); cright_rot(1:nof_variables) = rightv(1:nof_variables);
            end if

            select case (iriemann)

            case (1)                        !hllc
              call hllc_riemann_solver2d(n, cleft_rot, cright_rot, hllcflux, mp_source1, srf_speedrot)
              call rotateb2d(n, rhllcflux, hllcflux, angle1, angle2)
              if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
              rhllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=hllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)
              end if

            case (2)                        !rusanov
              call rusanov_riemann_solver2d(n, cleft_rot, cright_rot, hllcflux, mp_source1, srf_speedrot)
              call rotateb2d(n, rhllcflux, hllcflux, angle1, angle2)
              if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
              rhllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=hllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)
              end if

            case (3)                        !roe
              call rotateb2d(n, cleft, cleft_rot, angle1, angle2)
              call rotateb2d(n, cright, cright_rot, angle1, angle2)
              call roe_riemann_solver2d(n, cleft, cright, hllcflux, mp_source1, srf_speedrot, nx, ny)
              rhllcflux(1:nof_variables) = hllcflux(1:nof_variables)

            case (9)                        !hll
              call hll_riemann_solver2d(n, cleft_rot, cright_rot, hllcflux, mp_source1, srf_speedrot)
              call rotateb2d(n, rhllcflux, hllcflux, angle1, angle2)
              if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
              rhllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=hllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)
              end if

            case (4)                        !roe
              call rotateb2d(n, cleft, cleft_rot, angle1, angle2)
              call rotateb2d(n, cright, cright_rot, angle1, angle2)
              call rroe_riemann_solver2d(n, cleft, cright, hllcflux, mp_source1, srf_speedrot, nx, ny, b_code)
              rhllcflux = hllcflux

            end select

        godflux2(1:nof_variables) = godflux2(1:nof_variables) + (rhllcflux(1:nof_variables)*(weights_temp(ngp)*ielem(n, i)%surf(l)))
          if (multispecies .eq. 1) then
            mp_source2 = mp_source2 + mp_source1*(weights_temp(ngp)*ielem(n, i)%surf(l))
          end if
            if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
              if (icoupleturb .eq. 0) then        !first order upwind flux
                norms = (nx*(u_c(i)%val(1, 2)/u_c(i)%val(1, 1))) &
                        + (ny*(u_c(i)%val(1, 3)/u_c(i)%val(1, 1)))
                if (norms .ge. zero) then
                  rhllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=(norms)*cturbl(1:turbulenceequations+passivescalar)
                end if
                if (norms .lt. zero) then
                  rhllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=(norms)*cturbr(1:turbulenceequations+passivescalar)
                end if
              end if
                godflux2(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=godflux2(nof_variables+1:nof_variables+turbulenceequations+passivescalar)+&
                (rhllcflux(nof_variables + 1:nof_variables + turbulenceequations + passivescalar)*(weights_temp(ngp)*ielem(n, i)%surf(l)))
            end if

          end do
          rhs(i)%val(1:nof_variables) = rhs(i)%val(1:nof_variables) + godflux2(1:nof_variables)
          if (multispecies .eq. 1) then
            mp_source3 = mp_source3 + mp_source2
          end if

!              rhs(ielem(n,i)%ineigh(l))%val(1:nof_variables)=rhs(ielem(n,i)%ineigh(l))%val(1:nof_variables)-godflux2(1:nof_variables)
          if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
            rhst(i)%val(1:turbulenceequations + passivescalar) = rhst(i)%val(1:turbulenceequations + passivescalar) + &
            godflux2(nof_variables + 1:nof_variables + turbulenceequations + passivescalar)
!                                     rhst(ielem(n,i)%ineigh(l))%val(1:turbulenceequations+passivescalar)=rhst(ielem(n,i)%ineigh(l))%val(1:turbulenceequations+passivescalar)-&
!                                      godflux2(nof_variables+1:nof_variables+turbulenceequations+passivescalar)
!
          end if
!                                     end if
        end do
        if (multispecies .eq. 1) then
          rhs(i)%val(7) = rhs(i)%val(7) - (u_c(i)%val(1, 7)*mp_source3)!*ielem(n,i)%totvolume)

        end if

      end if
    end do
    !$omp end do

    !$omp do
    do ii = 1, nof_bounded
      i = el_bnd(ii)
      iconsidered = i
      mp_source3 = zero
      if (ielem(n, i)%recalc .eq. 1) then

        do l = 1, ielem(n, i)%ifca
          igoflux = 0
          if (ielem(n, i)%ineighb(l) .eq. n) then        !my cpu only
            if (ielem(n, i)%ibounds(l) .gt. 0) then        !check for boundaries
              if (ibound(n, ielem(n, i)%ibounds(l))%icode .eq. 5) then        !periodic in my cpu
                icaseb = 1        !periodic mine
              else
                icaseb = 3        !physical
              end if
            else
              icaseb = 2!no boundaries interior
            end if
          else
            if (ielem(n, i)%ibounds(l) .gt. 0) then        !check for boundaries
              if (ibound(n, ielem(n, i)%ibounds(l))%icode .eq. 5) then        !periodic in other cpu
                icaseb = 4
              end if
            else
              icaseb = 5
            end if
          end if

          b_code = 0
          nx = ielem(n, i)%faceanglex(l)
          ny = ielem(n, i)%faceangley(l)
          angle1 = nx
          angle2 = ny

          iqp = qp_line_n
          godflux2 = zero
          mp_source2 = zero
          do ngp = 1, iqp
            cleft(1:nof_variables) = ilocal_recon3(i)%uleft(1:nof_variables, l, ngp)
            if ((cascade .eq. 2) .and. (ielem(n, i)%mood .eq. 1)) then
              cleft(1:nof_variables) = u_c(i)%val(3, 1:nof_variables)
            end if
            if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
              if (icoupleturb .eq. 1) then
           cturbl(1:turbulenceequations + passivescalar) = ilocal_recon3(i)%uleftturb(1:turbulenceequations + passivescalar, l, ngp)
              else
                cturbl(1:turbulenceequations + passivescalar) = u_ct(i)%val(1, 1:turbulenceequations + passivescalar)
              end if
            end if

            if (ielem(n, i)%ineighb(l) .eq. n) then        !my cpu only
              if (ielem(n, i)%ibounds(l) .gt. 0) then        !check for boundaries
                if (ibound(n, ielem(n, i)%ibounds(l))%icode .eq. 5) then        !periodic in my cpu
                  cright(1:nof_variables) = ilocal_recon3(ielem(n, i)%ineigh(l))%uleft(1:nof_variables, ielem(n, i)%ineighn(l), ngp)

                  if ((cascade .eq. 2) .and. (ielem(n, (ielem(n, i)%ineigh(l)))%mood .eq. 1)) then
                    cright(1:nof_variables) = u_c(ielem(n, i)%ineigh(l))%val(3, 1:nof_variables)

                  end if

                  if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
                    if (icoupleturb .eq. 1) then
                      cturbr(1:turbulenceequations + passivescalar) = ilocal_recon3(ielem(n, i)%ineigh(l))%uleftturb &
                          (1:turbulenceequations + passivescalar, ielem(n, i)%ineighn(l), ngp)!right additional equations flow state
                    else
                      cturbr(1:turbulenceequations + passivescalar) = u_ct(ielem(n, i)%ineigh(l))%val(1, 1:turbulenceequations + passivescalar)
                    end if
                  end if

                  ikas = 1

                else
                  !not periodic ones in my cpu
                  facex = l; iconsidered = i
                  call coordinates_face_inner2dx(n, iconsidered, facex, vext, nodes_list)
                  cords(1:2) = zero
                  n_node = 2
                  cords(1:2) = cordinates2(n, nodes_list, n_node)

                  poy(1) = cords(2)
                  pox(1) = cords(1)

                  leftv(1:nof_variables) = cleft(1:nof_variables)
                  b_code = ibound(n, ielem(n, i)%ibounds(l))%icode
                  call boundarys2d(n,b_code,iconsidered,facex,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speed,srf_speedrot,ibfc)
                  cright(1:nof_variables) = rightv(1:nof_variables)
                  ikas = 2

                end if
              else
                cright(1:nof_variables) = ilocal_recon3(ielem(n, i)%ineigh(l))%uleft(1:nof_variables, ielem(n, i)%ineighn(l), ngp)
                if ((cascade .eq. 2) .and. (ielem(n, (ielem(n, i)%ineigh(l)))%mood .eq. 1)) then
                  cright(1:nof_variables) = u_c(ielem(n, i)%ineigh(l))%val(3, 1:nof_variables)
                end if
                if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
                  if (icoupleturb .eq. 1) then
                    cturbr(1:turbulenceequations + passivescalar) = ilocal_recon3(ielem(n, i)%ineigh(l))%uleftturb &
                          (1:turbulenceequations + passivescalar, ielem(n, i)%ineighn(l), ngp)!right additional equations flow state
                  else
                    cturbr(1:turbulenceequations + passivescalar) = u_ct(ielem(n, i)%ineigh(l))%val(1, 1:turbulenceequations + passivescalar)
                  end if
                end if
                ikas = 3
              end if
            else        !in other cpus they can only be periodic or mpi neighbours

              if (ielem(n, i)%ibounds(l) .gt. 0) then        !check for boundaries
                if (ibound(n, ielem(n, i)%ibounds(l))%icode .eq. 5) then        !periodic in other cpu
                  cright(1:nof_variables) = iexboundhir(ielem(n, i)%ineighn(l))%facesol(ielem(n, i)%q_face(l)%q_mapl(ngp), 1:nof_variables)

                if ((cascade .eq. 2) .and. (iexboundhir(ielem(n, i)%ineighn(l))%facesol_m(ielem(n, i)%q_face(l)%q_mapl(ngp), 1) .gt. 0.5)) then
                    cright(1:nof_variables) = iexsolhir(ilocal_recon3(i)%ihexn(1, ielem(n, i)%indexi(l)))%sol &
                      (ilocal_recon3(i)%ihexl(1, ielem(n, i)%indexi(l)), 1:nof_variables)
                  end if

                  if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
                    if (icoupleturb .eq. 1) then
                      cturbr(1:turbulenceequations + passivescalar) = iexboundhir(ielem(n, i)%ineighn(l))%facesol &
                      (ielem(n, i)%q_face(l)%q_mapl(ngp), nof_variables + 1:nof_variables + turbulenceequations + passivescalar)!right additional equations flow state
                    else
                      cturbr(1:turbulenceequations + passivescalar) = iexboundhir(ielem(n, i)%ineighn(l))%facesol &
                      (ielem(n, i)%q_face(l)%q_mapl(ngp), nof_variables + 1:nof_variables + turbulenceequations + passivescalar)!right additional equations flow state
                    end if
                  end if
                end if
              else

           cright(1:nof_variables) = iexboundhir(ielem(n, i)%ineighn(l))%facesol(ielem(n, i)%q_face(l)%q_mapl(ngp), 1:nof_variables)
             if ((cascade .eq. 2) .and. (iexboundhir(ielem(n, i)%ineighn(l))%facesol_m(ielem(n, i)%q_face(l)%q_mapl(ngp), 1) .gt. 0.5)) then
               cright(1:nof_variables) = iexsolhir(ilocal_recon3(i)%ihexn(1, ielem(n, i)%indexi(l)))%sol &
               (ilocal_recon3(i)%ihexl(1, ielem(n, i)%indexi(l)), 1:nof_variables)
                end if

                if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
                  if (icoupleturb .eq. 1) then
                    cturbr(1:turbulenceequations + passivescalar) = iexboundhir(ielem(n, i)%ineighn(l))%facesol &
                    (ielem(n, i)%q_face(l)%q_mapl(ngp), nof_variables + 1:nof_variables + turbulenceequations + passivescalar)!right additional equations flow state
                  else
                    cturbr(1:turbulenceequations + passivescalar) = iexboundhir(ielem(n, i)%ineighn(l))%facesol &
                    (ielem(n, i)%q_face(l)%q_mapl(ngp), nof_variables + 1:nof_variables + turbulenceequations + passivescalar)!right additional equations flow state
                  end if
                end if
                ikas = 4
              end if
            end if

            call rotatef2d(n, cright_rot, cright, angle1, angle2)        !rotate wrt to normalvector of face and solve 1d riemann problem
            call rotatef2d(n, cleft_rot, cleft, angle1, angle2)        !rotate wrt to normalvector of face and solve 1d riemann problem

            if ((lmach .eq. 1)) then    !application of the low mach number correction
              leftv(1:nof_variables) = cleft_rot(1:nof_variables); rightv(1:nof_variables) = cright_rot(1:nof_variables)
              call lmacht2d(n, leftv, rightv)
              cleft_rot(1:nof_variables) = leftv(1:nof_variables); cright_rot(1:nof_variables) = rightv(1:nof_variables);
            end if
            if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
              cleft_rot(nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = cturbl(1:turbulenceequations + passivescalar)
              cright_rot(nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = cturbr(1:turbulenceequations + passivescalar)
            end if

            select case (iriemann)

            case (1)                        !hllc
              call hllc_riemann_solver2d(n, cleft_rot, cright_rot, hllcflux, mp_source1, srf_speedrot)
              call rotateb2d(n, rhllcflux, hllcflux, angle1, angle2)
              if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
              rhllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=hllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)
              end if

            case (2)                        !rusanov
              call rusanov_riemann_solver2d(n, cleft_rot, cright_rot, hllcflux, mp_source1, srf_speedrot)
              call rotateb2d(n, rhllcflux, hllcflux, angle1, angle2)
              if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
              rhllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=hllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)
              end if

            case (9)                        !hll
              call hll_riemann_solver2d(n, cleft_rot, cright_rot, hllcflux, mp_source1, srf_speedrot)
              call rotateb2d(n, rhllcflux, hllcflux, angle1, angle2)
              if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
              rhllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=hllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)
              end if

            case (3)                        !roe
              call rotateb2d(n, cleft, cleft_rot, angle1, angle2)
              call rotateb2d(n, cright, cright_rot, angle1, angle2)
              call roe_riemann_solver2d(n, cleft, cright, hllcflux, mp_source1, srf_speedrot, nx, ny)
              rhllcflux = hllcflux

            case (4)                        !roe
              if ((b_code .le. 0)) then
                call rotateb2d(n, cleft, cleft_rot, angle1, angle2)
                call rotateb2d(n, cright, cright_rot, angle1, angle2)
                call rroe_riemann_solver2d(n, cleft, cright, hllcflux, mp_source1, srf_speedrot, nx, ny, b_code)
                rhllcflux = hllcflux
              else
                call rusanov_riemann_solver2d(n, cleft, cright, hllcflux, mp_source1, srf_speedrot)
                call rotateb2d(n, rhllcflux, hllcflux, angle1, angle2)
              end if

            end select

        godflux2(1:nof_variables) = godflux2(1:nof_variables) + (rhllcflux(1:nof_variables)*(weights_temp(ngp)*ielem(n, i)%surf(l)))
            if (multispecies .eq. 1) then
              mp_source2 = mp_source2 + mp_source1*(weights_temp(ngp)*ielem(n, i)%surf(l))
            end if

            if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
              if (icoupleturb .eq. 0) then        !first order upwind flux

                norms = (nx*(u_c(i)%val(1, 2)/u_c(i)%val(1, 1))) &
                        + (ny*(u_c(i)%val(1, 3)/u_c(i)%val(1, 1)))

                if (norms .ge. zero) then
                  rhllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=(norms)*cturbl(1:turbulenceequations+passivescalar)
                end if
                if (norms .lt. zero) then
                  rhllcflux(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=(norms)*cturbr(1:turbulenceequations+passivescalar)
                end if
              end if

              if ((b_code .eq. 4) .or. (b_code .eq. 3)) then
                rhllcflux(nof_variables + 1:nof_variables + turbulenceequations + passivescalar) = zero
              end if
                godflux2(nof_variables+1:nof_variables+turbulenceequations+passivescalar)=godflux2(nof_variables+1:nof_variables+turbulenceequations+passivescalar)+&
                (rhllcflux(nof_variables + 1:nof_variables + turbulenceequations + passivescalar)*(weights_temp(ngp)*ielem(n, i)%surf(l)))
            end if
          end do

          rhs(i)%val(1:nof_variables) = rhs(i)%val(1:nof_variables) + godflux2(1:nof_variables)
          if (multispecies .eq. 1) then
            mp_source3 = mp_source3 + mp_source2
          end if
!                                     if ((igoflux.eq.1))then
!                                     rhs(ielem(n,i)%ineigh(l))%val(1:nof_variables)=rhs(ielem(n,i)%ineigh(l))%val(1:nof_variables)-godflux2(1:nof_variables)
!                                     end if

          if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
            rhst(i)%val(1:turbulenceequations + passivescalar) = rhst(i)%val(1:turbulenceequations + passivescalar) + &
            godflux2(nof_variables + 1:nof_variables + turbulenceequations + passivescalar)
!                                      if ((igoflux.eq.1))then
!                                      rhst(ielem(n,i)%ineigh(l))%val(1:turbulenceequations+passivescalar)=rhst(ielem(n,i)%ineigh(l))%val(1:turbulenceequations+passivescalar)-&
!                                     godflux2(nof_variables+1:nof_variables+turbulenceequations+passivescalar)
!                                     end if
!                                     end if
          end if
        end do
        if (multispecies .eq. 1) then
          rhs(i)%val(7) = rhs(i)%val(7) - (u_c(i)%val(1, 7)*mp_source3)!*ielem(n,i)%totvolume)
        end if
      end if
    end do
    !$omp end do

  end subroutine calculate_fluxeshi_convective2d_mood

end module fluxes
