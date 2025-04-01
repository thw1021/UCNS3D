module implicit_fluxes
  use declaration
  use library
  use transform
  use local
  use riemann
  use flow_operations
  use source
  implicit none

contains
  subroutine calculate_jacobian(n)
    !> @brief
    !> this subroutine computes the approximate jacobian for implicit time stepping in 3d
    implicit none
    integer, intent(in)::n
    real, dimension(1:nof_variables + turbulenceequations + passivescalar)::godflux2
    integer::i, l, ngp, kmaxe, iqp, ii, nvar, n_node, ibfc
    real::sum_detect, norms, vpp, asound1, asound2, mul1, dxb, tempxx, viscots
    real, dimension(nof_variables, nof_variables)::identity1
    real, dimension(nof_variables, nof_variables)::convj, diffj
    integer::iconsidered, facex, pointx, igoflux
    integer::b_code, srf
    real::angle1, angle2, nx, ny, nz
    real, dimension(1:nof_variables)::cleft, cright, cright_rot, cleft_rot
    real, dimension(1:turbulenceequations + passivescalar)::cturbl, cturbr
    real, dimension(1:nof_variables)::leftv, srf_speedrot, srf_speed
    real, dimension(1:nof_variables)::rightv
    real, dimension(1:dimensiona)::pox, poy, poz
    real, dimension(1:8, 1:dimensiona)::vext
    real, dimension(1:8, 1:dimensiona)::nodes_list
    real, dimension(1:dimensiona)::cords
    real, dimension(1:4)::viscl, laml
    real, dimension(1:2)::turbmv
    real, dimension(1)::etvm
    real, dimension(1:20)::eddyfl, eddyfr
    real::mp_pinfl, gammal
    real::mp_pinfr, gammar
    real, dimension(1:nof_variables, 1:nof_variables)::eigvl

    kmaxe = xmpielrank(n)
    identity1(:, :) = zero
    do l = 1, nof_variables
      identity1(l, l) = 1.0d0
    end do

    !$omp do
    do ii = 1, nof_interior        !for all the interior elements
      i = el_int(ii)
      iconsidered = i
      impdiag(i, :, :) = zero
      impoff(i, :, :, :) = zero
      if (turbulence .eq. 1) then
        impdiagt(i, :) = zero
        impofft(i, :, :) = zero
      end if

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
        mul1 = ielem(n, i)%surf(l)
        b_code = 0
        cleft(1:nof_variables) = u_c(i)%val(1, 1:nof_variables)
        cright(1:nof_variables) = u_c(ielem(n, i)%ineigh(l))%val(1, 1:nof_variables)

        if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

          cturbl(1:turbulenceequations + passivescalar) = u_ct(i)%val(1, 1:turbulenceequations + passivescalar)!left additional equations flow state
          cturbr(1:turbulenceequations + passivescalar) = u_ct(ielem(n, i)%ineigh(l))%val(1, 1:turbulenceequations + passivescalar)

        end if

        call rotatef(n, cright_rot, cright, angle1, angle2)        !rotate wrt to normalvector of face and solve 1d riemann problem
        call rotatef(n, cleft_rot, cleft, angle1, angle2)        !rotate wrt to normalvector of face and solve 1d riemann problem

        if ((lmach .eq. 1)) then    !application of the low mach number correction
          leftv(1:nof_variables) = cleft_rot(1:nof_variables); rightv(1:nof_variables) = cright_rot(1:nof_variables)
          call lmacht(n, leftv, rightv)
          cleft_rot(1:nof_variables) = leftv(1:nof_variables); cright_rot(1:nof_variables) = rightv(1:nof_variables);
          call rotateb(n, cright, cright_rot, angle1, angle2)
          call rotateb(n, cleft, cleft_rot, angle1, angle2)

        end if

        leftv(1:nof_variables) = cleft(1:nof_variables); rightv(1:nof_variables) = cright(1:nof_variables)
        call cons2prim2(n, leftv, rightv, mp_pinfl, mp_pinfr, gammal, gammar)

        asound1 = sqrt(leftv(5)*gamma/leftv(1)) + abs(cleft_rot(2)/cleft_rot(1))
        asound2 = sqrt(rightv(5)*gamma/rightv(1)) + abs(cright_rot(2)/cright_rot(1))
        if (ilocal_recon3(i)%mrf .eq. 1) then
          !retrieve the rotational velocity (at the gaussian point just for second order)
          srf_speed(2:4) = ilocal_recon3(i)%rotvel(l, 1, 1:3)
          call rotatef(n, srf_speedrot, srf_speed, angle1, angle2)
          !calculate the new eigenvalue for rotating reference frame
          asound1 = sqrt(leftv(5)*gamma/leftv(1)) + abs(cleft_rot(2)/cleft_rot(1) - srf_speedrot(2))
          asound2 = sqrt(rightv(5)*gamma/rightv(1)) + abs(cright_rot(2)/cright_rot(1) - srf_speedrot(2))
        end if
        vpp = max(asound1, asound2)

        if (itestcase .eq. 4) then
          call sutherland(n, leftv, rightv, viscl, laml)

          viscots = (viscl(1) + viscl(2))*oo2
          mul1 = ielem(n, i)%surf(l)

          viscots = max((4.0/(3.0*0.5*(cleft(1) + cright(1))))*viscots*mul1 &
                        /ielem(n, i)%dih(l), ((gamma/(0.5*(cleft(1) + cright(1)))) &
                                              *viscots*mul1/(ielem(n, i)%dih(l)*prandtl)))
          vpp = max(asound1, asound2) + viscots

          if (turbulence .eq. 1) then
            if (turbulencemodel .eq. 1) then
              turbmv(1) = cturbl(1); turbmv(2) = cturbr(1); eddyfl(2) = turbmv(1); eddyfr(2) = turbmv(2)
              call eddyvisco(n, viscl, laml, turbmv, etvm, eddyfl, eddyfr, leftv, rightv)
            end if
            if (turbulencemodel .eq. 2) then
              eddyfl(1) = ielem(n, i)%walldist; eddyfl(2) = cturbl(1); eddyfl(3) = cturbl(2)

              eddyfl(4:6) = ilocal_recon3(i)%grads(1, 1:3); eddyfl(7:9) = ilocal_recon3(i)%grads(2, 1:3)
              eddyfl(10:12) = ilocal_recon3(i)%grads(3, 1:3); eddyfl(13:15) = ilocal_recon3(i)%grads(5, 1:3)
              eddyfl(16:18) = ilocal_recon3(i)%grads(6, 1:3)

              eddyfr = eddyfl
              call eddyvisco(n, viscl, laml, turbmv, etvm, eddyfl, eddyfr, leftv, rightv)
            end if

            viscots = oo2*((viscl(1) + viscl(3)) + (viscl(2) + viscl(4)))

            viscots = max((4.0/(3.0*0.5*(cleft(1) + cright(1))))*viscots* &
                          mul1/ielem(n, i)%dih(l), ((gamma/(0.5*(cleft(1) + cright(1))))* &
                                                    viscots*mul1/(ielem(n, i)%dih(l)*(prandtl + prtu))))

            vpp = max(asound1, asound2) + viscots
          end if

          impdiag(i, 1:nof_variables, 1:nof_variables) = impdiag(i, 1:nof_variables, 1:nof_variables) + (oo2*((vpp*identity1))*mul1)
          call compute_jacobianse(n, iconsidered, eigvl, cright, gamma, angle1, angle2, srf_speedrot, nx, ny, nz)
          convj = eigvl
                                                  impoff(i,l,1:nof_variables,1:nof_variables)=impoff(i,l,1:nof_variables,1:nof_variables)+(((oo2*convj(1:nof_variables,1:nof_variables))&
                                                                                                    - ((oo2*vpp)*identity1))*mul1)

          if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
            if (turbulence .eq. 1) then
            do nvar = 1, turbulenceequations

              vpp = max(asound1, asound2) + viscots
              impdiagt(i, nvar) = impdiagt(i, nvar) + (oo2*((vpp))*mul1)
              impofft(i, l, nvar) = impofft(i, l, nvar) - (((oo2*vpp))*mul1)
            end do
            end if
            if (passivescalar .gt. 0) then
            do nvar = turbulenceequations + 1, turbulenceequations + passivescalar
              viscl(1) = viscl(1)/schmidt_lam
              viscl(2) = viscl(2)/schmidt_lam
              if (turbulence .eq. 1) then
                viscl(3) = viscl(3)/schmidt_turb
                viscl(4) = viscl(4)/schmidt_turb
                viscots = 0.5*((viscl(1) + viscl(3)) + (viscl(2) + viscl(4)))
              else
                viscots = 0.5*((viscl(1) + viscl(2)))

              end if

              viscots = (2.0*viscots)/((cleft(1) + cright(1))*ielem(n, i)%dih(l))
              vpp = max(asound1, asound2) + viscots
              impdiagt(i, nvar) = impdiagt(i, nvar) + (oo2*((vpp))*mul1)
              impofft(i, l, nvar) = impofft(i, l, nvar) - (oo2*((vpp))*mul1)
            end do
            end if
          end if
        else

          impdiag(i, 1:nof_variables, 1:nof_variables) = impdiag(i, 1:nof_variables, 1:nof_variables) + (oo2*((vpp*identity1))*mul1)
          call compute_jacobianse(n, iconsidered, eigvl, cright, gamma, angle1, angle2, srf_speedrot, nx, ny, nz)
          convj = eigvl
          impoff(i,l,1:nof_variables,1:nof_variables)=impoff(i,l,1:nof_variables,1:nof_variables)+(((oo2*convj(1:nof_variables,1:nof_variables))&
                                                                                                    - ((oo2*vpp)*identity1))*mul1)

        end if

      end do
    end do
    !$omp end do

    !$omp do
    do ii = 1, nof_bounded
      i = el_bnd(ii)
      iconsidered = i

      impdiag(i, :, :) = 0.0
      impoff(i, :, :, :) = 0.0
      if (turbulence .eq. 1) then
        impdiagt(i, :) = 0.0
        impofft(i, :, :) = 0.0
      end if
      if (mrf .eq. 1) then
        srf = ilocal_recon3(i)%mrf
      end if
      do l = 1, ielem(n, i)%ifca
        b_code = 0
        angle1 = ielem(n, i)%faceanglex(l)
        angle2 = ielem(n, i)%faceangley(l)
        nx = (cos(angle1)*sin(angle2))
        ny = (sin(angle1)*sin(angle2))
        nz = (cos(angle2))
        mul1 = ielem(n, i)%surf(l)
        b_code = 0
        cleft(1:nof_variables) = u_c(i)%val(1, 1:nof_variables)
        if (ilocal_recon3(i)%mrf .eq. 1) then
          !retrieve rotational velocity in case of rotating reference frame to calculate the correct value of the boundary condition
          srf_speed(2:4) = ilocal_recon3(i)%rotvel(l, 1, 1:3)
          call rotatef(n, srf_speedrot, srf_speed, angle1, angle2)
        end if
        if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

          cturbl(1:turbulenceequations + passivescalar) = u_ct(i)%val(1, 1:turbulenceequations + passivescalar)

        end if

        if (ielem(n, i)%ineighb(l) .eq. n) then        !my cpu only
          if (ielem(n, i)%ibounds(l) .gt. 0) then        !check for boundaries
            if ((ibound(n, ielem(n, i)%ibounds(l))%icode .eq. 5) .or. (ibound(n, ielem(n, i)%ibounds(l))%icode .eq. 50)) then        !periodic in my cpu
              cright(1:nof_variables) = u_c(ielem(n, i)%ineigh(l))%val(1, 1:nof_variables)
              if (per_rot .eq. 1) then
                cright(2:4) = rotate_per_1(cright(2:4), ibound(n, ielem(n, i)%ibounds(l))%icode, angle_per)
              end if
              if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

				cturbr(1:turbulenceequations + passivescalar) = u_ct(ielem(n, i)%ineigh(l))%val(1, 1:turbulenceequations + passivescalar)

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
            cright(1:nof_variables) = u_c(ielem(n, i)%ineigh(l))%val(1, 1:nof_variables)

            if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

				cturbr(1:turbulenceequations + passivescalar) = u_ct(ielem(n, i)%ineigh(l))%val(1, 1:turbulenceequations + passivescalar)

            end if

          end if
        else        !in other cpus they can only be periodic or mpi neighbours

          if (ielem(n, i)%ibounds(l) .gt. 0) then        !check for boundaries
            if ((ibound(n, ielem(n, i)%ibounds(l))%icode .eq. 5) .or. (ibound(n, ielem(n, i)%ibounds(l))%icode .eq. 50)) then        !periodic in other cpu

              if (fastest .eq. 1) then
                cright(1:nof_variables) = solchanger(ielem(n, i)%ineighn(l))%sol(ielem(n, i)%q_face(l)%q_mapl(1), 1:nof_variables)
              else

                cright(1:nof_variables) = iexsolhir(ilocal_recon3(i)%ihexn(1, ielem(n, i)%indexi(l)))%sol &
                                          (ilocal_recon3(i)%ihexl(1, ielem(n, i)%indexi(l)), 1:nof_variables)
              end if
              if (per_rot .eq. 1) then
                cright(2:4) = rotate_per_1(cright(2:4), ibound(n, ielem(n, i)%ibounds(l))%icode, angle_per)
              end if
              if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
                if (fastest .eq. 1) then
                  cturbr(1:turbulenceequations + passivescalar) = solchanger(ielem(n, i)%ineighn(l))%sol &
                                                        (ielem(n, i)%q_face(l)%q_mapl(1), 6:5 + turbulenceequations + passivescalar)
                else

                  cturbr(1:turbulenceequations + passivescalar) = iexsolhir(ilocal_recon3(i)%ihexn(1, ielem(n, i)%indexi(l)))%sol &
                                       (ilocal_recon3(i)%ihexl(1, ielem(n, i)%indexi(l)), 6:5 + turbulenceequations + passivescalar)
                end if
              end if

            end if
          else

            if (fastest .eq. 1) then
              cright(1:nof_variables) = solchanger(ielem(n, i)%ineighn(l))%sol(ielem(n, i)%q_face(l)%q_mapl(1), 1:nof_variables)
            else

              cright(1:nof_variables) = iexsolhir(ilocal_recon3(i)%ihexn(1, ielem(n, i)%indexi(l)))%sol &
                                        (ilocal_recon3(i)%ihexl(1, ielem(n, i)%indexi(l)), 1:nof_variables)
            end if

            if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
              if (fastest .eq. 1) then
                cturbr(1:turbulenceequations + passivescalar) = solchanger(ielem(n, i)%ineighn(l))%sol &
                                                        (ielem(n, i)%q_face(l)%q_mapl(1), 6:5 + turbulenceequations + passivescalar)
              else
                cturbr(1:turbulenceequations + passivescalar) = iexsolhir(ilocal_recon3(i)%ihexn(1, ielem(n, i)%indexi(l)))%sol &
                                       (ilocal_recon3(i)%ihexl(1, ielem(n, i)%indexi(l)), 6:5 + turbulenceequations + passivescalar)
              end if
            end if

!
          end if
        end if

        call rotatef(n, cright_rot, cright, angle1, angle2)        !rotate wrt to normalvector of face and solve 1d riemann problem
        call rotatef(n, cleft_rot, cleft, angle1, angle2)

        if ((lmach .eq. 1)) then    !application of the low mach number correction
          leftv(1:nof_variables) = cleft_rot(1:nof_variables); rightv(1:nof_variables) = cright_rot(1:nof_variables)
          call lmacht(n, leftv, rightv)
          cleft_rot(1:nof_variables) = leftv(1:nof_variables); cright_rot(1:nof_variables) = rightv(1:nof_variables);
          call rotateb(n, cright, cright_rot, angle1, angle2)
          call rotateb(n, cleft, cleft_rot, angle1, angle2)

        end if

        leftv(1:nof_variables) = cleft(1:nof_variables); rightv(1:nof_variables) = cright(1:nof_variables)
        call cons2prim2(n, leftv, rightv, mp_pinfl, mp_pinfr, gammal, gammar)
        asound1 = sqrt(leftv(5)*gamma/leftv(1)) + abs(cleft_rot(2)/cleft_rot(1))
        asound2 = sqrt(rightv(5)*gamma/rightv(1)) + abs(cright_rot(2)/cright_rot(1))
        if (ilocal_recon3(i)%mrf .eq. 1) then
          asound1 = sqrt(leftv(5)*gamma/leftv(1)) + abs(cleft_rot(2)/cleft_rot(1) - srf_speedrot(2))
          asound2 = sqrt(leftv(5)*gamma/leftv(1)) + abs(cright_rot(2)/cright_rot(1) - srf_speedrot(2))
        end if

        vpp = max(asound1, asound2)

        if (itestcase .eq. 4) then
          call sutherland(n, leftv, rightv, viscl, laml)

          viscots = (viscl(1) + viscl(2))*oo2
          mul1 = ielem(n, i)%surf(l)

          viscots = max((4.0/(3.0*0.5*(cleft(1) + cright(1))))*viscots*mul1 &
                        /ielem(n, i)%dih(l), ((gamma/(0.5*(cleft(1) + cright(1)))) &
                                              *viscots*mul1/(ielem(n, i)%dih(l)*prandtl)))
          vpp = max(asound1, asound2) + viscots

          if (turbulence .eq. 1) then
            if (turbulencemodel .eq. 1) then
              turbmv(1) = cturbl(1); turbmv(2) = cturbr(1); eddyfl(2) = turbmv(1); eddyfr(2) = turbmv(2)
              call eddyvisco(n, viscl, laml, turbmv, etvm, eddyfl, eddyfr, leftv, rightv)
            end if
            if (turbulencemodel .eq. 2) then
              eddyfl(1) = ielem(n, i)%walldist; eddyfl(2) = cturbl(1); eddyfl(3) = cturbl(2)
              eddyfl(4:6) = ilocal_recon3(i)%grads(1, 1:3); eddyfl(7:9) = ilocal_recon3(i)%grads(2, 1:3)
              eddyfl(10:12) = ilocal_recon3(i)%grads(3, 1:3); eddyfl(13:15) = ilocal_recon3(i)%grads(5, 1:3)
              eddyfl(16:18) = ilocal_recon3(i)%grads(6, 1:3)

              eddyfr = eddyfl

              call eddyvisco(n, viscl, laml, turbmv, etvm, eddyfl, eddyfr, leftv, rightv)
            end if

            viscots = oo2*((viscl(1) + viscl(3)) + (viscl(2) + viscl(4)))

            viscots = max((4.0/(3.0*0.5*(cleft(1) + cright(1))))*viscots* &
                          mul1/ielem(n, i)%dih(l), ((gamma/(0.5*(cleft(1) + cright(1))))* &
                                                    viscots*mul1/(ielem(n, i)%dih(l)*(prandtl + prtu))))

            vpp = max(asound1, asound2) + viscots
          end if

          impdiag(i, 1:nof_variables, 1:nof_variables) = impdiag(i, 1:nof_variables, 1:nof_variables) + (oo2*((vpp*identity1))*mul1)
          call compute_jacobianse(n, iconsidered, eigvl, cright, gamma, angle1, angle2, srf_speedrot, nx, ny, nz)
          convj = eigvl
                                                  impoff(i,l,1:nof_variables,1:nof_variables)=impoff(i,l,1:nof_variables,1:nof_variables)+(((oo2*convj(1:nof_variables,1:nof_variables))&
                                                                                                    - ((oo2*vpp)*identity1))*mul1)

          if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
            if (turbulence .eq. 1) then
            do nvar = 1, turbulenceequations

              vpp = max(asound1, asound2) + viscots
              impdiagt(i, nvar) = impdiagt(i, nvar) + (oo2*((vpp))*mul1)
              impofft(i, l, nvar) = impofft(i, l, nvar) - (((oo2*vpp))*mul1)
            end do
            end if
            if (passivescalar .gt. 0) then
            do nvar = turbulenceequations + 1, turbulenceequations + passivescalar
              viscl(1) = viscl(1)/schmidt_lam
              viscl(2) = viscl(2)/schmidt_lam

              if (turbulence .eq. 1) then
                viscl(3) = viscl(3)/schmidt_turb
                viscl(4) = viscl(4)/schmidt_turb
                viscots = 0.5*((viscl(1) + viscl(3)) + (viscl(2) + viscl(4)))
              else
                viscots = 0.5*((viscl(1) + viscl(2)))

              end if

              viscots = (2.0*viscots)/((cleft(1) + cright(1))*ielem(n, i)%dih(l))
              vpp = max(asound1, asound2) + viscots
              impdiagt(i, nvar) = impdiagt(i, nvar) + (oo2*((vpp))*mul1)
              impofft(i, l, nvar) = impofft(i, l, nvar) - (oo2*((vpp))*mul1)
            end do
            end if
          end if
        else

          impdiag(i, 1:nof_variables, 1:nof_variables) = impdiag(i, 1:nof_variables, 1:nof_variables) + (oo2*((vpp*identity1))*mul1)
          call compute_jacobianse(n, iconsidered, eigvl, cright, gamma, angle1, angle2, srf_speedrot, nx, ny, nz)
          convj = eigvl
                                                  impoff(i,l,1:nof_variables,1:nof_variables)=impoff(i,l,1:nof_variables,1:nof_variables)+(((oo2*convj(1:nof_variables,1:nof_variables))&
                                                                                                    - ((oo2*vpp)*identity1))*mul1)

        end if

      end do
    end do
    !$omp end do

    !add the contribution of the source term to the jacobian of the diagonal matrix
    if (srfg .eq. 1) then
      !$omp do
      do i = 1, kmaxe
        impdiag(i, 2, 3) = -srf_velocity(3)*ielem(n, i)%totvolume
        impdiag(i, 2, 4) = srf_velocity(2)*ielem(n, i)%totvolume
        impdiag(i, 3, 2) = srf_velocity(3)*ielem(n, i)%totvolume
        impdiag(i, 3, 4) = -srf_velocity(1)*ielem(n, i)%totvolume
        impdiag(i, 4, 2) = -srf_velocity(2)*ielem(n, i)%totvolume
        impdiag(i, 4, 3) = srf_velocity(1)*ielem(n, i)%totvolume
      end do
      !$omp end do
    end if
    if (mrf .eq. 1) then
      !$omp do
      do i = 1, kmaxe
        srf = ilocal_recon3(i)%mrf
        if (ilocal_recon3(i)%mrf .eq. 1) then
          impdiag(i, 2, 3) = -ilocal_recon3(i)%mrf_velocity(3)*ielem(n, i)%totvolume
          impdiag(i, 2, 4) = ilocal_recon3(i)%mrf_velocity(2)*ielem(n, i)%totvolume
          impdiag(i, 3, 2) = ilocal_recon3(i)%mrf_velocity(3)*ielem(n, i)%totvolume
          impdiag(i, 3, 4) = -ilocal_recon3(i)%mrf_velocity(1)*ielem(n, i)%totvolume
          impdiag(i, 4, 2) = -ilocal_recon3(i)%mrf_velocity(2)*ielem(n, i)%totvolume
          impdiag(i, 4, 3) = ilocal_recon3(i)%mrf_velocity(1)*ielem(n, i)%totvolume
        end if
        srf = 0
      end do
      !$omp end do
    end if

    if (rungekutta .eq. 10) then
      !$omp do
      do i = 1, kmaxe
        do l = 1, nof_variables
          impdiag(i, l, l) = impdiag(i, l, l) + (ielem(n, i)%totvolume/ielem(n, i)%dtl)
        end do
      end do
      !$omp end do
    else
      !$omp do
!
      do i = 1, kmaxe
        do l = 1, nof_variables
          impdiag(i, l, l) = ielem(n, i)%totvolume*((1.0d0/ielem(n, i)%dtl) + (1.5d0/dt)) + (impdiag(i, l, l))
        end do
      end do
      !$omp end do
    end if

    if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
      if (turbulence .eq. 1) call sources_derivatives_computation(n)
      if (rungekutta .eq. 10) then
!$omp do
        do i = 1, kmaxe
          if (turbulence .eq. 1) then
          do nvar = 1, turbulenceequations
            impdiagt(i, nvar) = impdiagt(i, nvar) + (ielem(n, i)%totvolume/(ielem(n, i)%dtl)) - sht(i, nvar)
          end do
          end if
          if (passivescalar .gt. 0) then
          do nvar = turbulenceequations + 1, turbulenceequations + passivescalar
            impdiagt(i, nvar) = impdiagt(i, nvar) + (ielem(n, i)%totvolume/ielem(n, i)%dtl)
          end do
          end if
        end do
!$omp end do
      else
!$omp do
        do i = 1, kmaxe
          if (turbulence .eq. 1) then
          do nvar = 1, turbulenceequations

            impdiagt(i, nvar) = ielem(n, i)%totvolume*((1.0d0/ielem(n, i)%dtl) + (1.5d0/dt)) + (impdiagt(i, nvar)) - sht(i, nvar)
          end do
          end if
          if (passivescalar .gt. 0) then
          do nvar = turbulenceequations + 1, turbulenceequations + passivescalar
            impdiagt(i, nvar) = ielem(n, i)%totvolume*((1.0d0/ielem(n, i)%dtl) + (1.5d0/dt)) + (impdiagt(i, 1))
          end do
          end if
        end do
!$omp end do

      end if
    end if

  end subroutine calculate_jacobian

  subroutine calculate_jacobian_2d(n)
    !> @brief
!> this subroutine computes the approximate jacobian for implicit time stepping in 2d
    implicit none
    integer, intent(in)::n
    real, dimension(1:nof_variables + turbulenceequations + passivescalar)::godflux2
    integer::i, l, ngp, kmaxe, iqp, ii, nvar, n_node, ibfc
    real::sum_detect, norms, vpp, asound1, asound2, mul1, dxb, tempxx, viscots
    real, dimension(nof_variables, nof_variables)::identity1
    real, dimension(nof_variables, nof_variables)::convj, diffj
    integer::iconsidered, facex, pointx, igoflux, kas
    integer::b_code
    real::angle1, angle2, nx, ny, nz
    real, dimension(1:nof_variables + turbulenceequations + passivescalar)::cleft, cright, cright_rot, cleft_rot
    real, dimension(1:turbulenceequations + passivescalar)::cturbl, cturbr
    real, dimension(1:nof_variables)::leftv, srf_speedrot, srf_speed
    real, dimension(1:nof_variables)::rightv
    real, dimension(1:dimensiona)::pox, poy, poz
    real, dimension(1:8, 1:dimensiona)::vext, nodes_list
    real, dimension(1:dimensiona)::cords
    real, dimension(1:4)::viscl, laml
    real, dimension(1:2)::turbmv
    real, dimension(1)::etvm
    real, dimension(1:20)::eddyfl, eddyfr

    real::mp_pinfl, gammal
    real::mp_pinfr, gammar
    real, dimension(1:nof_variables, 1:nof_variables)::eigvl

    kmaxe = xmpielrank(n)
    identity1(:, :) = zero
    identity1(1, 1) = 1.0d0
    identity1(2, 2) = 1.0d0
    identity1(3, 3) = 1.0d0
    identity1(4, 4) = 1.0d0

    !$omp do
    do ii = 1, nof_interior        !for all the interior elements
      i = el_int(ii)
      iconsidered = i
      impdiag(i, :, :) = zero
      impoff(i, :, :, :) = zero
      if (turbulence .eq. 1) then
        impdiagt(i, :) = zero
        impofft(i, :, :) = zero
      end if

      do l = 1, ielem(n, i)%ifca !for all their faces
        b_code = 0
        angle1 = ielem(n, i)%faceanglex(l)
        angle2 = ielem(n, i)%faceangley(l)
        nx = angle1
        ny = angle2
        mul1 = ielem(n, i)%surf(l)

        cleft(1:nof_variables) = u_c(i)%val(1, 1:nof_variables)
        cright(1:nof_variables) = u_c(ielem(n, i)%ineigh(l))%val(1, 1:nof_variables)

        if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

          cturbl(1:turbulenceequations + passivescalar) = u_ct(i)%val(1, 1:turbulenceequations + passivescalar)!left additional equations flow state
          cturbr(1:turbulenceequations + passivescalar) = u_ct(ielem(n, i)%ineigh(l))%val(1, 1:turbulenceequations + passivescalar)

        end if

        call rotatef2d(n, cright_rot, cright, angle1, angle2)        !rotate wrt to normalvector of face and solve 1d riemann problem
        call rotatef2d(n, cleft_rot, cleft, angle1, angle2)        !rotate wrt to normalvector of face and solve 1d riemann problem

        if ((lmach .eq. 1)) then    !application of the low mach number correction
          leftv(1:nof_variables) = cleft_rot(1:nof_variables); rightv(1:nof_variables) = cright_rot(1:nof_variables)
          call lmacht2d(n, leftv, rightv)
          cleft_rot(1:nof_variables) = leftv(1:nof_variables); cright_rot(1:nof_variables) = rightv(1:nof_variables);
          call rotateb2d(n, cright, cright_rot, angle1, angle2)        !rotate wrt to normalvector of face and solve 1d riemann problem
          call rotateb2d(n, cleft, cleft_rot, angle1, angle2)

        end if

        leftv(1:nof_variables) = cleft(1:nof_variables); rightv(1:nof_variables) = cright(1:nof_variables)
        call cons2prim2(n, leftv, rightv, mp_pinfl, mp_pinfr, gammal, gammar)

        asound1 = sqrt(leftv(4)*gamma/leftv(1)) + abs(cleft_rot(2)/cleft_rot(1))
        asound2 = sqrt(rightv(4)*gamma/rightv(1)) + abs(cright_rot(2)/cright_rot(1))

        vpp = max(asound1, asound2)

        if (itestcase .eq. 4) then
          call sutherland2d(n, leftv, rightv, viscl, laml)

          viscots = (viscl(1) + viscl(2))*oo2
          mul1 = ielem(n, i)%surf(l)
!                                                   viscots=viscots/((0.5*(cleft(1)+cright(1)))*ielem(n,i)%dih(l))
          viscots = max((4.0/(3.0*0.5*(cleft(1) + cright(1))))*viscots*mul1 &
                        /ielem(n, i)%dih(l), ((gamma/(0.5*(cleft(1) + cright(1)))) &
                                              *viscots*mul1/(ielem(n, i)%dih(l)*prandtl)))
          vpp = max(asound1, asound2) + viscots

          if (turbulence .eq. 1) then
            if (turbulencemodel .eq. 1) then
              turbmv(1) = cturbl(1); turbmv(2) = cturbr(1); eddyfl(2) = turbmv(1); eddyfr(2) = turbmv(2)
              call eddyvisco2d(n, viscl, laml, turbmv, etvm, eddyfl, eddyfr, leftv, rightv)
            end if
            if (turbulencemodel .eq. 2) then
              eddyfl(1) = ielem(n, i)%walldist; eddyfl(2) = cturbl(1); eddyfl(3) = cturbl(2)
              eddyfl(4:5) = ilocal_recon3(i)%grads(1, 1:2); eddyfl(6:7) = ilocal_recon3(i)%grads(2, 1:2)
              eddyfl(8:9) = ilocal_recon3(i)%grads(4, 1:2)
              eddyfl(10:11) = ilocal_recon3(i)%grads(5, 1:2)

              eddyfr = eddyfl
              call eddyvisco2d(n, viscl, laml, turbmv, etvm, eddyfl, eddyfr, leftv, rightv)
            end if

            viscots = oo2*((viscl(1) + viscl(3)) + (viscl(2) + viscl(4)))

!                                                       viscots=viscots/((0.5*(cleft(1)+cright(1)))*ielem(n,i)%dih(l))
            viscots = max((4.0/(3.0*0.5*(cleft(1) + cright(1))))*viscots* &
                          mul1/ielem(n, i)%dih(l), ((gamma/(0.5*(cleft(1) + cright(1))))* &
                                                    viscots*mul1/(ielem(n, i)%dih(l)*(prandtl + prtu))))

            vpp = max(asound1, asound2) + viscots
          end if

          impdiag(i, 1:nof_variables, 1:nof_variables) = impdiag(i, 1:nof_variables, 1:nof_variables) + (oo2*((vpp*identity1))*mul1)
          call compute_jacobianse2d(n, eigvl, cright, gamma, angle1, angle2, nx, ny, nz)
          convj = eigvl
          impoff(i,l,1:nof_variables,1:nof_variables)=impoff(i,l,1:nof_variables,1:nof_variables)+(((oo2*convj(1:nof_variables,1:nof_variables))&
                - ((oo2*vpp)*identity1))*mul1)

          if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
            if (turbulence .eq. 1) then
            do nvar = 1, turbulenceequations

              vpp = max(asound1, asound2) + viscots

              impdiagt(i, nvar) = impdiagt(i, nvar) + (oo2*((vpp))*mul1)
              impofft(i, l, nvar) = impofft(i, l, nvar) - (((oo2*vpp))*mul1)
            end do
            end if
            if (passivescalar .gt. 0) then

              do nvar = turbulenceequations + 1, turbulenceequations + passivescalar
                viscl(1) = viscl(1)/schmidt_lam
                viscl(2) = viscl(2)/schmidt_lam

                if (turbulence .eq. 1) then
                  viscl(3) = viscl(3)/schmidt_turb
                  viscl(4) = viscl(4)/schmidt_turb
                  viscots = 0.5*((viscl(1) + viscl(3)) + (viscl(2) + viscl(4)))
                else
                  viscots = 0.5*((viscl(1) + viscl(2)))

                end if

                viscots = (2.0*viscots)/((cleft(1) + cright(1))*ielem(n, i)%dih(l))
                vpp = max(asound1, asound2) + viscots
                impdiagt(i, nvar) = impdiagt(i, nvar) + (oo2*((vpp))*mul1)
                impofft(i, l, nvar) = impofft(i, l, nvar) - (oo2*((vpp))*mul1)
              end do
            end if
          end if
        else

          impdiag(i, 1:nof_variables, 1:nof_variables) = impdiag(i, 1:nof_variables, 1:nof_variables) + (oo2*((vpp*identity1))*mul1)
          call compute_jacobianse2d(n, eigvl, cright, gamma, angle1, angle2, nx, ny, nz)
          convj = eigvl
                                                  impoff(i,l,1:nof_variables,1:nof_variables)=impoff(i,l,1:nof_variables,1:nof_variables)+(((oo2*convj(1:nof_variables,1:nof_variables))&
                                                                                                    - ((oo2*vpp)*identity1))*mul1)

        end if
      end do
    end do
    !$omp end do

    !$omp do
    do ii = 1, nof_bounded
      i = el_bnd(ii)
      iconsidered = i

      impdiag(i, :, :) = zero
      impoff(i, :, :, :) = zero
      if (turbulence .eq. 1) then
        impdiagt(i, :) = zero
        impofft(i, :, :) = zero
      end if

      do l = 1, ielem(n, i)%ifca
        mul1 = ielem(n, i)%surf(l)
        angle1 = ielem(n, i)%faceanglex(l)
        angle2 = ielem(n, i)%faceangley(l)
        nx = angle1
        ny = angle2

        b_code = 0
        cleft(1:nof_variables) = u_c(i)%val(1, 1:nof_variables)
        if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

          cturbl(1:turbulenceequations + passivescalar) = u_ct(i)%val(1, 1:turbulenceequations + passivescalar)

        end if

        if (ielem(n, i)%ineighb(l) .eq. n) then        !my cpu only
          if (ielem(n, i)%ibounds(l) .gt. 0) then        !check for boundaries
            if (ibound(n, ielem(n, i)%ibounds(l))%icode .eq. 5) then        !periodic in my cpu
              cright(1:nof_variables) = u_c(ielem(n, i)%ineigh(l))%val(1, 1:nof_variables)

              if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

           cturbr(1:turbulenceequations + passivescalar) = u_ct(ielem(n, i)%ineigh(l))%val(1, 1:turbulenceequations + passivescalar)

              end if

              kas = 1

            else
              !not periodic ones in my cpu

              facex = l; iconsidered = i
              call coordinates_face_inner2dx(n, iconsidered, facex, vext, nodes_list)
              n_node = 2
              cords(1:2) = zero
              cords(1:2) = cordinates2(n, nodes_list, n_node)

              poy(1) = cords(2)
              pox(1) = cords(1)

              leftv(1:nof_variables) = cleft(1:nof_variables)
              b_code = ibound(n, ielem(n, i)%ibounds(l))%icode

              call boundarys2d(n,b_code,iconsidered,facex,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speed,srf_speedrot,ibfc)
              cright(1:nof_variables) = rightv(1:nof_variables)

              kas = 2

            end if
          else
            cright(1:nof_variables) = u_c(ielem(n, i)%ineigh(l))%val(1, 1:nof_variables)

            if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

           cturbr(1:turbulenceequations + passivescalar) = u_ct(ielem(n, i)%ineigh(l))%val(1, 1:turbulenceequations + passivescalar)

            end if

            kas = 3

          end if
        else        !in other cpus they can only be periodic or mpi neighbours

          if (ielem(n, i)%ibounds(l) .gt. 0) then        !check for boundaries
            if (ibound(n, ielem(n, i)%ibounds(l))%icode .eq. 5) then        !periodic in other cpu

              if (fastest .eq. 1) then
                cright(1:nof_variables) = solchanger(ielem(n, i)%ineighn(l))%sol(ielem(n, i)%q_face(l)%q_mapl(1), 1:nof_variables)
              else

                cright(1:nof_variables) = iexsolhir(ilocal_recon3(i)%ihexn(1, ielem(n, i)%indexi(l)))%sol &
                                          (ilocal_recon3(i)%ihexl(1, ielem(n, i)%indexi(l)), 1:nof_variables)
              end if
              kas = 4

              if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
                if (fastest .eq. 1) then
                  cturbr(1:turbulenceequations + passivescalar) = solchanger(ielem(n, i)%ineighn(l))%sol &
                                                        (ielem(n, i)%q_face(l)%q_mapl(1), 5:4 + turbulenceequations + passivescalar)
                else

                  cturbr(1:turbulenceequations + passivescalar) = iexsolhir(ilocal_recon3(i)%ihexn(1, ielem(n, i)%indexi(l)))%sol &
                                       (ilocal_recon3(i)%ihexl(1, ielem(n, i)%indexi(l)), 5:4 + turbulenceequations + passivescalar)
                end if
              end if

            end if
          else
            kas = 5
            if (fastest .eq. 1) then
              cright(1:nof_variables) = solchanger(ielem(n, i)%ineighn(l))%sol(ielem(n, i)%q_face(l)%q_mapl(1), 1:nof_variables)
            else

              cright(1:nof_variables) = iexsolhir(ilocal_recon3(i)%ihexn(1, ielem(n, i)%indexi(l)))%sol &
                                        (ilocal_recon3(i)%ihexl(1, ielem(n, i)%indexi(l)), 1:nof_variables)
            end if

            if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
              if (fastest .eq. 1) then
                cturbr(1:turbulenceequations + passivescalar) = solchanger(ielem(n, i)%ineighn(l))%sol &
                    (ielem(n, i)%q_face(l)%q_mapl(1), 5:4 + turbulenceequations + passivescalar)
              else

                cturbr(1:turbulenceequations + passivescalar) = iexsolhir(ilocal_recon3(i)%ihexn(1, ielem(n, i)%indexi(l)))%sol &
                                       (ilocal_recon3(i)%ihexl(1, ielem(n, i)%indexi(l)), 5:4 + turbulenceequations + passivescalar)
              end if
            end if

!
          end if
        end if

        call rotatef2d(n, cright_rot, cright, angle1, angle2)        !rotate wrt to normalvector of face and solve 1d riemann problem
        call rotatef2d(n, cleft_rot, cleft, angle1, angle2)        !rotate wrt to normalvector of face and solve 1d riemann problem

        if ((lmach .eq. 1)) then    !application of the low mach number correction
          leftv(1:nof_variables) = cleft_rot(1:nof_variables); rightv(1:nof_variables) = cright_rot(1:nof_variables)
          call lmacht2d(n, leftv, rightv)
          cleft_rot(1:nof_variables) = leftv(1:nof_variables); cright_rot(1:nof_variables) = rightv(1:nof_variables);
          call rotateb2d(n, cright, cright_rot, angle1, angle2)        !rotate wrt to normalvector of face and solve 1d riemann problem
          call rotateb2d(n, cleft, cleft_rot, angle1, angle2)

        end if

        leftv(1:nof_variables) = cleft(1:nof_variables); rightv(1:nof_variables) = cright(1:nof_variables)
        call cons2prim2(n, leftv, rightv, mp_pinfl, mp_pinfr, gammal, gammar)

        asound1 = sqrt(leftv(4)*gamma/leftv(1)) + abs(cleft_rot(2)/cleft_rot(1))
        asound2 = sqrt(rightv(4)*gamma/rightv(1)) + abs(cright_rot(2)/cright_rot(1))

        vpp = max(asound1, asound2)

        if (itestcase .eq. 4) then
          call sutherland2d(n, leftv, rightv, viscl, laml)

          viscots = (viscl(1) + viscl(2))*oo2
          mul1 = ielem(n, i)%surf(l)
!                                                   viscots=viscots/((0.5*(cleft(1)+cright(1)))*ielem(n,i)%dih(l))
          viscots = max((4.0/(3.0*0.5*(cleft(1) + cright(1))))*viscots*mul1 &
                        /ielem(n, i)%dih(l), ((gamma/(0.5*(cleft(1) + cright(1)))) &
                                              *viscots*mul1/(ielem(n, i)%dih(l)*prandtl)))
          vpp = max(asound1, asound2) + viscots

          if (turbulence .eq. 1) then
            if (turbulencemodel .eq. 1) then
              turbmv(1) = cturbl(1); turbmv(2) = cturbr(1); eddyfl(2) = turbmv(1); eddyfr(2) = turbmv(2)
              call eddyvisco2d(n, viscl, laml, turbmv, etvm, eddyfl, eddyfr, leftv, rightv)
            end if
            if (turbulencemodel .eq. 2) then
              eddyfl(1) = ielem(n, i)%walldist; eddyfl(2) = cturbl(1); eddyfl(3) = cturbl(2)
              eddyfl(4:5) = ilocal_recon3(i)%grads(1, 1:2); eddyfl(6:7) = ilocal_recon3(i)%grads(2, 1:2)
              eddyfl(8:9) = ilocal_recon3(i)%grads(4, 1:2)
              eddyfl(10:11) = ilocal_recon3(i)%grads(5, 1:2)

              eddyfr = eddyfl
              call eddyvisco2d(n, viscl, laml, turbmv, etvm, eddyfl, eddyfr, leftv, rightv)
            end if

            viscots = oo2*((viscl(1) + viscl(3)) + (viscl(2) + viscl(4)))
!                                                       viscots=viscots/((0.5*(cleft(1)+cright(1)))*ielem(n,i)%dih(l))
            viscots = max((4.0/(3.0*0.5*(cleft(1) + cright(1))))*viscots* &
                          mul1/ielem(n, i)%dih(l), ((gamma/(0.5*(cleft(1) + cright(1))))* &
                                                    viscots*mul1/(ielem(n, i)%dih(l)*(prandtl + prtu))))

            vpp = max(asound1, asound2) + viscots
          end if

          impdiag(i, 1:nof_variables, 1:nof_variables) = impdiag(i, 1:nof_variables, 1:nof_variables) + (oo2*((vpp*identity1))*mul1)
          call compute_jacobianse2d(n, eigvl, cright, gamma, angle1, angle2, nx, ny, nz)
          convj = eigvl
                                                  impoff(i,l,1:nof_variables,1:nof_variables)=impoff(i,l,1:nof_variables,1:nof_variables)+(((oo2*convj(1:nof_variables,1:nof_variables))&
                                                                                                    - ((oo2*vpp)*identity1))*mul1)

          if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
            if (turbulence .eq. 1) then
            do nvar = 1, turbulenceequations

              vpp = max(asound1, asound2) + viscots
!
              impdiagt(i, nvar) = impdiagt(i, nvar) + (oo2*((vpp))*mul1)
              impofft(i, l, nvar) = impofft(i, l, nvar) - (((oo2*vpp))*mul1)
            end do
            end if
            if (passivescalar .gt. 0) then
            do nvar = turbulenceequations + 1, turbulenceequations + passivescalar
              viscl(1) = viscl(1)/schmidt_lam
              viscl(2) = viscl(2)/schmidt_lam

              if (turbulence .eq. 1) then
                viscl(3) = viscl(3)/schmidt_turb
                viscl(4) = viscl(4)/schmidt_turb
                viscots = 0.5*((viscl(1) + viscl(3)) + (viscl(2) + viscl(4)))
              else
                viscots = 0.5*((viscl(1) + viscl(2)))

              end if

              viscots = (2.0*viscots)/((cleft(1) + cright(1))*ielem(n, i)%dih(l))
              vpp = max(asound1, asound2) + viscots
              impdiagt(i, nvar) = impdiagt(i, nvar) + (oo2*((vpp))*mul1)
              impofft(i, l, nvar) = impofft(i, l, nvar) - (oo2*((vpp))*mul1)
            end do
            end if
          end if
        else

          impdiag(i, 1:nof_variables, 1:nof_variables) = impdiag(i, 1:nof_variables, 1:nof_variables) + (oo2*((vpp*identity1))*mul1)
          call compute_jacobianse2d(n, eigvl, cright, gamma, angle1, angle2, nx, ny, nz)
          convj = eigvl
                                                  impoff(i,l,1:nof_variables,1:nof_variables)=impoff(i,l,1:nof_variables,1:nof_variables)+(((oo2*convj(1:nof_variables,1:nof_variables))&
                                                                                                    - ((oo2*vpp)*identity1))*mul1)

        end if

      end do
    end do
    !$omp end do

    if (rungekutta .eq. 10) then
      !$omp do
      do i = 1, kmaxe

        impdiag(i, 1, 1) = impdiag(i, 1, 1) + (ielem(n, i)%totvolume/ielem(n, i)%dtl)
        impdiag(i, 2, 2) = impdiag(i, 2, 2) + (ielem(n, i)%totvolume/ielem(n, i)%dtl)
        impdiag(i, 3, 3) = impdiag(i, 3, 3) + (ielem(n, i)%totvolume/ielem(n, i)%dtl)
        impdiag(i, 4, 4) = impdiag(i, 4, 4) + (ielem(n, i)%totvolume/ielem(n, i)%dtl)

      end do
      !$omp end do
    else
      !$omp do
      do i = 1, kmaxe
        impdiag(i, 1, 1) = ielem(n, i)%totvolume*((1.0d0/ielem(n, i)%dtl) + (1.5d0/dt)) + (impdiag(i, 1, 1))
        impdiag(i, 2, 2) = ielem(n, i)%totvolume*((1.0d0/ielem(n, i)%dtl) + (1.5d0/dt)) + (impdiag(i, 2, 2))
        impdiag(i, 3, 3) = ielem(n, i)%totvolume*((1.0d0/ielem(n, i)%dtl) + (1.5d0/dt)) + (impdiag(i, 3, 3))
        impdiag(i, 4, 4) = ielem(n, i)%totvolume*((1.0d0/ielem(n, i)%dtl) + (1.5d0/dt)) + (impdiag(i, 4, 4))

      end do
      !$omp end do
    end if

    if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then

      if (turbulence .eq. 1) call sources_derivatives_computation2d(n)
      if (rungekutta .eq. 10) then
!$omp do
        do i = 1, kmaxe
          if (turbulence .eq. 1) then
          do nvar = 1, turbulenceequations
!
            impdiagt(i, nvar) = impdiagt(i, nvar) + ((ielem(n, i)%totvolume/(ielem(n, i)%dtl))) - sht(i, nvar)
          end do
          end if
          if (passivescalar .gt. 0) then
          do nvar = turbulenceequations + 1, turbulenceequations + passivescalar
            impdiagt(i, nvar) = impdiagt(i, nvar) + (ielem(n, i)%totvolume/(ielem(n, i)%dtl))
          end do
          end if
        end do
!$omp end do
      else
!$omp do
        do i = 1, kmaxe
          if (turbulence .eq. 1) then
          do nvar = 1, turbulenceequations
            impdiagt(i, nvar) = ielem(n, i)%totvolume*((1.0d0/(ielem(n, i)%dtl)) + (1.5d0/dt)) + (impdiagt(i, 1)) - sht(i, nvar)

          end do
          end if
          if (passivescalar .gt. 0) then
          do nvar = turbulenceequations + 1, turbulenceequations + passivescalar
            impdiagt(i, nvar) = ielem(n, i)%totvolume*((1.0d0/(ielem(n, i)%dtl)) + (1.5d0/dt)) + (impdiagt(i, 1))

          end do
          end if
        end do
!$omp end do

      end if
    end if

  end subroutine calculate_jacobian_2d

  subroutine calculate_jacobianlm(n, iconsidered, impdiag, impdiagt, impoff, impofft)
    implicit none
    integer, intent(in)::n, iconsidered
    real, dimension(1:nof_variables + turbulenceequations + passivescalar)::godflux2
    integer::i, l, ngp, kmaxe, iqp, ii, nvar, n_node, ibfc
    real::sum_detect, norms, vpp, asound1, asound2, mul1, dxb, tempxx, viscots
    real, dimension(nof_variables, nof_variables)::identity1
    real, dimension(nof_variables, nof_variables)::convj, diffj
    integer::facex, pointx, igoflux
    integer::b_code
    real::angle1, angle2, nx, ny, nz
    real, dimension(1:nof_variables)::cleft, cright, cright_rot, cleft_rot
    real, dimension(1:turbulenceequations + passivescalar)::cturbl, cturbr
    real, dimension(1:nof_variables)::leftv, srf_speedrot, srf_speed
    real, dimension(1:nof_variables)::rightv
    real, dimension(1:dimensiona)::pox, poy, poz

    real, dimension(1:8, 1:dimensiona)::vext, nodes_list
    real, dimension(1:4)::viscl, laml
    real, dimension(1:2)::turbmv
    real, dimension(1)::etvm
    real, dimension(1:20)::eddyfl, eddyfr

    real, dimension(1:dimensiona)::cords
    real::mp_pinfl, gammal
    real::mp_pinfr, gammar
    real, allocatable, dimension(:, :), intent(inout)::impdiagt
    real, allocatable, dimension(:, :, :), intent(inout)::impdiag, impofft
    real, allocatable, dimension(:, :, :, :), intent(inout)::impoff
    real, dimension(1:nof_variables, 1:nof_variables)::eigvl
    real, dimension(turbulenceequations)::source_t

    identity1(:, :) = zero
    identity1(1, 1) = 1.0d0
    identity1(2, 2) = 1.0d0
    identity1(3, 3) = 1.0d0
    identity1(4, 4) = 1.0d0
    identity1(5, 5) = 1.0d0

    if (ielem(n, iconsidered)%interior .eq. 0) then

      i = iconsidered
      impdiag(1, :, :) = zero
      impoff(1, :, :, :) = zero
      if (turbulence .eq. 1) then
        impdiagt(1, :) = zero
        impofft(1, :, :) = zero
      end if

      do l = 1, ielem(n, i)%ifca !for all their faces
        godflux2 = zero
        angle1 = ielem(n, i)%faceanglex(l)
        angle2 = ielem(n, i)%faceangley(l)
        nx = (cos(angle1)*sin(angle2))
        ny = (sin(angle1)*sin(angle2))
        nz = (cos(angle2))
        mul1 = ielem(n, i)%surf(l)

        cleft(1:nof_variables) = u_c(i)%val(1, 1:nof_variables)
        cright(1:nof_variables) = u_c(ielem(n, i)%ineigh(l))%val(1, 1:nof_variables)

        if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

          cturbl(1:turbulenceequations + passivescalar) = u_ct(i)%val(1, 1:turbulenceequations + passivescalar)!left additional equations flow state
          cturbr(1:turbulenceequations + passivescalar) = u_ct(ielem(n, i)%ineigh(l))%val(1, 1:turbulenceequations + passivescalar)

        end if

        call rotatef(n, cright_rot, cright, angle1, angle2)        !rotate wrt to normalvector of face and solve 1d riemann problem
        call rotatef(n, cleft_rot, cleft, angle1, angle2)        !rotate wrt to normalvector of face and solve 1d riemann problem

        if ((lmach .eq. 1)) then    !application of the low mach number correction
          leftv(1:nof_variables) = cleft_rot(1:nof_variables); rightv(1:nof_variables) = cright_rot(1:nof_variables)
          call lmacht(n, leftv, rightv)
          cleft_rot(1:nof_variables) = leftv(1:nof_variables); cright_rot(1:nof_variables) = rightv(1:nof_variables);
          call rotatef(n, cright, cright_rot, angle1, angle2)        !rotate wrt to normalvector of face and solve 1d riemann problem
          call rotatef(n, cleft, cleft_rot, angle1, angle2)

        end if

        leftv(1:nof_variables) = cleft(1:nof_variables); rightv(1:nof_variables) = cright(1:nof_variables)
        call cons2prim2(n, leftv, rightv, mp_pinfl, mp_pinfr, gammal, gammar)

        asound1 = sqrt(leftv(5)*gamma/leftv(1)) + abs(cleft_rot(2)/cleft_rot(1))
        asound2 = sqrt(rightv(5)*gamma/rightv(1)) + abs(cright_rot(2)/cright_rot(1))

        vpp = max(asound1, asound2)

        if (itestcase .eq. 4) then
          call sutherland(n, leftv, rightv, viscl, laml)

          viscots = (viscl(1) + viscl(2))*oo2
          mul1 = ielem(n, i)%surf(l)

          viscots = max((4.0/(3.0*0.5*(cleft(1) + cright(1))))*viscots*mul1 &
                        /ielem(n, i)%dih(l), ((gamma/(0.5*(cleft(1) + cright(1)))) &
                                              *viscots*mul1/(ielem(n, i)%dih(l)*prandtl)))
          vpp = max(asound1, asound2) + viscots

          if (turbulence .eq. 1) then
            if (turbulencemodel .eq. 1) then
              turbmv(1) = cturbl(1); turbmv(2) = cturbr(1); eddyfl(2) = turbmv(1); eddyfr(2) = turbmv(2)
              call eddyvisco(n, viscl, laml, turbmv, etvm, eddyfl, eddyfr, leftv, rightv)
            end if
            if (turbulencemodel .eq. 2) then

              eddyfl(1) = ielem(n, i)%walldist; eddyfl(2) = cturbl(1); eddyfl(3) = cturbl(2)

              eddyfl(4:6) = ilocal_recon3(i)%grads(1, 1:3); eddyfl(7:9) = ilocal_recon3(i)%grads(2, 1:3)
              eddyfl(10:12) = ilocal_recon3(i)%grads(3, 1:3); eddyfl(13:15) = ilocal_recon3(i)%grads(5, 1:3)
              eddyfl(16:18) = ilocal_recon3(i)%grads(6, 1:3)

              eddyfr = eddyfl
              call eddyvisco(n, viscl, laml, turbmv, etvm, eddyfl, eddyfr, leftv, rightv)
            end if

            viscots = oo2*((viscl(1) + viscl(3)) + (viscl(2) + viscl(4)))

            viscots = max((4.0/(3.0*0.5*(cleft(1) + cright(1))))*viscots* &
                          mul1/ielem(n, i)%dih(l), ((gamma/(0.5*(cleft(1) + cright(1))))* &
                                                    viscots*mul1/(ielem(n, i)%dih(l)*(prandtl + prtu))))

            vpp = max(asound1, asound2) + viscots
          end if

          impdiag(1, 1:nof_variables, 1:nof_variables) = impdiag(1, 1:nof_variables, 1:nof_variables) + (oo2*((vpp*identity1))*mul1)
          call compute_jacobianse(n, iconsidered, eigvl, cright, gamma, angle1, angle2, srf_speedrot, nx, ny, nz)
          convj = eigvl
                                                  impoff(1,l,1:nof_variables,1:nof_variables)=impoff(1,l,1:nof_variables,1:nof_variables)+(((oo2*convj(1:nof_variables,1:nof_variables))&
                                                                                                    - ((oo2*vpp)*identity1))*mul1)

          if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
            if (turbulence .eq. 1) then
            do nvar = 1, turbulenceequations

              vpp = max(asound1, asound2) + viscots
              impdiagt(1, nvar) = impdiagt(1, nvar) + (oo2*((vpp))*mul1)
              impofft(1, l, nvar) = impofft(1, l, nvar) - (((oo2*vpp))*mul1)
            end do
            end if
            if (passivescalar .gt. 0) then
            do nvar = turbulenceequations + 1, turbulenceequations + passivescalar
              viscl(1) = viscl(1)/schmidt_lam
              viscl(2) = viscl(2)/schmidt_lam
              viscl(3) = viscl(3)/schmidt_turb
              viscl(4) = viscl(4)/schmidt_turb
              viscots = 0.5*((viscl(1) + viscl(3)) + (viscl(2) + viscl(4)))

              viscots = (2.0*viscots)/((cleft(1) + cright(1))*ielem(n, i)%dih(l))
              vpp = max(asound1, asound2) + viscots
              impdiagt(1, nvar) = impdiagt(1, nvar) + (oo2*((vpp))*mul1)
              impofft(1, l, nvar) = impofft(1, l, nvar) - (oo2*((vpp))*mul1)
            end do
            end if
          end if
        else

          impdiag(1, 1:nof_variables, 1:nof_variables) = impdiag(1, 1:nof_variables, 1:nof_variables) + (oo2*((vpp*identity1))*mul1)
          call compute_jacobianse(n, iconsidered, eigvl, cright, gamma, angle1, angle2, srf_speedrot, nx, ny, nz)
          convj = eigvl
                                                  impoff(1,l,1:nof_variables,1:nof_variables)=impoff(1,l,1:nof_variables,1:nof_variables)+(((oo2*convj(1:nof_variables,1:nof_variables))&
                                                                                                    - ((oo2*vpp)*identity1))*mul1)

        end if
      end do

    else
      i = iconsidered
      impdiag(1, :, :) = zero
      impoff(1, :, :, :) = zero
      if (turbulence .eq. 1) then
        impdiagt(1, :) = zero
        impofft(1, :, :) = zero
      end if

      do l = 1, ielem(n, i)%ifca
        mul1 = ielem(n, i)%surf(l)
        angle1 = ielem(n, i)%faceanglex(l)
        angle2 = ielem(n, i)%faceangley(l)
        nx = (cos(angle1)*sin(angle2))
        ny = (sin(angle1)*sin(angle2))
        nz = (cos(angle2))

        b_code = 0
        cleft(1:nof_variables) = u_c(i)%val(1, 1:nof_variables)
        if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

          cturbl(1:turbulenceequations + passivescalar) = u_ct(i)%val(1, 1:turbulenceequations + passivescalar)

        end if

        if (ielem(n, i)%ineighb(l) .eq. n) then        !my cpu only
          if (ielem(n, i)%ibounds(l) .gt. 0) then        !check for boundaries
            if (ibound(n, ielem(n, i)%ibounds(l))%icode .eq. 5) then        !periodic in my cpu
              cright(1:nof_variables) = u_c(ielem(n, i)%ineigh(l))%val(1, 1:nof_variables)

              if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

           cturbr(1:turbulenceequations + passivescalar) = u_ct(ielem(n, i)%ineigh(l))%val(1, 1:turbulenceequations + passivescalar)

              end if

            else

              facex = l;
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
            cright(1:nof_variables) = u_c(ielem(n, i)%ineigh(l))%val(1, 1:nof_variables)

            if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

           cturbr(1:turbulenceequations + passivescalar) = u_ct(ielem(n, i)%ineigh(l))%val(1, 1:turbulenceequations + passivescalar)

            end if

          end if
        else        !in other cpus they can only be periodic or mpi neighbours

          if (ielem(n, i)%ibounds(l) .gt. 0) then        !check for boundaries
            if (ibound(n, ielem(n, i)%ibounds(l))%icode .eq. 5) then        !periodic in other cpu

              if (fastest .eq. 1) then
                cright(1:nof_variables) = solchanger(ielem(n, i)%ineighn(l))%sol(ielem(n, i)%q_face(l)%q_mapl(1), 1:nof_variables)
              else

                cright(1:nof_variables) = iexsolhir(ilocal_recon3(i)%ihexn(1, ielem(n, i)%indexi(l)))%sol &
                                          (ilocal_recon3(i)%ihexl(1, ielem(n, i)%indexi(l)), 1:nof_variables)
              end if

              if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
                if (fastest .eq. 1) then
                  cturbr(1:turbulenceequations + passivescalar) = solchanger(ielem(n, i)%ineighn(l))%sol &
                                                        (ielem(n, i)%q_face(l)%q_mapl(1), 6:5 + turbulenceequations + passivescalar)
                else

                  cturbr(1:turbulenceequations + passivescalar) = iexsolhir(ilocal_recon3(i)%ihexn(1, ielem(n, i)%indexi(l)))%sol &
                                       (ilocal_recon3(i)%ihexl(1, ielem(n, i)%indexi(l)), 6:5 + turbulenceequations + passivescalar)
                end if
              end if

            end if
          else

            if (fastest .eq. 1) then
              cright(1:nof_variables) = solchanger(ielem(n, i)%ineighn(l))%sol(ielem(n, i)%q_face(l)%q_mapl(1), 1:nof_variables)
            else

              cright(1:nof_variables) = iexsolhir(ilocal_recon3(i)%ihexn(1, ielem(n, i)%indexi(l)))%sol &
                                        (ilocal_recon3(i)%ihexl(1, ielem(n, i)%indexi(l)), 1:nof_variables)
            end if

            if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
              if (fastest .eq. 1) then
                cturbr(1:turbulenceequations + passivescalar) = solchanger(ielem(n, i)%ineighn(l))%sol &
                                                        (ielem(n, i)%q_face(l)%q_mapl(1), 6:5 + turbulenceequations + passivescalar)
              else

                cturbr(1:turbulenceequations + passivescalar) = iexsolhir(ilocal_recon3(i)%ihexn(1, ielem(n, i)%indexi(l)))%sol &
                                       (ilocal_recon3(i)%ihexl(1, ielem(n, i)%indexi(l)), 6:5 + turbulenceequations + passivescalar)
              end if
            end if

!
          end if
        end if

        call rotatef(n, cright_rot, cright, angle1, angle2)        !rotate wrt to normalvector of face and solve 1d riemann problem
        call rotatef(n, cleft_rot, cleft, angle1, angle2)        !rotate wrt to normalvector of face and solve 1d riemann problem

        if ((lmach .eq. 1)) then    !application of the low mach number correction
          leftv(1:nof_variables) = cleft_rot(1:nof_variables); rightv(1:nof_variables) = cright_rot(1:nof_variables)
          call lmacht(n, leftv, rightv)
          cleft_rot(1:nof_variables) = leftv(1:nof_variables); cright_rot(1:nof_variables) = rightv(1:nof_variables);
          call rotatef(n, cright, cright_rot, angle1, angle2)        !rotate wrt to normalvector of face and solve 1d riemann problem
          call rotatef(n, cleft, cleft_rot, angle1, angle2)

        end if

        leftv(1:nof_variables) = cleft(1:nof_variables); rightv(1:nof_variables) = cright(1:nof_variables)
        call cons2prim2(n, leftv, rightv, mp_pinfl, mp_pinfr, gammal, gammar)

        asound1 = sqrt(leftv(5)*gamma/leftv(1)) + abs(cleft_rot(2)/cleft_rot(1))
        asound2 = sqrt(rightv(5)*gamma/rightv(1)) + abs(cright_rot(2)/cright_rot(1))

        vpp = max(asound1, asound2)

        if (itestcase .eq. 4) then
          call sutherland(n, leftv, rightv, viscl, laml)

          viscots = (viscl(1) + viscl(2))*oo2
          mul1 = ielem(n, i)%surf(l)

          viscots = max((4.0/(3.0*0.5*(cleft(1) + cright(1))))*viscots*mul1 &
                        /ielem(n, i)%dih(l), ((gamma/(0.5*(cleft(1) + cright(1)))) &
                                              *viscots*mul1/(ielem(n, i)%dih(l)*prandtl)))
          vpp = max(asound1, asound2) + viscots

          if (turbulence .eq. 1) then
            if (turbulencemodel .eq. 1) then
              turbmv(1) = cturbl(1); turbmv(2) = cturbr(1); eddyfl(2) = turbmv(1); eddyfr(2) = turbmv(2)
              call eddyvisco(n, viscl, laml, turbmv, etvm, eddyfl, eddyfr, leftv, rightv)
            end if
            if (turbulencemodel .eq. 2) then
              eddyfl(1) = ielem(n, i)%walldist; eddyfl(2) = cturbl(1); eddyfl(3) = cturbl(2)

              eddyfl(4:6) = ilocal_recon3(i)%grads(1, 1:3); eddyfl(7:9) = ilocal_recon3(i)%grads(2, 1:3)
              eddyfl(10:12) = ilocal_recon3(i)%grads(3, 1:3); eddyfl(13:15) = ilocal_recon3(i)%grads(5, 1:3)
              eddyfl(16:18) = ilocal_recon3(i)%grads(6, 1:3)

              eddyfr = eddyfl
              call eddyvisco(n, viscl, laml, turbmv, etvm, eddyfl, eddyfr, leftv, rightv)
            end if

            viscots = oo2*((viscl(1) + viscl(3)) + (viscl(2) + viscl(4)))

            viscots = max((4.0/(3.0*0.5*(cleft(1) + cright(1))))*viscots* &
                          mul1/ielem(n, i)%dih(l), ((gamma/(0.5*(cleft(1) + cright(1))))* &
                                                    viscots*mul1/(ielem(n, i)%dih(l)*(prandtl + prtu))))

            vpp = max(asound1, asound2) + viscots
          end if

          impdiag(1, 1:nof_variables, 1:nof_variables) = impdiag(1, 1:nof_variables, 1:nof_variables) + (oo2*((vpp*identity1))*mul1)
          call compute_jacobianse(n, iconsidered, eigvl, cright, gamma, angle1, angle2, srf_speedrot, nx, ny, nz)
          convj = eigvl
                                                  impoff(1,l,1:nof_variables,1:nof_variables)=impoff(1,l,1:nof_variables,1:nof_variables)+(((oo2*convj(1:nof_variables,1:nof_variables))&
                                                                                                    - ((oo2*vpp)*identity1))*mul1)

          if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
            if (turbulence .eq. 1) then
            do nvar = 1, turbulenceequations

              vpp = max(asound1, asound2) + viscots
              impdiagt(1, nvar) = impdiagt(1, nvar) + (oo2*((vpp))*mul1)
              impofft(1, l, nvar) = impofft(1, l, nvar) - (((oo2*vpp))*mul1)
            end do
            end if
            if (passivescalar .gt. 0) then
            do nvar = turbulenceequations + 1, turbulenceequations + passivescalar
              viscl(1) = viscl(1)/schmidt_lam
              viscl(2) = viscl(2)/schmidt_lam
              viscl(3) = viscl(3)/schmidt_turb
              viscl(4) = viscl(4)/schmidt_turb
              viscots = 0.5*((viscl(1) + viscl(3)) + (viscl(2) + viscl(4)))

              viscots = (2.0*viscots)/((cleft(1) + cright(1))*ielem(n, i)%dih(l))
              vpp = max(asound1, asound2) + viscots
              impdiagt(1, nvar) = impdiagt(1, nvar) + (oo2*((vpp))*mul1)
              impofft(1, l, nvar) = impofft(1, l, nvar) - (oo2*((vpp))*mul1)
            end do
            end if
          end if
        else

          impdiag(1, 1:nof_variables, 1:nof_variables) = impdiag(1, 1:nof_variables, 1:nof_variables) + (oo2*((vpp*identity1))*mul1)
          call compute_jacobianse(n, iconsidered, eigvl, cright, gamma, angle1, angle2, srf_speedrot, nx, ny, nz)
          convj = eigvl
                                                  impoff(1,l,1:nof_variables,1:nof_variables)=impoff(1,l,1:nof_variables,1:nof_variables)+(((oo2*convj(1:nof_variables,1:nof_variables))&
                                                                                                    - ((oo2*vpp)*identity1))*mul1)

        end if

      end do
    end if

    if (rungekutta .eq. 10) then

      impdiag(1, 1, 1) = impdiag(1, 1, 1) + (ielem(n, i)%totvolume/ielem(n, i)%dtl)
      impdiag(1, 2, 2) = impdiag(1, 2, 2) + (ielem(n, i)%totvolume/ielem(n, i)%dtl)
      impdiag(1, 3, 3) = impdiag(1, 3, 3) + (ielem(n, i)%totvolume/ielem(n, i)%dtl)
      impdiag(1, 4, 4) = impdiag(1, 4, 4) + (ielem(n, i)%totvolume/ielem(n, i)%dtl)
      impdiag(1, 5, 5) = impdiag(1, 5, 5) + (ielem(n, i)%totvolume/ielem(n, i)%dtl)

    else

!             impdiag(1,1,1)=ielem(n,i)%totvolume*( ((dt+1.5d0*ielem(n,i)%dtl)/dt) +(impdiag(1,1,1)*ielem(n,i)%dtl/ielem(n,i)%totvolume))/ielem(n,i)%dtl
!             impdiag(1,2,2)=ielem(n,i)%totvolume*(((dt+1.5d0*ielem(n,i)%dtl)/dt)+(impdiag(1,2,2)*ielem(n,i)%dtl/ielem(n,i)%totvolume))/ielem(n,i)%dtl
!             impdiag(1,3,3)=ielem(n,i)%totvolume*(((dt+1.5d0*ielem(n,i)%dtl)/dt)+(impdiag(1,3,3)*ielem(n,i)%dtl/ielem(n,i)%totvolume))/ielem(n,i)%dtl
!             impdiag(1,4,4)=ielem(n,i)%totvolume*(((dt+1.5d0*ielem(n,i)%dtl)/dt)+(impdiag(1,4,4)*ielem(n,i)%dtl/ielem(n,i)%totvolume))/ielem(n,i)%dtl
!             impdiag(1,5,5)=ielem(n,i)%totvolume*(((dt+1.5d0*ielem(n,i)%dtl)/dt)+(impdiag(1,5,5)*ielem(n,i)%dtl/ielem(n,i)%totvolume))/ielem(n,i)%dtl

      impdiag(1, 1, 1) = ielem(n, i)%totvolume*((1.0d0/ielem(n, i)%dtl) + (1.5d0/dt)) + (impdiag(1, 1, 1))
      impdiag(1, 2, 2) = ielem(n, i)%totvolume*((1.0d0/ielem(n, i)%dtl) + (1.5d0/dt)) + (impdiag(1, 2, 2))
      impdiag(1, 3, 3) = ielem(n, i)%totvolume*((1.0d0/ielem(n, i)%dtl) + (1.5d0/dt)) + (impdiag(1, 3, 3))
      impdiag(1, 4, 4) = ielem(n, i)%totvolume*((1.0d0/ielem(n, i)%dtl) + (1.5d0/dt)) + (impdiag(1, 4, 4))
      impdiag(1, 5, 5) = ielem(n, i)%totvolume*((1.0d0/ielem(n, i)%dtl) + (1.5d0/dt)) + (impdiag(1, 5, 5))

    end if

    if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
      call sources_derivatives(n, iconsidered, source_t)
      sht(i, 1:turbulenceequations) = (source_t(1:turbulenceequations)*ielem(n, i)%totvolume)
      if (rungekutta .eq. 10) then

        if (turbulence .eq. 1) then
        do nvar = 1, turbulenceequations
          impdiagt(1, nvar) = impdiagt(1, nvar) + (ielem(n, i)%totvolume/ielem(n, i)%dtl) - sht(i, nvar)
        end do
        end if
        if (passivescalar .gt. 0) then
        do nvar = turbulenceequations + 1, turbulenceequations + passivescalar
          impdiagt(1, nvar) = impdiagt(1, nvar) + (ielem(n, i)%totvolume/ielem(n, i)%dtl)
        end do
        end if

      else

        if (turbulence .eq. 1) then
        do nvar = 1, turbulenceequations
          impdiagt(1, nvar) = (ielem(n, i)%totvolume*((1.0d0/ielem(n, i)%dtl) + ((1.5d0/dt)*impdiagt(1, nvar)))) - sht(i, nvar)
        end do
        end if
        if (passivescalar .gt. 0) then
        do nvar = turbulenceequations + 1, turbulenceequations + passivescalar
          impdiagt(1, nvar) = (ielem(n, i)%totvolume*((1.0d0/ielem(n, i)%dtl) + ((1.5d0/dt)*impdiagt(1, nvar))))
        end do
        end if

      end if
    end if

  end subroutine calculate_jacobianlm

  subroutine calculate_jacobian_2dlm(n, iconsidered, impdiag, impdiagt, impoff, impofft)
    implicit none
    integer, intent(in)::n, iconsidered
    real, dimension(1:nof_variables + turbulenceequations + passivescalar)::godflux2
    integer::i, l, ngp, kmaxe, iqp, ii, nvar, n_node, ibfc
    real::sum_detect, norms, vpp, asound1, asound2, mul1, dxb, tempxx, viscots
    real, dimension(nof_variables, nof_variables)::identity1
    real, dimension(nof_variables, nof_variables)::convj, diffj
    integer::facex, pointx, igoflux
    integer::b_code
    real::angle1, angle2, nx, ny, nz
    real, dimension(1:nof_variables)::cleft, cright, cright_rot, cleft_rot
    real, dimension(1:turbulenceequations + passivescalar)::cturbl, cturbr
    real, dimension(1:nof_variables)::leftv, srf_speedrot, srf_speed
    real, dimension(1:nof_variables)::rightv
    real, dimension(1:dimensiona)::pox, poy, poz
    real, dimension(1:8, 1:dimensiona)::vext, nodes_list
    real, dimension(1:4)::viscl, laml
    real, dimension(1:2)::turbmv
    real, dimension(1)::etvm
    real, dimension(1:20)::eddyfl, eddyfr
    real, dimension(1:dimensiona)::cords
    real::mp_pinfl, gammal
    real::mp_pinfr, gammar
    real, allocatable, dimension(:, :), intent(inout)::impdiagt
    real, allocatable, dimension(:, :, :), intent(inout)::impdiag, impofft
    real, allocatable, dimension(:, :, :, :), intent(inout)::impoff
    real, dimension(1:nof_variables, 1:nof_variables)::eigvl
    real, dimension(turbulenceequations)::source_t

    identity1(:, :) = zero
    identity1(1, 1) = 1.0d0
    identity1(2, 2) = 1.0d0
    identity1(3, 3) = 1.0d0
    identity1(4, 4) = 1.0d0

    if (ielem(n, iconsidered)%interior .eq. 0) then

      i = iconsidered
      impdiag(1, :, :) = zero
      impoff(1, :, :, :) = zero
      if (turbulence .eq. 1) then
        impdiagt(1, :) = zero
        impofft(1, :, :) = zero
      end if

      do l = 1, ielem(n, i)%ifca !for all their faces
        godflux2 = zero
        angle1 = ielem(n, i)%faceanglex(l)
        angle2 = ielem(n, i)%faceangley(l)
        nx = angle1
        ny = angle2
        mul1 = ielem(n, i)%surf(l)

        cleft(1:nof_variables) = u_c(i)%val(1, 1:nof_variables)
        cright(1:nof_variables) = u_c(ielem(n, i)%ineigh(l))%val(1, 1:nof_variables)

        if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

          cturbl(1:turbulenceequations + passivescalar) = u_ct(i)%val(1, 1:turbulenceequations + passivescalar)!left additional equations flow state
          cturbr(1:turbulenceequations + passivescalar) = u_ct(ielem(n, i)%ineigh(l))%val(1, 1:turbulenceequations + passivescalar)

        end if

        call rotatef2d(n, cright_rot, cright, angle1, angle2)        !rotate wrt to normalvector of face and solve 1d riemann problem
        call rotatef2d(n, cleft_rot, cleft, angle1, angle2)        !rotate wrt to normalvector of face and solve 1d riemann problem

        if ((lmach .eq. 1)) then    !application of the low mach number correction
          leftv(1:nof_variables) = cleft_rot(1:nof_variables); rightv(1:nof_variables) = cright_rot(1:nof_variables)
          call lmacht2d(n, leftv, rightv)
          cleft_rot(1:nof_variables) = leftv(1:nof_variables); cright_rot(1:nof_variables) = rightv(1:nof_variables);
          call rotateb2d(n, cright, cright_rot, angle1, angle2)        !rotate wrt to normalvector of face and solve 1d riemann problem
          call rotateb2d(n, cleft, cleft_rot, angle1, angle2)        !rotate wrt to normalvector of face and so
        end if

        leftv(1:nof_variables) = cleft(1:nof_variables); rightv(1:nof_variables) = cright(1:nof_variables)
        call cons2prim2(n, leftv, rightv, mp_pinfl, mp_pinfr, gammal, gammar)

        asound1 = sqrt(leftv(4)*gamma/leftv(1)) + abs(cleft_rot(2)/cleft_rot(1))
        asound2 = sqrt(rightv(4)*gamma/rightv(1)) + abs(cright_rot(2)/cright_rot(1))

        vpp = max(asound1, asound2)

        if (itestcase .eq. 4) then
          call sutherland2d(n, leftv, rightv, viscl, laml)

          viscots = (viscl(1) + viscl(2))*oo2
          mul1 = ielem(n, i)%surf(l)

          viscots = max((4.0/(3.0*0.5*(cleft(1) + cright(1))))*viscots*mul1 &
                        /ielem(n, i)%dih(l), ((gamma/(0.5*(cleft(1) + cright(1)))) &
                                              *viscots*mul1/(ielem(n, i)%dih(l)*prandtl)))
          vpp = max(asound1, asound2) + viscots

          if (turbulence .eq. 1) then
            if (turbulencemodel .eq. 1) then
              turbmv(1) = cturbl(1); turbmv(2) = cturbr(1); eddyfl(2) = turbmv(1); eddyfr(2) = turbmv(2)
              call eddyvisco2d(n, viscl, laml, turbmv, etvm, eddyfl, eddyfr, leftv, rightv)
            end if
            if (turbulencemodel .eq. 2) then
              eddyfl(1) = ielem(n, i)%walldist; eddyfl(2) = cturbl(1); eddyfl(3) = cturbl(2)
              eddyfl(4:5) = ilocal_recon3(i)%grads(1, 1:2); eddyfl(6:7) = ilocal_recon3(i)%grads(2, 1:2)
              eddyfl(8:9) = ilocal_recon3(i)%grads(4, 1:2)
              eddyfl(10:11) = ilocal_recon3(i)%grads(5, 1:2)

              eddyfr = eddyfl
              call eddyvisco2d(n, viscl, laml, turbmv, etvm, eddyfl, eddyfr, leftv, rightv)
            end if

            viscots = oo2*((viscl(1) + viscl(3)) + (viscl(2) + viscl(4)))

            viscots = max((4.0/(3.0*0.5*(cleft(1) + cright(1))))*viscots* &
                          mul1/ielem(n, i)%dih(l), ((gamma/(0.5*(cleft(1) + cright(1))))* &
                                                    viscots*mul1/(ielem(n, i)%dih(l)*(prandtl + prtu))))

            vpp = max(asound1, asound2) + viscots
          end if

          impdiag(1, 1:nof_variables, 1:nof_variables) = impdiag(1, 1:nof_variables, 1:nof_variables) + (oo2*((vpp*identity1))*mul1)
          call compute_jacobianse2d(n, eigvl, cright, gamma, angle1, angle2, nx, ny, nz)
          convj = eigvl
                                                  impoff(1,l,1:nof_variables,1:nof_variables)=impoff(1,l,1:nof_variables,1:nof_variables)+(((oo2*convj(1:nof_variables,1:nof_variables))&
                                                                                                    - ((oo2*vpp)*identity1))*mul1)

          if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
            if (turbulence .eq. 1) then
            do nvar = 1, turbulenceequations

              vpp = max(asound1, asound2) + viscots
              impdiagt(1, nvar) = impdiagt(1, nvar) + (oo2*((vpp))*mul1)
              impofft(1, l, nvar) = impofft(1, l, nvar) - (((oo2*vpp))*mul1)
            end do
            end if
            if (passivescalar .gt. 0) then
            do nvar = turbulenceequations + 1, turbulenceequations + passivescalar
              viscl(1) = viscl(1)/schmidt_lam
              viscl(2) = viscl(2)/schmidt_lam
              viscl(3) = viscl(3)/schmidt_turb
              viscl(4) = viscl(4)/schmidt_turb
              viscots = 0.5*((viscl(1) + viscl(3)) + (viscl(2) + viscl(4)))

              viscots = (2.0*viscots)/((cleft(1) + cright(1))*ielem(n, i)%dih(l))
              vpp = max(asound1, asound2) + viscots
              impdiagt(1, nvar) = impdiagt(1, nvar) + (oo2*((vpp))*mul1)
              impofft(1, l, nvar) = impofft(1, l, nvar) - (oo2*((vpp))*mul1)
            end do
            end if
          end if
        else

          impdiag(1, 1:nof_variables, 1:nof_variables) = impdiag(1, 1:nof_variables, 1:nof_variables) + (oo2*((vpp*identity1))*mul1)
          call compute_jacobianse2d(n, eigvl, cright, gamma, angle1, angle2, nx, ny, nz)
          convj = eigvl
                                                  impoff(1,l,1:nof_variables,1:nof_variables)=impoff(1,l,1:nof_variables,1:nof_variables)+(((oo2*convj(1:nof_variables,1:nof_variables))&
                                                                                                    - ((oo2*vpp)*identity1))*mul1)

        end if
      end do
    else

      impdiag(1, :, :) = 0.0
      impoff(1, :, :, :) = 0.0
      if (turbulence .eq. 1) then
        impdiagt(1, :) = 0.0
        impofft(1, :, :) = 0.0
      end if

      do l = 1, ielem(n, i)%ifca
        mul1 = ielem(n, i)%surf(l)
        angle1 = ielem(n, i)%faceanglex(l)
        angle2 = ielem(n, i)%faceangley(l)
        nx = angle1
        ny = angle2

        b_code = 0
        cleft(1:nof_variables) = u_c(i)%val(1, 1:nof_variables)
        if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

          cturbl(1:turbulenceequations + passivescalar) = u_ct(i)%val(1, 1:turbulenceequations + passivescalar)

        end if

        if (ielem(n, i)%ineighb(l) .eq. n) then        !my cpu only
          if (ielem(n, i)%ibounds(l) .gt. 0) then        !check for boundaries
            if (ibound(n, ielem(n, i)%ibounds(l))%icode .eq. 5) then        !periodic in my cpu
              cright(1:nof_variables) = u_c(ielem(n, i)%ineigh(l))%val(1, 1:nof_variables)

              if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

           cturbr(1:turbulenceequations + passivescalar) = u_ct(ielem(n, i)%ineigh(l))%val(1, 1:turbulenceequations + passivescalar)

              end if

            else


              facex = l;
              call coordinates_face_inner2dx(n, iconsidered, facex, vext, nodes_list)
              n_node = 2
              cords(1:2) = zero
              cords(1:2) = cordinates2(n, nodes_list, n_node)

              poy(1) = cords(2)
              pox(1) = cords(1)

              leftv(1:nof_variables) = cleft(1:nof_variables)
              b_code = ibound(n, ielem(n, i)%ibounds(l))%icode

                                                                    call boundarys2d(n,b_code,iconsidered,facex,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speed,srf_speedrot,ibfc)
              cright(1:nof_variables) = rightv(1:nof_variables)

            end if
          else
            cright(1:nof_variables) = u_c(ielem(n, i)%ineigh(l))%val(1, 1:nof_variables)

            if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

           cturbr(1:turbulenceequations + passivescalar) = u_ct(ielem(n, i)%ineigh(l))%val(1, 1:turbulenceequations + passivescalar)

            end if

          end if
        else        !in other cpus they can only be periodic or mpi neighbours

          if (ielem(n, i)%ibounds(l) .gt. 0) then        !check for boundaries
            if (ibound(n, ielem(n, i)%ibounds(l))%icode .eq. 5) then        !periodic in other cpu

              if (fastest .eq. 1) then
                cright(1:nof_variables) = solchanger(ielem(n, i)%ineighn(l))%sol(ielem(n, i)%q_face(l)%q_mapl(1), 1:nof_variables)
              else

                cright(1:nof_variables) = iexsolhir(ilocal_recon3(i)%ihexn(1, ielem(n, i)%indexi(l)))%sol &
                                          (ilocal_recon3(i)%ihexl(1, ielem(n, i)%indexi(l)), 1:nof_variables)
              end if

              if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
                if (fastest .eq. 1) then
                  cturbr(1:turbulenceequations + passivescalar) = solchanger(ielem(n, i)%ineighn(l))%sol &
                                                        (ielem(n, i)%q_face(l)%q_mapl(1), 5:4 + turbulenceequations + passivescalar)
                else

                  cturbr(1:turbulenceequations + passivescalar) = iexsolhir(ilocal_recon3(i)%ihexn(1, ielem(n, i)%indexi(l)))%sol &
                                       (ilocal_recon3(i)%ihexl(1, ielem(n, i)%indexi(l)), 5:4 + turbulenceequations + passivescalar)
                end if
              end if

            end if
          else

            if (fastest .eq. 1) then
              cright(1:nof_variables) = solchanger(ielem(n, i)%ineighn(l))%sol(ielem(n, i)%q_face(l)%q_mapl(1), 1:nof_variables)
            else

              cright(1:nof_variables) = iexsolhir(ilocal_recon3(i)%ihexn(1, ielem(n, i)%indexi(l)))%sol &
                                        (ilocal_recon3(i)%ihexl(1, ielem(n, i)%indexi(l)), 1:nof_variables)
            end if

            if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
              if (fastest .eq. 1) then
                cturbr(1:turbulenceequations + passivescalar) = solchanger(ielem(n, i)%ineighn(l))%sol &
                                                        (ielem(n, i)%q_face(l)%q_mapl(1), 5:4 + turbulenceequations + passivescalar)
              else

                cturbr(1:turbulenceequations + passivescalar) = iexsolhir(ilocal_recon3(i)%ihexn(1, ielem(n, i)%indexi(l)))%sol &
                                       (ilocal_recon3(i)%ihexl(1, ielem(n, i)%indexi(l)), 5:4 + turbulenceequations + passivescalar)
              end if
            end if

!
          end if
        end if

        call rotatef2d(n, cright_rot, cright, angle1, angle2)        !rotate wrt to normalvector of face and solve 1d riemann problem
        call rotatef2d(n, cleft_rot, cleft, angle1, angle2)        !rotate wrt to normalvector of face and solve 1d riemann problem

        if ((lmach .eq. 1)) then    !application of the low mach number correction
          leftv(1:nof_variables) = cleft_rot(1:nof_variables); rightv(1:nof_variables) = cright_rot(1:nof_variables)
          call lmacht2d(n, leftv, rightv)
          cleft_rot(1:nof_variables) = leftv(1:nof_variables); cright_rot(1:nof_variables) = rightv(1:nof_variables);
          call rotateb2d(n, cright, cright_rot, angle1, angle2)        !rotate wrt to normalvector of face and solve 1d riemann problem
          call rotateb2d(n, cleft, cleft_rot, angle1, angle2)        !rotate wrt to normalvector of face and so
        end if

        leftv(1:nof_variables) = cleft(1:nof_variables); rightv(1:nof_variables) = cright(1:nof_variables)
        call cons2prim2(n, leftv, rightv, mp_pinfl, mp_pinfr, gammal, gammar)

        asound1 = sqrt(leftv(4)*gamma/leftv(1)) + abs(cleft_rot(2)/cleft_rot(1))
        asound2 = sqrt(rightv(4)*gamma/rightv(1)) + abs(cright_rot(2)/cright_rot(1))

        vpp = max(asound1, asound2)

        if (itestcase .eq. 4) then
          call sutherland2d(n, leftv, rightv, viscl, laml)

          viscots = (viscl(1) + viscl(2))*oo2
          mul1 = ielem(n, i)%surf(l)

          viscots = max((4.0/(3.0*0.5*(cleft(1) + cright(1))))*viscots*mul1 &
                        /ielem(n, i)%dih(l), ((gamma/(0.5*(cleft(1) + cright(1)))) &
                                              *viscots*mul1/(ielem(n, i)%dih(l)*prandtl)))
          vpp = max(asound1, asound2) + viscots

          if (turbulence .eq. 1) then
            if (turbulencemodel .eq. 1) then
              turbmv(1) = cturbl(1); turbmv(2) = cturbr(1); eddyfl(2) = turbmv(1); eddyfr(2) = turbmv(2)
              call eddyvisco2d(n, viscl, laml, turbmv, etvm, eddyfl, eddyfr, leftv, rightv)
            end if
            if (turbulencemodel .eq. 2) then
              eddyfl(1) = ielem(n, i)%walldist; eddyfl(2) = cturbl(1); eddyfl(3) = cturbl(2)
              eddyfl(4:5) = ilocal_recon3(i)%grads(1, 1:2); eddyfl(6:7) = ilocal_recon3(i)%grads(2, 1:2)
              eddyfl(8:9) = ilocal_recon3(i)%grads(4, 1:2)
              eddyfl(10:11) = ilocal_recon3(i)%grads(5, 1:2)

              eddyfr = eddyfl
              call eddyvisco2d(n, viscl, laml, turbmv, etvm, eddyfl, eddyfr, leftv, rightv)
            end if

            viscots = oo2*((viscl(1) + viscl(3)) + (viscl(2) + viscl(4)))

            viscots = max((4.0/(3.0*0.5*(cleft(1) + cright(1))))*viscots* &
                          mul1/ielem(n, i)%dih(l), ((gamma/(0.5*(cleft(1) + cright(1))))* &
                                                    viscots*mul1/(ielem(n, i)%dih(l)*(prandtl + prtu))))

            vpp = max(asound1, asound2) + viscots
          end if

          impdiag(1, 1:nof_variables, 1:nof_variables) = impdiag(1, 1:nof_variables, 1:nof_variables) + (oo2*((vpp*identity1))*mul1)
          call compute_jacobianse2d(n, eigvl, cright, gamma, angle1, angle2, nx, ny, nz)
          convj = eigvl
                                                  impoff(1,l,1:nof_variables,1:nof_variables)=impoff(1,l,1:nof_variables,1:nof_variables)+(((oo2*convj(1:nof_variables,1:nof_variables))&
                                                                                                    - ((oo2*vpp)*identity1))*mul1)

          if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
            if (turbulence .eq. 1) then
            do nvar = 1, turbulenceequations

              vpp = max(asound1, asound2) + viscots
              impdiagt(i, nvar) = impdiagt(i, nvar) + (oo2*((vpp))*mul1)
              impofft(i, l, nvar) = impofft(i, l, nvar) - (((oo2*vpp))*mul1)
            end do
            end if
            if (passivescalar .gt. 0) then
            do nvar = turbulenceequations + 1, turbulenceequations + passivescalar
              viscl(1) = viscl(1)/schmidt_lam
              viscl(2) = viscl(2)/schmidt_lam
              viscl(3) = viscl(3)/schmidt_turb
              viscl(4) = viscl(4)/schmidt_turb
              viscots = 0.5*((viscl(1) + viscl(3)) + (viscl(2) + viscl(4)))

              viscots = (2.0*viscots)/((cleft(1) + cright(1))*ielem(n, i)%dih(l))
              vpp = max(asound1, asound2) + viscots
              impdiagt(1, nvar) = impdiagt(1, nvar) + (oo2*((vpp))*mul1)
              impofft(1, l, nvar) = impofft(1, l, nvar) - (oo2*((vpp))*mul1)
            end do
            end if
          end if
        else

          impdiag(1, 1:nof_variables, 1:nof_variables) = impdiag(1, 1:nof_variables, 1:nof_variables) + (oo2*((vpp*identity1))*mul1)
          call compute_jacobianse2d(n, eigvl, cright, gamma, angle1, angle2, nx, ny, nz)
          convj = eigvl
                                                  impoff(1,l,1:nof_variables,1:nof_variables)=impoff(1,l,1:nof_variables,1:nof_variables)+(((oo2*convj(1:nof_variables,1:nof_variables))&
                                                                                                    - ((oo2*vpp)*identity1))*mul1)

        end if

      end do
    end if

    if (rungekutta .eq. 10) then

      impdiag(1, 1, 1) = impdiag(1, 1, 1) + (ielem(n, i)%totvolume/ielem(n, i)%dtl)
      impdiag(1, 2, 2) = impdiag(1, 2, 2) + (ielem(n, i)%totvolume/ielem(n, i)%dtl)
      impdiag(1, 3, 3) = impdiag(1, 3, 3) + (ielem(n, i)%totvolume/ielem(n, i)%dtl)
      impdiag(1, 4, 4) = impdiag(1, 4, 4) + (ielem(n, i)%totvolume/ielem(n, i)%dtl)

    else
      impdiag(1, 1, 1) = ielem(n, i)%totvolume*((1.0d0/ielem(n, i)%dtl) + (1.5d0/dt)) + (impdiag(1, 1, 1))
      impdiag(1, 2, 2) = ielem(n, i)%totvolume*((1.0d0/ielem(n, i)%dtl) + (1.5d0/dt)) + (impdiag(1, 2, 2))
      impdiag(1, 3, 3) = ielem(n, i)%totvolume*((1.0d0/ielem(n, i)%dtl) + (1.5d0/dt)) + (impdiag(1, 3, 3))
      impdiag(1, 4, 4) = ielem(n, i)%totvolume*((1.0d0/ielem(n, i)%dtl) + (1.5d0/dt)) + (impdiag(1, 4, 4))

    end if

    if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then
      call sources_derivatives2d(n, iconsidered, source_t)
      sht(i, 1:turbulenceequations) = (source_t(1:turbulenceequations)*ielem(n, i)%totvolume)
      if (rungekutta .eq. 10) then

        if (turbulence .eq. 1) then
        do nvar = 1, turbulenceequations
          impdiagt(1, nvar) = impdiagt(1, nvar) + (ielem(n, i)%totvolume/ielem(n, i)%dtl) - sht(i, nvar)
        end do
        end if
        if (passivescalar .gt. 0) then
        do nvar = turbulenceequations + 1, turbulenceequations + passivescalar

          impdiagt(1, nvar) = impdiagt(1, nvar) + (ielem(n, i)%totvolume/ielem(n, i)%dtl)
        end do
        end if

      else

        if (turbulence .eq. 1) then
        do nvar = 1, turbulenceequations
          impdiagt(1, nvar) = ielem(n, i)%totvolume*((1.0d0/ielem(n, i)%dtl) + (1.5d0/dt)) + (impdiagt(1, 1)) - sht(i, nvar)

        end do
        end if
        if (passivescalar .gt. 0) then
        do nvar = turbulenceequations + 1, turbulenceequations + passivescalar
          impdiagt(1, nvar) = ielem(n, i)%totvolume*((1.0d0/ielem(n, i)%dtl) + (1.5d0/dt)) + (impdiagt(1, 1))
        end do
        end if

      end if
    end if

  end subroutine calculate_jacobian_2dlm

  subroutine calculate_jacobian_2d_mf(n)
    implicit none
    integer, intent(in)::n
    real, dimension(1:nof_variables + turbulenceequations + passivescalar)::godflux2
    integer::i, l, ngp, kmaxe, iqp, ii, nvar, n_node, ibfc, kas
    real::sum_detect, norms, vpp, asound1, asound2, mul1, dxb, tempxx, viscots
    real, dimension(nof_variables, nof_variables)::identity1
    real, dimension(nof_variables, nof_variables)::convj, diffj
    integer::iconsidered, facex, pointx, igoflux
    integer::b_code
    real::angle1, angle2, nx, ny, nz
    real, dimension(1:nof_variables)::cleft, cright, cright_rot, cleft_rot
    real, dimension(1:turbulenceequations + passivescalar)::cturbl, cturbr
    real, dimension(1:nof_variables)::leftv, srf_speedrot, srf_speed
    real, dimension(1:nof_variables)::rightv
    real, dimension(1:dimensiona)::pox, poy, poz

    real, dimension(1:8, 1:dimensiona)::vext, nodes_list
    real, dimension(1:dimensiona)::cords
    real, dimension(1:4)::viscl, laml
    real, dimension(1:2)::turbmv
    real, dimension(1)::etvm
    real, dimension(1:20)::eddyfl, eddyfr

    real::mp_pinfl, gammal
    real::mp_pinfr, gammar
    real, dimension(1:nof_variables, 1:nof_variables)::eigvl
    kmaxe = xmpielrank(n)

    !$omp do
    do ii = 1, nof_interior        !for all the interior elements
      i = el_int(ii)
      iconsidered = i
      impdiag_mf(i) = zero
      impoff_mf(i, :) = zero
      if (turbulence .eq. 1) then
        impdiagt(i, 1:turbulenceequations + passivescalar) = 0.0
        impofft(i, :, :) = zero
      end if
      do l = 1, ielem(n, i)%ifca !for all their faces
        b_code = 0
        angle1 = ielem(n, i)%faceanglex(l)
        angle2 = ielem(n, i)%faceangley(l)
        nx = angle1
        ny = angle2
        mul1 = ielem(n, i)%surf(l)

        cleft(1:nof_variables) = u_c(i)%val(1, 1:nof_variables)
        cright(1:nof_variables) = u_c(ielem(n, i)%ineigh(l))%val(1, 1:nof_variables)

        if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

          cturbl(1:turbulenceequations + passivescalar) = u_ct(i)%val(1, 1:turbulenceequations + passivescalar)!left additional equations flow state
          cturbr(1:turbulenceequations + passivescalar) = u_ct(ielem(n, i)%ineigh(l))%val(1, 1:turbulenceequations + passivescalar)

        end if

        call rotatef2d(n, cright_rot, cright, angle1, angle2)        !rotate wrt to normalvector of face and solve 1d riemann problem
        call rotatef2d(n, cleft_rot, cleft, angle1, angle2)        !rotate wrt to normalvector of face and solve 1d riemann problem

        if ((lmach .eq. 1)) then    !application of the low mach number correction
          leftv(1:nof_variables) = cleft_rot(1:nof_variables); rightv(1:nof_variables) = cright_rot(1:nof_variables)
          call lmacht2d(n, leftv, rightv)
          cleft_rot(1:nof_variables) = leftv(1:nof_variables); cright_rot(1:nof_variables) = rightv(1:nof_variables);
          call rotateb2d(n, cright, cright_rot, angle1, angle2)        !rotate wrt to normalvector of face and solve 1d riemann problem
          call rotateb2d(n, cleft, cleft_rot, angle1, angle2)

        end if

        leftv(1:nof_variables) = cleft(1:nof_variables); rightv(1:nof_variables) = cright(1:nof_variables)
        call cons2prim2(n, leftv, rightv, mp_pinfl, mp_pinfr, gammal, gammar)

        asound1 = sqrt(leftv(4)*gamma/leftv(1)) + abs(cleft_rot(2)/cleft_rot(1))
        asound2 = sqrt(rightv(4)*gamma/rightv(1)) + abs(cright_rot(2)/cright_rot(1))

        vpp = max(asound1, asound2)

        if (itestcase .eq. 4) then
          call sutherland2d(n, leftv, rightv, viscl, laml)

          viscots = (viscl(1) + viscl(2))*oo2

!                                                   viscots=viscots/((0.5*(cleft(1)+cright(1)))*ielem(n,i)%dih(l))
          viscots = max((4.0/(3.0*0.5*(cleft(1) + cright(1))))*viscots*mul1 &
                        /ielem(n, i)%dih(l), ((gamma/(0.5*(cleft(1) + cright(1)))) &
                                              *viscots*mul1/(ielem(n, i)%dih(l)*prandtl)))
          vpp = max(asound1, asound2) + viscots

          if (turbulence .eq. 1) then
            if (turbulencemodel .eq. 1) then
              turbmv(1) = cturbl(1); turbmv(2) = cturbr(1); eddyfl(2) = turbmv(1); eddyfr(2) = turbmv(2)
              call eddyvisco2d(n, viscl, laml, turbmv, etvm, eddyfl, eddyfr, leftv, rightv)
            end if
            if (turbulencemodel .eq. 2) then
              eddyfl(1) = ielem(n, i)%walldist; eddyfl(2) = cturbl(1); eddyfl(3) = cturbl(2)
              eddyfl(4:5) = ilocal_recon3(i)%grads(1, 1:2); eddyfl(6:7) = ilocal_recon3(i)%grads(2, 1:2)
              eddyfl(8:9) = ilocal_recon3(i)%grads(4, 1:2)
              eddyfl(10:11) = ilocal_recon3(i)%grads(5, 1:2)

              eddyfr = eddyfl
              call eddyvisco2d(n, viscl, laml, turbmv, etvm, eddyfl, eddyfr, leftv, rightv)
            end if

            viscots = oo2*((viscl(1) + viscl(3)) + (viscl(2) + viscl(4)))

            !                                                       viscots=viscots/((0.5*(cleft(1)+cright(1)))*ielem(n,i)%dih(l))
            viscots = max((4.0/(3.0*0.5*(cleft(1) + cright(1))))*viscots* &
                          mul1/ielem(n, i)%dih(l), ((gamma/(0.5*(cleft(1) + cright(1))))* &
                                                    viscots*mul1/(ielem(n, i)%dih(l)*(prandtl + prtu))))

            vpp = max(asound1, asound2) + viscots
          end if
        end if

        impdiag_mf(i) = impdiag_mf(i) + (oo2*vpp*mul1)
        impoff_mf(i, l) = vpp

        if (turbulence .eq. 1) then
          impdiagt(i, :) = impdiagt(i, :) + (oo2*vpp*mul1)
          impofft(i, l, :) = impofft(i, l, :) - (oo2*((vpp))*mul1)
        end if

      end do
    end do
    !$omp end do

    !$omp do
    do ii = 1, nof_bounded
      i = el_bnd(ii)
      iconsidered = i

      impdiag_mf(i) = zero
      impoff_mf(i, :) = zero
      if (turbulence .eq. 1) then
        impdiagt(i, 1:turbulenceequations + passivescalar) = zero
        impofft(i, :, :) = zero
      end if
      do l = 1, ielem(n, i)%ifca
        mul1 = ielem(n, i)%surf(l)
        angle1 = ielem(n, i)%faceanglex(l)
        angle2 = ielem(n, i)%faceangley(l)
        nx = angle1
        ny = angle2

        b_code = 0
        cleft(1:nof_variables) = u_c(i)%val(1, 1:nof_variables)
!                                       cleft(1:nof_variables)=ilocal_recon3(i)%uleft(1:nof_variables,l,1)
        if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

          cturbl(1:turbulenceequations + passivescalar) = u_ct(i)%val(1, 1:turbulenceequations + passivescalar)

        end if

        if (ielem(n, i)%ineighb(l) .eq. n) then        !my cpu only
          if (ielem(n, i)%ibounds(l) .gt. 0) then        !check for boundaries
            if (ibound(n, ielem(n, i)%ibounds(l))%icode .eq. 5) then        !periodic in my cpu
              cright(1:nof_variables) = u_c(ielem(n, i)%ineigh(l))%val(1, 1:nof_variables)
!                                                                   cright(1:nof_variables)=ilocal_recon3(ielem(n,i)%ineigh(l))%uleft(1:nof_variables,ielem(n,i)%ineighn(l),1)

              if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

           cturbr(1:turbulenceequations + passivescalar) = u_ct(ielem(n, i)%ineigh(l))%val(1, 1:turbulenceequations + passivescalar)

              end if

              kas = 1

            else
              !not periodic ones in my cpu

              facex = l; iconsidered = i
              call coordinates_face_inner2dx(n, iconsidered, facex, vext, nodes_list)
              n_node = 2
              cords(1:2) = zero
              cords(1:2) = cordinates2(n, nodes_list, n_node)

              poy(1) = cords(2)
              pox(1) = cords(1)

              leftv(1:nof_variables) = cleft(1:nof_variables)
              b_code = ibound(n, ielem(n, i)%ibounds(l))%icode

                                                                    call boundarys2d(n,b_code,iconsidered,facex,leftv,rightv,pox,poy,poz,angle1,angle2,nx,ny,nz,cturbl,cturbr,cright_rot,cleft_rot,srf_speed,srf_speedrot,ibfc)
              cright(1:nof_variables) = rightv(1:nof_variables)

              kas = 2

            end if
          else
            cright(1:nof_variables) = u_c(ielem(n, i)%ineigh(l))%val(1, 1:nof_variables)
!                                                               cright(1:nof_variables)=ilocal_recon3(ielem(n,i)%ineigh(l))%uleft(1:nof_variables,ielem(n,i)%ineighn(l),1)

            if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

           cturbr(1:turbulenceequations + passivescalar) = u_ct(ielem(n, i)%ineigh(l))%val(1, 1:turbulenceequations + passivescalar)

            end if

            kas = 3

          end if
        else        !in other cpus they can only be periodic or mpi neighbours

          if (ielem(n, i)%ibounds(l) .gt. 0) then        !check for boundaries
            if (ibound(n, ielem(n, i)%ibounds(l))%icode .eq. 5) then        !periodic in other cpu

              if (fastest .eq. 1) then
                cright(1:nof_variables) = solchanger(ielem(n, i)%ineighn(l))%sol(ielem(n, i)%q_face(l)%q_mapl(1), 1:nof_variables)
              else

                cright(1:nof_variables) = iexsolhir(ilocal_recon3(i)%ihexn(1, ielem(n, i)%indexi(l)))%sol &
                                          (ilocal_recon3(i)%ihexl(1, ielem(n, i)%indexi(l)), 1:nof_variables)

!                                                               cright(1:nof_variables)=iexboundhir(ielem(n,i)%ineighn(l))%facesol(ielem(n,i)%q_face(l)%q_mapl(1),1:nof_variables)

              end if
              kas = 4

              if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
                if (fastest .eq. 1) then
                  cturbr(1:turbulenceequations + passivescalar) = solchanger(ielem(n, i)%ineighn(l))%sol &
                                                        (ielem(n, i)%q_face(l)%q_mapl(1), 5:4 + turbulenceequations + passivescalar)
                else

                  cturbr(1:turbulenceequations + passivescalar) = iexsolhir(ilocal_recon3(i)%ihexn(1, ielem(n, i)%indexi(l)))%sol &
                                       (ilocal_recon3(i)%ihexl(1, ielem(n, i)%indexi(l)), 5:4 + turbulenceequations + passivescalar)
                end if
              end if

            end if
          else
            kas = 5
            if (fastest .eq. 1) then
              cright(1:nof_variables) = solchanger(ielem(n, i)%ineighn(l))%sol(ielem(n, i)%q_face(l)%q_mapl(1), 1:nof_variables)
            else

              cright(1:nof_variables) = iexsolhir(ilocal_recon3(i)%ihexn(1, ielem(n, i)%indexi(l)))%sol &
                                        (ilocal_recon3(i)%ihexl(1, ielem(n, i)%indexi(l)), 1:nof_variables)

!                                                               cright(1:nof_variables)=iexboundhir(ielem(n,i)%ineighn(l))%facesol(ielem(n,i)%q_face(l)%q_mapl(1),1:nof_variables)

            end if

            if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
              if (fastest .eq. 1) then
                cturbr(1:turbulenceequations + passivescalar) = solchanger(ielem(n, i)%ineighn(l))%sol &
                                                        (ielem(n, i)%q_face(l)%q_mapl(1), 5:4 + turbulenceequations + passivescalar)
              else

                cturbr(1:turbulenceequations + passivescalar) = iexsolhir(ilocal_recon3(i)%ihexn(1, ielem(n, i)%indexi(l)))%sol &
                                       (ilocal_recon3(i)%ihexl(1, ielem(n, i)%indexi(l)), 5:4 + turbulenceequations + passivescalar)
              end if
            end if

!
          end if
        end if

        call rotatef2d(n, cright_rot, cright, angle1, angle2)        !rotate wrt to normalvector of face and solve 1d riemann problem
        call rotatef2d(n, cleft_rot, cleft, angle1, angle2)        !rotate wrt to normalvector of face and solve 1d riemann problem

        if ((lmach .eq. 1)) then    !application of the low mach number correction
          leftv(1:nof_variables) = cleft_rot(1:nof_variables); rightv(1:nof_variables) = cright_rot(1:nof_variables)
          call lmacht2d(n, leftv, rightv)
          cleft_rot(1:nof_variables) = leftv(1:nof_variables); cright_rot(1:nof_variables) = rightv(1:nof_variables);
          call rotateb2d(n, cright, cright_rot, angle1, angle2)        !rotate wrt to normalvector of face and solve 1d riemann problem
          call rotateb2d(n, cleft, cleft_rot, angle1, angle2)

        end if

        leftv(1:nof_variables) = cleft(1:nof_variables); rightv(1:nof_variables) = cright(1:nof_variables)
        call cons2prim2(n, leftv, rightv, mp_pinfl, mp_pinfr, gammal, gammar)

        asound1 = sqrt(leftv(4)*gamma/leftv(1)) + abs(cleft_rot(2)/cleft_rot(1))
        asound2 = sqrt(rightv(4)*gamma/rightv(1)) + abs(cright_rot(2)/cright_rot(1))

        vpp = max(asound1, asound2)

        if (itestcase .eq. 4) then
          call sutherland2d(n, leftv, rightv, viscl, laml)

          viscots = (viscl(1) + viscl(2))*oo2

!                                                   viscots=viscots/((0.5*(cleft(1)+cright(1)))*ielem(n,i)%dih(l))
          viscots = max((4.0/(3.0*0.5*(cleft(1) + cright(1))))*viscots*mul1 &
                        /ielem(n, i)%dih(l), ((gamma/(0.5*(cleft(1) + cright(1)))) &
                                              *viscots*mul1/(ielem(n, i)%dih(l)*prandtl)))
          vpp = max(asound1, asound2) + viscots

          if (turbulence .eq. 1) then
            if (turbulencemodel .eq. 1) then
              turbmv(1) = cturbl(1); turbmv(2) = cturbr(1); eddyfl(2) = turbmv(1); eddyfr(2) = turbmv(2)
              call eddyvisco2d(n, viscl, laml, turbmv, etvm, eddyfl, eddyfr, leftv, rightv)
            end if
            if (turbulencemodel .eq. 2) then
              eddyfl(1) = ielem(n, i)%walldist; eddyfl(2) = cturbl(1); eddyfl(3) = cturbl(2)
              eddyfl(4:5) = ilocal_recon3(i)%grads(1, 1:2); eddyfl(6:7) = ilocal_recon3(i)%grads(2, 1:2)
              eddyfl(8:9) = ilocal_recon3(i)%grads(4, 1:2)
              eddyfl(10:11) = ilocal_recon3(i)%grads(5, 1:2)

              eddyfr = eddyfl
              call eddyvisco2d(n, viscl, laml, turbmv, etvm, eddyfl, eddyfr, leftv, rightv)
            end if

            viscots = oo2*((viscl(1) + viscl(3)) + (viscl(2) + viscl(4)))
!                                                       viscots=viscots/((0.5*(cleft(1)+cright(1)))*ielem(n,i)%dih(l))
            viscots = max((4.0/(3.0*0.5*(cleft(1) + cright(1))))*viscots* &
                          mul1/ielem(n, i)%dih(l), ((gamma/(0.5*(cleft(1) + cright(1))))* &
                                                    viscots*mul1/(ielem(n, i)%dih(l)*(prandtl + prtu))))

            vpp = max(asound1, asound2) + viscots
          end if

          if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
            if (turbulence .eq. 1) then
            do nvar = 1, turbulenceequations

              vpp = max(asound1, asound2) + viscots
!

            end do
            end if
            if (passivescalar .gt. 0) then
            do nvar = turbulenceequations + 1, turbulenceequations + passivescalar
              viscl(1) = viscl(1)/schmidt_lam
              viscl(2) = viscl(2)/schmidt_lam

              if (turbulence .eq. 1) then
                viscl(3) = viscl(3)/schmidt_turb
                viscl(4) = viscl(4)/schmidt_turb
                viscots = 0.5*((viscl(1) + viscl(3)) + (viscl(2) + viscl(4)))
              else
                viscots = 0.5*((viscl(1) + viscl(2)))

              end if

              viscots = (2.0*viscots)/((cleft(1) + cright(1))*ielem(n, i)%dih(l))
              vpp = max(asound1, asound2) + viscots

            end do
            end if
          end if
        else

        end if

        impdiag_mf(i) = impdiag_mf(i) + (oo2*vpp*mul1)
        impoff_mf(i, l) = vpp
        if (turbulence .eq. 1) then
          impdiagt(i, :) = impdiagt(i, :) + (oo2*vpp*mul1)

          impofft(i, l, :) = impofft(i, l, :) - (oo2*((vpp))*mul1)
        end if

      end do
    end do
    !$omp end do

    if (rungekutta .eq. 10) then
      !$omp do
      do i = 1, kmaxe

        impdiag_mf(i) = (impdiag_mf(i)) + (ielem(n, i)%totvolume/ielem(n, i)%dtl)

      end do
      !$omp end do
    else
      !$omp do
      do i = 1, kmaxe
        impdiag_mf(i) = ielem(n, i)%totvolume*((1.0d0/ielem(n, i)%dtl) + (1.5d0/dt)) + (impdiag_mf(i))
      end do
      !$omp end do
    end if

    if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then

      if (turbulence .eq. 1) call sources_derivatives_computation2d(n)
      if (rungekutta .eq. 10) then
!$omp do
        do i = 1, kmaxe
          if (turbulence .eq. 1) then
          do nvar = 1, turbulenceequations
            impdiagt(i, nvar) = impdiagt(i, nvar) + (ielem(n, i)%totvolume/ielem(n, i)%dtl) - sht(i, nvar)
          end do
          end if
          if (passivescalar .gt. 0) then
          do nvar = turbulenceequations + 1, turbulenceequations + passivescalar
            impdiagt(i, nvar) = impdiagt(i, nvar) + (ielem(n, i)%totvolume/ielem(n, i)%dtl) - sht(i, nvar)
          end do
          end if
        end do
!$omp end do
      else
!$omp do
        do i = 1, kmaxe
          if (turbulence .eq. 1) then
          do nvar = 1, turbulenceequations
!
            impdiagt(i, nvar) = ielem(n, i)%totvolume*((1.0d0/ielem(n, i)%dtl) + (1.5d0/dt)) + (impdiagt(i, nvar)) - sht(i, nvar)
          end do
          end if
          if (passivescalar .gt. 0) then
          do nvar = turbulenceequations + 1, turbulenceequations + passivescalar
            impdiagt(i, nvar) = impdiagt(i, nvar) + (ielem(n, i)%totvolume/ielem(n, i)%dtl) - sht(i, nvar)
          end do
          end if

        end do
!$omp end do

      end if
    end if

  end subroutine calculate_jacobian_2d_mf

  subroutine calculate_jacobian_3d_mf(n)
    !> @brief
!> this subroutine computes the approximate jacobian for implicit time stepping in 2d
    implicit none
    integer, intent(in)::n
    real, dimension(1:nof_variables + turbulenceequations + passivescalar)::godflux2
    integer::i, l, ngp, kmaxe, iqp, ii, nvar, n_node, ibfc, kas
    real::sum_detect, norms, vpp, asound1, asound2, mul1, dxb, tempxx, viscots
    real, dimension(nof_variables, nof_variables)::identity1
    real, dimension(nof_variables, nof_variables)::convj, diffj
    integer::iconsidered, facex, pointx, igoflux
    integer::b_code
    real::angle1, angle2, nx, ny, nz
    real, dimension(1:nof_variables)::cleft, cright, cright_rot, cleft_rot
    real, dimension(1:turbulenceequations + passivescalar)::cturbl, cturbr
    real, dimension(1:nof_variables)::leftv, srf_speedrot, srf_speed
    real, dimension(1:nof_variables)::rightv
    real, dimension(1:dimensiona)::pox, poy, poz
    real, dimension(1:8, 1:dimensiona)::vext, nodes_list

    real, dimension(1:4)::viscl, laml
    real, dimension(1:2)::turbmv
    real, dimension(1)::etvm
    real, dimension(1:20)::eddyfl, eddyfr

    real, dimension(1:dimensiona)::cords
    real::mp_pinfl, gammal
    real::mp_pinfr, gammar
    real, dimension(1:nof_variables, 1:nof_variables)::eigvl
    kmaxe = xmpielrank(n)

    !$omp do
    do ii = 1, nof_interior        !for all the interior elements
      i = el_int(ii)
      iconsidered = i
      impdiag_mf(i) = zero
      impoff_mf(i, :) = zero
      if (turbulence .eq. 1) then
        impdiagt(i, 1:turbulenceequations + passivescalar) = 0.0
        impofft(i, :, :) = zero
      end if
      do l = 1, ielem(n, i)%ifca !for all their faces
        b_code = 0
        angle1 = ielem(n, i)%faceanglex(l)
        angle2 = ielem(n, i)%faceangley(l)
        nx = angle1
        ny = angle2
        mul1 = ielem(n, i)%surf(l)

        cleft(1:nof_variables) = u_c(i)%val(1, 1:nof_variables)
        cright(1:nof_variables) = u_c(ielem(n, i)%ineigh(l))%val(1, 1:nof_variables)


        if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

          cturbl(1:turbulenceequations + passivescalar) = u_ct(i)%val(1, 1:turbulenceequations + passivescalar)!left additional equations flow state
          cturbr(1:turbulenceequations + passivescalar) = u_ct(ielem(n, i)%ineigh(l))%val(1, 1:turbulenceequations + passivescalar)

        end if

        call rotatef(n, cright_rot, cright, angle1, angle2)        !rotate wrt to normalvector of face and solve 1d riemann problem
        call rotatef(n, cleft_rot, cleft, angle1, angle2)        !rotate wrt to normalvector of face and solve 1d riemann problem

        if ((lmach .eq. 1)) then    !application of the low mach number correction
          leftv(1:nof_variables) = cleft_rot(1:nof_variables); rightv(1:nof_variables) = cright_rot(1:nof_variables)
          call lmacht(n, leftv, rightv)
          cleft_rot(1:nof_variables) = leftv(1:nof_variables); cright_rot(1:nof_variables) = rightv(1:nof_variables);
          call rotateb(n, cright, cright_rot, angle1, angle2)        !rotate wrt to normalvector of face and solve 1d riemann problem
          call rotateb(n, cleft, cleft_rot, angle1, angle2)

        end if

        leftv(1:nof_variables) = cleft(1:nof_variables); rightv(1:nof_variables) = cright(1:nof_variables)
        call cons2prim2(n, leftv, rightv, mp_pinfl, mp_pinfr, gammal, gammar)

        asound1 = sqrt(leftv(5)*gamma/leftv(1)) + abs(cleft_rot(2)/cleft_rot(1))
        asound2 = sqrt(rightv(5)*gamma/rightv(1)) + abs(cright_rot(2)/cright_rot(1))

        vpp = max(asound1, asound2)

        if (itestcase .eq. 4) then
          call sutherland(n, leftv, rightv, viscl, laml)

          viscots = (viscl(1) + viscl(2))*oo2

!                                                   viscots=viscots/((0.5*(cleft(1)+cright(1)))*ielem(n,i)%dih(l))
          viscots = max((4.0/(3.0*0.5*(cleft(1) + cright(1))))*viscots*mul1 &
                        /ielem(n, i)%dih(l), ((gamma/(0.5*(cleft(1) + cright(1)))) &
                                              *viscots*mul1/(ielem(n, i)%dih(l)*prandtl)))
          vpp = max(asound1, asound2) + viscots

          if (turbulence .eq. 1) then
            if (turbulencemodel .eq. 1) then
              turbmv(1) = cturbl(1); turbmv(2) = cturbr(1); eddyfl(2) = turbmv(1); eddyfr(2) = turbmv(2)
              call eddyvisco(n, viscl, laml, turbmv, etvm, eddyfl, eddyfr, leftv, rightv)
            end if
            if (turbulencemodel .eq. 2) then
              eddyfl(1) = ielem(n, i)%walldist; eddyfl(2) = cturbl(1); eddyfl(3) = cturbl(2)
              eddyfl(4:6) = ilocal_recon3(i)%grads(1, 1:3); eddyfl(7:9) = ilocal_recon3(i)%grads(2, 1:3)
              eddyfl(10:12) = ilocal_recon3(i)%grads(3, 1:3); eddyfl(13:15) = ilocal_recon3(i)%grads(5, 1:3)
              eddyfl(16:18) = ilocal_recon3(i)%grads(6, 1:3)

              eddyfr = eddyfl
              call eddyvisco(n, viscl, laml, turbmv, etvm, eddyfl, eddyfr, leftv, rightv)
            end if

            viscots = oo2*((viscl(1) + viscl(3)) + (viscl(2) + viscl(4)))

            !                                                       viscots=viscots/((0.5*(cleft(1)+cright(1)))*ielem(n,i)%dih(l))
            viscots = max((4.0/(3.0*0.5*(cleft(1) + cright(1))))*viscots* &
                          mul1/ielem(n, i)%dih(l), ((gamma/(0.5*(cleft(1) + cright(1))))* &
                                                    viscots*mul1/(ielem(n, i)%dih(l)*(prandtl + prtu))))

            vpp = max(asound1, asound2) + viscots
          end if
        end if

        impdiag_mf(i) = impdiag_mf(i) + (oo2*vpp*mul1)
        impoff_mf(i, l) = vpp

        if (turbulence .eq. 1) then
          impdiagt(i, :) = impdiagt(i, :) + (oo2*vpp*mul1)
          impofft(i, l, :) = impofft(i, l, :) - (oo2*((vpp))*mul1)
        end if

      end do
    end do
    !$omp end do

    !$omp do
    do ii = 1, nof_bounded
      i = el_bnd(ii)
      iconsidered = i

      impdiag_mf(i) = zero
      impoff_mf(i, :) = zero
      if (turbulence .eq. 1) then
        impdiagt(i, 1:turbulenceequations + passivescalar) = zero
        impofft(i, :, :) = zero
      end if
      do l = 1, ielem(n, i)%ifca
        mul1 = ielem(n, i)%surf(l)
        angle1 = ielem(n, i)%faceanglex(l)
        angle2 = ielem(n, i)%faceangley(l)
        nx = angle1
        ny = angle2

        b_code = 0
        cleft(1:nof_variables) = u_c(i)%val(1, 1:nof_variables)
!                                       cleft(1:nof_variables)=ilocal_recon3(i)%uleft(1:nof_variables,l,1)
        if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

          cturbl(1:turbulenceequations + passivescalar) = u_ct(i)%val(1, 1:turbulenceequations + passivescalar)

        end if

        if (ielem(n, i)%ineighb(l) .eq. n) then        !my cpu only
          if (ielem(n, i)%ibounds(l) .gt. 0) then        !check for boundaries
            if (ibound(n, ielem(n, i)%ibounds(l))%icode .eq. 5) then        !periodic in my cpu
              cright(1:nof_variables) = u_c(ielem(n, i)%ineigh(l))%val(1, 1:nof_variables)
!                                                                   cright(1:nof_variables)=ilocal_recon3(ielem(n,i)%ineigh(l))%uleft(1:nof_variables,ielem(n,i)%ineighn(l),1)

              if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

           cturbr(1:turbulenceequations + passivescalar) = u_ct(ielem(n, i)%ineigh(l))%val(1, 1:turbulenceequations + passivescalar)

              end if

              kas = 1

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

              kas = 2

            end if
          else
            cright(1:nof_variables) = u_c(ielem(n, i)%ineigh(l))%val(1, 1:nof_variables)
!                                                               cright(1:nof_variables)=ilocal_recon3(ielem(n,i)%ineigh(l))%uleft(1:nof_variables,ielem(n,i)%ineighn(l),1)

            if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then

           cturbr(1:turbulenceequations + passivescalar) = u_ct(ielem(n, i)%ineigh(l))%val(1, 1:turbulenceequations + passivescalar)

            end if

            kas = 3

          end if
        else        !in other cpus they can only be periodic or mpi neighbours

          if (ielem(n, i)%ibounds(l) .gt. 0) then        !check for boundaries
            if (ibound(n, ielem(n, i)%ibounds(l))%icode .eq. 5) then        !periodic in other cpu

              if (fastest .eq. 1) then
                cright(1:nof_variables) = solchanger(ielem(n, i)%ineighn(l))%sol(ielem(n, i)%q_face(l)%q_mapl(1), 1:nof_variables)
              else

                cright(1:nof_variables) = iexsolhir(ilocal_recon3(i)%ihexn(1, ielem(n, i)%indexi(l)))%sol &
                                          (ilocal_recon3(i)%ihexl(1, ielem(n, i)%indexi(l)), 1:nof_variables)

!                                                               cright(1:nof_variables)=iexboundhir(ielem(n,i)%ineighn(l))%facesol(ielem(n,i)%q_face(l)%q_mapl(1),1:nof_variables)

              end if
              kas = 4

              if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
                if (fastest .eq. 1) then
                  cturbr(1:turbulenceequations + passivescalar) = solchanger(ielem(n, i)%ineighn(l))%sol &
                                                        (ielem(n, i)%q_face(l)%q_mapl(1), 5:4 + turbulenceequations + passivescalar)
                else

                  cturbr(1:turbulenceequations + passivescalar) = iexsolhir(ilocal_recon3(i)%ihexn(1, ielem(n, i)%indexi(l)))%sol &
                                       (ilocal_recon3(i)%ihexl(1, ielem(n, i)%indexi(l)), 5:4 + turbulenceequations + passivescalar)
                end if
              end if

            end if
          else
            kas = 5
            if (fastest .eq. 1) then
              cright(1:nof_variables) = solchanger(ielem(n, i)%ineighn(l))%sol(ielem(n, i)%q_face(l)%q_mapl(1), 1:nof_variables)
            else

              cright(1:nof_variables) = iexsolhir(ilocal_recon3(i)%ihexn(1, ielem(n, i)%indexi(l)))%sol &
                                        (ilocal_recon3(i)%ihexl(1, ielem(n, i)%indexi(l)), 1:nof_variables)

            end if

            if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
              if (fastest .eq. 1) then
                cturbr(1:turbulenceequations + passivescalar) = solchanger(ielem(n, i)%ineighn(l))%sol &
                                                        (ielem(n, i)%q_face(l)%q_mapl(1), 5:4 + turbulenceequations + passivescalar)
              else

                cturbr(1:turbulenceequations + passivescalar) = iexsolhir(ilocal_recon3(i)%ihexn(1, ielem(n, i)%indexi(l)))%sol &
                                       (ilocal_recon3(i)%ihexl(1, ielem(n, i)%indexi(l)), 5:4 + turbulenceequations + passivescalar)
              end if
            end if

!
          end if
        end if

        call rotatef(n, cright_rot, cright, angle1, angle2)        !rotate wrt to normalvector of face and solve 1d riemann problem
        call rotatef(n, cleft_rot, cleft, angle1, angle2)        !rotate wrt to normalvector of face and solve 1d riemann problem

        if ((lmach .eq. 1)) then    !application of the low mach number correction
          leftv(1:nof_variables) = cleft_rot(1:nof_variables); rightv(1:nof_variables) = cright_rot(1:nof_variables)
          call lmacht(n, leftv, rightv)
          cleft_rot(1:nof_variables) = leftv(1:nof_variables); cright_rot(1:nof_variables) = rightv(1:nof_variables);
          call rotateb(n, cright, cright_rot, angle1, angle2)        !rotate wrt to normalvector of face and solve 1d riemann problem
          call rotateb(n, cleft, cleft_rot, angle1, angle2)

        end if

        leftv(1:nof_variables) = cleft(1:nof_variables); rightv(1:nof_variables) = cright(1:nof_variables)
        call cons2prim2(n, leftv, rightv, mp_pinfl, mp_pinfr, gammal, gammar)

        asound1 = sqrt(leftv(5)*gamma/leftv(1)) + abs(cleft_rot(2)/cleft_rot(1))
        asound2 = sqrt(rightv(5)*gamma/rightv(1)) + abs(cright_rot(2)/cright_rot(1))

        vpp = max(asound1, asound2)

        if (itestcase .eq. 4) then
          call sutherland(n, leftv, rightv, viscl, laml)

          viscots = (viscl(1) + viscl(2))*oo2

!                                                   viscots=viscots/((0.5*(cleft(1)+cright(1)))*ielem(n,i)%dih(l))
          viscots = max((4.0/(3.0*0.5*(cleft(1) + cright(1))))*viscots*mul1 &
                        /ielem(n, i)%dih(l), ((gamma/(0.5*(cleft(1) + cright(1)))) &
                                              *viscots*mul1/(ielem(n, i)%dih(l)*prandtl)))
          vpp = max(asound1, asound2) + viscots

          if (turbulence .eq. 1) then
            if (turbulencemodel .eq. 1) then
              turbmv(1) = cturbl(1); turbmv(2) = cturbr(1); eddyfl(2) = turbmv(1); eddyfr(2) = turbmv(2)
              call eddyvisco(n, viscl, laml, turbmv, etvm, eddyfl, eddyfr, leftv, rightv)
            end if
            if (turbulencemodel .eq. 2) then
              eddyfl(1) = ielem(n, i)%walldist; eddyfl(2) = cturbl(1); eddyfl(3) = cturbl(2)
              eddyfl(4:6) = ilocal_recon3(i)%grads(1, 1:3); eddyfl(7:9) = ilocal_recon3(i)%grads(2, 1:3)
              eddyfl(10:12) = ilocal_recon3(i)%grads(3, 1:3); eddyfl(13:15) = ilocal_recon3(i)%grads(5, 1:3)
              eddyfl(16:18) = ilocal_recon3(i)%grads(6, 1:3)

              eddyfr = eddyfl
              call eddyvisco(n, viscl, laml, turbmv, etvm, eddyfl, eddyfr, leftv, rightv)
            end if

            viscots = oo2*((viscl(1) + viscl(3)) + (viscl(2) + viscl(4)))
!                                                       viscots=viscots/((0.5*(cleft(1)+cright(1)))*ielem(n,i)%dih(l))
            viscots = max((4.0/(3.0*0.5*(cleft(1) + cright(1))))*viscots* &
                          mul1/ielem(n, i)%dih(l), ((gamma/(0.5*(cleft(1) + cright(1))))* &
                                                    viscots*mul1/(ielem(n, i)%dih(l)*(prandtl + prtu))))

            vpp = max(asound1, asound2) + viscots
          end if

          if ((turbulence .eq. 1) .or. (passivescalar .gt. 0)) then
            if (turbulence .eq. 1) then
            do nvar = 1, turbulenceequations

              vpp = max(asound1, asound2) + viscots
!

            end do
            end if
            if (passivescalar .gt. 0) then
            do nvar = turbulenceequations + 1, turbulenceequations + passivescalar
              viscl(1) = viscl(1)/schmidt_lam
              viscl(2) = viscl(2)/schmidt_lam

              if (turbulence .eq. 1) then
                viscl(3) = viscl(3)/schmidt_turb
                viscl(4) = viscl(4)/schmidt_turb
                viscots = 0.5*((viscl(1) + viscl(3)) + (viscl(2) + viscl(4)))
              else
                viscots = 0.5*((viscl(1) + viscl(2)))

              end if

              viscots = (2.0*viscots)/((cleft(1) + cright(1))*ielem(n, i)%dih(l))
              vpp = max(asound1, asound2) + viscots

            end do
            end if
          end if
        else

        end if

        impdiag_mf(i) = impdiag_mf(i) + (oo2*vpp*mul1)
        impoff_mf(i, l) = vpp
        if (turbulence .eq. 1) then
          impdiagt(i, :) = impdiagt(i, :) + (oo2*vpp*mul1)

          impofft(i, l, :) = impofft(i, l, :) - (oo2*((vpp))*mul1)
        end if

      end do
    end do
    !$omp end do

    if (rungekutta .eq. 10) then
      !$omp do
      do i = 1, kmaxe

        impdiag_mf(i) = (impdiag_mf(i)) + (ielem(n, i)%totvolume/ielem(n, i)%dtl)

      end do
      !$omp end do
    else
      !$omp do
      do i = 1, kmaxe
        impdiag_mf(i) = ielem(n, i)%totvolume*((1.0d0/ielem(n, i)%dtl) + (1.5d0/dt)) + (impdiag_mf(i))
      end do
      !$omp end do
    end if

    if ((turbulence .gt. 0) .or. (passivescalar .gt. 0)) then

      if (turbulence .eq. 1) call sources_derivatives_computation(n)
      if (rungekutta .eq. 10) then
!$omp do
        do i = 1, kmaxe
          if (turbulence .eq. 1) then
          do nvar = 1, turbulenceequations
            impdiagt(i, nvar) = impdiagt(i, nvar) + (ielem(n, i)%totvolume/ielem(n, i)%dtl) - sht(i, nvar)
          end do
          end if
          if (passivescalar .gt. 0) then
          do nvar = turbulenceequations + 1, turbulenceequations + passivescalar
            impdiagt(i, nvar) = impdiagt(i, nvar) + (ielem(n, i)%totvolume/ielem(n, i)%dtl) - sht(i, nvar)
          end do
          end if
        end do
!$omp end do
      else
!$omp do
        do i = 1, kmaxe
          if (turbulence .eq. 1) then
          do nvar = 1, turbulenceequations
!
            impdiagt(i, nvar) = ielem(n, i)%totvolume*((1.0d0/ielem(n, i)%dtl) + (1.5d0/dt)) + (impdiagt(i, nvar)) - sht(i, nvar)
          end do
          end if
          if (passivescalar .gt. 0) then
          do nvar = turbulenceequations + 1, turbulenceequations + passivescalar
            impdiagt(i, nvar) = impdiagt(i, nvar) + (ielem(n, i)%totvolume/ielem(n, i)%dtl) - sht(i, nvar)
          end do
          end if

        end do
!$omp end do

      end if
    end if

  end subroutine calculate_jacobian_3d_mf

end module implicit_fluxes
