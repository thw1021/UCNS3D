module lapck
  use mpiinfo
  use declaration
  implicit none

contains

  real function lxnorm(xqr, pdim)
    implicit none
    integer, intent(in) :: pdim
    real, allocatable, dimension(:), intent(in)::xqr
    integer i
    lxnorm = 0.0d0
    do i = 1, pdim
      lxnorm = lxnorm + xqr(i)**2
    end do
    lxnorm = sqrt(lxnorm)
  end function

  subroutine house(xqr, vqr1, pdim)
    implicit none
    integer, intent(in) :: pdim
    real, allocatable, dimension(:), intent(inout)::xqr
    real, allocatable, dimension(:), intent(inout)::vqr1
    vqr1 = xqr
    vqr1(1) = xqr(1) + sign(1.0, xqr(1))*lxnorm(xqr, pdim)
  end subroutine

  subroutine computehousematrix(pqr, vqr, ideg)
    implicit none
    integer, intent(in) :: ideg
    real, allocatable, dimension(:, :), intent(inout)::pqr
    real, allocatable, dimension(:), intent(inout)::vqr
    real ::vnorm
    integer:: i, j
    pqr = 0.0d0
    do i = 1, ideg
      pqr(i, i) = 1.0d0
    end do
    vnorm = lxnorm(vqr, ideg)
    vqr = vqr/vnorm
    do i = 1, ideg
    do j = 1, ideg
      pqr(i, j) = pqr(i, j) - 2.0d0*vqr(i)*vqr(j)
    end do
    end do
  end subroutine

  subroutine transposematrix(qff, qtff, ideg)
    implicit none
    integer, intent(in) :: ideg
    real, allocatable, dimension(:, :), intent(in) :: qff
    real, allocatable, dimension(:, :), intent(inout) :: qtff
    integer i, j
    do i = 1, ideg
    do j = 1, ideg
      qtff(i, j) = qff(j, i)
    end do
    end do
  end subroutine

  subroutine transposematrix_dg(qff_dg, qtff_dg, ideg)
    implicit none
    integer, intent(in) :: ideg
    real, allocatable, dimension(:, :), intent(in) :: qff_dg
    real, allocatable, dimension(:, :), intent(inout) :: qtff_dg
    integer i, j
    do i = 1, ideg
    do j = 1, ideg
      qtff_dg(i, j) = qff_dg(j, i)
    end do
    end do
  end subroutine

  subroutine qrdecomposition(lscqm, qff, rff, ideg)
    implicit none
    integer, intent(in) :: ideg
    real, allocatable, dimension(:, :), intent(in) :: lscqm
    real, allocatable, dimension(:, :), intent(inout) ::  qff, rff
    real::identity(ideg, ideg)
    real:: test(ideg)
    integer:: i, j, l, pdim
    real, allocatable, dimension(:)::xqr
    real, allocatable, dimension(:, :)::pqr
    real, allocatable, dimension(:)::vqr
    real, allocatable, dimension(:)::vqr1
    identity(1:ideg, 1:ideg) = 0.0d0
    qff(1:ideg, 1:ideg) = 0.0d0
    rff(1:ideg, 1:ideg) = 0.0d0
    do i = 1, ideg
      identity(i, i) = 1.0d0
    end do
    qff(1:ideg, 1:ideg) = identity(1:ideg, 1:ideg)
    rff(1:ideg, 1:ideg) = lscqm(1:ideg, 1:ideg)
    allocate (pqr(1:ideg, 1:ideg), vqr(1:ideg)); pqr = 0.0d0; vqr = 0.0d0
    do l = 1, ideg
      pdim = ideg - l + 1
      allocate (vqr1(1:pdim), xqr(1:pdim))
      xqr = rff(l:ideg, l)
      call house(xqr, vqr1, pdim)
      vqr(1:ideg) = 0.0d0
      vqr(l:ideg) = vqr1
      call computehousematrix(pqr, vqr, ideg)
      rff(1:ideg, 1:ideg) = matmul(pqr(1:ideg, 1:ideg), rff(1:ideg, 1:ideg))
      qff(1:ideg, 1:ideg) = matmul(qff(1:ideg, 1:ideg), pqr(1:ideg, 1:ideg))
      deallocate (vqr1, xqr)
    end do
    deallocate (pqr, vqr)
  end subroutine
  subroutine qrdecomposition_dg(lscqm_dg, qff_dg, rff_dg, ideg)
    implicit none
    integer, intent(in) :: ideg
    real, allocatable, dimension(:, :), intent(in) :: lscqm_dg
    real, allocatable, dimension(:, :), intent(inout) ::  qff_dg, rff_dg
    real::identity(ideg, ideg)
    real:: test(ideg)
    integer:: i, j, l, pdim
    real, allocatable, dimension(:)::xqr
    real, allocatable, dimension(:, :)::pqr
    real, allocatable, dimension(:)::vqr
    real, allocatable, dimension(:)::vqr1
    identity(1:ideg, 1:ideg) = 0.0d0
    qff_dg(1:ideg, 1:ideg) = 0.0d0
    rff_dg(1:ideg, 1:ideg) = 0.0d0
    do i = 1, ideg
      identity(i, i) = 1.0d0
    end do
    qff_dg(1:ideg, 1:ideg) = identity(1:ideg, 1:ideg)
    rff_dg(1:ideg, 1:ideg) = lscqm_dg(1:ideg, 1:ideg)
    allocate (pqr(1:ideg, 1:ideg), vqr(1:ideg)); pqr = 0.0d0; vqr = 0.0d0
    do l = 1, ideg
      pdim = ideg - l + 1
      allocate (vqr1(1:pdim), xqr(1:pdim))
      xqr = rff_dg(l:ideg, l)
      call house(xqr, vqr1, pdim)
      vqr(1:ideg) = 0.0d0
      vqr(l:ideg) = vqr1
      call computehousematrix(pqr, vqr, ideg)
      rff_dg(1:ideg, 1:ideg) = matmul(pqr(1:ideg, 1:ideg), rff_dg(1:ideg, 1:ideg))
      qff_dg(1:ideg, 1:ideg) = matmul(qff_dg(1:ideg, 1:ideg), pqr(1:ideg, 1:ideg))
      deallocate (vqr1, xqr)
    end do
    deallocate (pqr, vqr)
  end subroutine

end module
