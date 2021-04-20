! This file is part of multicharge.
! SPDX-Identifier: Apache-2.0
!
! Licensed under the Apache License, Version 2.0 (the "License");
! you may not use this file except in compliance with the License.
! You may obtain a copy of the License at
!
!     http://www.apache.org/licenses/LICENSE-2.0
!
! Unless required by applicable law or agreed to in writing, software
! distributed under the License is distributed on an "AS IS" BASIS,
! WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
! See the License for the specific language governing permissions and
! limitations under the License.

module multicharge_model
   use mctc_env, only : error_type, wp
   use mctc_io, only : structure_type
   use mctc_io_constants, only : pi
   use mctc_io_math, only : matdet_3x3, matinv_3x3
   use multicharge_blas, only : gemv, symv, gemm
   use multicharge_coulomb, only : mchrg_coulomb_type, mchrg_coulomb_cache
   use multicharge_cutoff, only : get_lattice_points
   use multicharge_ewald, only : get_alpha
   use multicharge_lapack, only : sytrf, sytrs, sytri
   use multicharge_wignerseitz, only : wignerseitz_cell_type, new_wignerseitz_cell
   implicit none
   private

   public :: mchrg_model_type, new_mchrg_model


   type :: mchrg_model_type
      class(mchrg_coulomb_type), allocatable :: coulomb
      real(wp), allocatable :: chi(:)
      real(wp), allocatable :: eta(:)
      real(wp), allocatable :: kcn(:)
   contains
      procedure :: solve
   end type mchrg_model_type


contains


subroutine new_mchrg_model(self, chi, eta, kcn)
   type(mchrg_model_type), intent(out) :: self
   real(wp), intent(in) :: chi(:)
   real(wp), intent(in) :: eta(:)
   real(wp), intent(in) :: kcn(:)

   self%chi = chi
   self%eta = eta
   self%kcn = kcn

end subroutine new_mchrg_model


subroutine get_vrhs(self, mol, cn, xvec, dxdcn)
   type(mchrg_model_type), intent(in) :: self
   type(structure_type), intent(in) :: mol
   real(wp), intent(in) :: cn(:)
   real(wp), intent(out) :: xvec(:)
   real(wp), intent(out), optional :: dxdcn(:)
   real(wp), parameter :: reg = 1.0e-14_wp

   integer :: iat, izp
   real(wp) :: tmp

   if (present(dxdcn)) then
      !$omp parallel do default(none) schedule(runtime) &
      !$omp shared(mol, self, cn, xvec, dxdcn) private(iat, izp, tmp)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         tmp = self%kcn(izp) / sqrt(cn(iat) + reg)
         xvec(iat) = -self%chi(izp) + tmp*cn(iat)
         dxdcn(iat) = 0.5_wp*tmp
      end do
      dxdcn(mol%nat+1) = 0.0_wp
   else
      !$omp parallel do default(none) schedule(runtime) &
      !$omp shared(mol, self, cn, xvec) private(iat, izp, tmp)
      do iat = 1, mol%nat
         izp = mol%id(iat)
         tmp = self%kcn(izp) / sqrt(cn(iat) + reg)
         xvec(iat) = -self%chi(izp) + tmp*cn(iat)
      end do
   end if
   xvec(mol%nat+1) = mol%charge

end subroutine get_vrhs


subroutine solve(self, mol, cn, dcndr, dcndL, energy, gradient, sigma, qvec, dqdr, dqdL)
   class(mchrg_model_type), intent(in) :: self
   type(structure_type), intent(in) :: mol
   real(wp), intent(in), contiguous :: cn(:)
   real(wp), intent(in), contiguous, optional :: dcndr(:, :, :)
   real(wp), intent(in), contiguous, optional :: dcndL(:, :, :)
   real(wp), intent(out), contiguous, optional :: qvec(:)
   real(wp), intent(out), contiguous, optional :: dqdr(:, :, :)
   real(wp), intent(out), contiguous, optional :: dqdL(:, :, :)
   real(wp), intent(inout), contiguous, optional :: energy(:)
   real(wp), intent(inout), contiguous, optional :: gradient(:, :)
   real(wp), intent(inout), contiguous, optional :: sigma(:, :)

   integer :: ic, jc, iat, ndim, info
   logical :: grad, cpq, dcn
   integer, allocatable :: ipiv(:)
   real(wp), allocatable :: xvec(:), vrhs(:), amat(:, :), ainv(:, :)
   real(wp), allocatable :: dxdcn(:), atrace(:, :), dadr(:, :, :), dadL(:, :, :)
   type(mchrg_coulomb_cache) :: cache

   ndim = mol%nat + 1
   call cache%update(mol)

   dcn = present(dcndr) .and. present(dcndL)
   grad = present(gradient) .and. present(sigma) .and. dcn
   cpq = present(dqdr) .and. present(dqdL) .and. dcn

   allocate(amat(ndim, ndim), xvec(ndim))
   allocate(ipiv(ndim))
   if (grad.or.cpq) then
      allocate(dxdcn(ndim))
   end if

   call get_vrhs(self, mol, cn, xvec, dxdcn)

   amat(:, :) = 0.0_wp
   amat(ndim, :mol%nat) = 1.0_wp
   amat(:mol%nat, ndim) = 1.0_wp
   call self%coulomb%get_coulomb_matrix(mol, cache, amat)
   do iat = 1, mol%nat
      amat(iat, iat) = amat(iat, iat) + self%eta(mol%id(iat))
   end do

   vrhs = xvec
   ainv = amat

   call sytrf(ainv, ipiv, info=info, uplo='l')

   if (info == 0) then
      if (cpq) then
         call sytri(ainv, ipiv, info=info, uplo='l')
         if (info == 0) then
            call symv(ainv, xvec, vrhs, uplo='l')
            do ic = 1, ndim
               do jc = ic+1, ndim
                  ainv(ic, jc) = ainv(jc, ic)
               end do
            end do
         end if
      else
         call sytrs(ainv, vrhs, ipiv, info=info, uplo='l')
      end if
   end if

   if (present(qvec)) then
      qvec(:) = vrhs(:mol%nat)
   end if

   if (present(energy)) then
      call symv(amat, vrhs, xvec, alpha=0.5_wp, beta=-1.0_wp, uplo='l')
      energy(:) = energy(:) + vrhs(:mol%nat) * xvec(:mol%nat)
   end if

   if (grad.or.cpq) then
      allocate(dadr(3, mol%nat, ndim), dadL(3, 3, ndim), atrace(3, mol%nat))
      dadr(:, :, :) = 0.0_wp
      dadL(:, :, :) = 0.0_wp
      atrace(:, :) = 0.0_wp
      call self%coulomb%get_coulomb_derivs(mol, cache, vrhs, dadr, dadL, atrace)
      xvec(:) = -dxdcn * vrhs
   end if

   if (grad) then
      call gemv(dadr, vrhs, gradient, beta=1.0_wp)
      call gemv(dcndr, xvec(:mol%nat), gradient, beta=1.0_wp)
      call gemv(dadL, vrhs, sigma, beta=1.0_wp, alpha=0.5_wp)
      call gemv(dcndL, xvec(:mol%nat), sigma, beta=1.0_wp)
   end if

   if (cpq) then
      do iat = 1, mol%nat
         dadr(:, iat, iat) = atrace(:, iat) + dadr(:, iat, iat)
         dadr(:, :, iat) = -dcndr(:, :, iat) * dxdcn(iat) + dadr(:, :, iat)
         dadL(:, :, iat) = -dcndL(:, :, iat) * dxdcn(iat) + dadL(:, :, iat)
      end do

      call gemm(dadr, ainv(:, :mol%nat), dqdr, alpha=-1.0_wp)
      call gemm(dadL, ainv(:, :mol%nat), dqdL, alpha=-1.0_wp)
   end if

end subroutine solve


end module multicharge_model
