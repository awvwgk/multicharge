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

module multicharge_coulomb_type
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use multicharge_blas, only : gemv, symv
   use multicharge_coulomb_cache, only : mchrg_coulomb_cache
   implicit none
   private

   public :: mchrg_coulomb_type

   type, abstract :: mchrg_coulomb_type
   contains
      procedure :: get_energy
      procedure :: get_gradient
      procedure(get_coulomb_matrix), deferred :: get_coulomb_matrix
      procedure(get_coulomb_derivs), deferred :: get_coulomb_derivs
   end type mchrg_coulomb_type


   abstract interface
      subroutine get_coulomb_matrix(self, mol, cache, amat)
         import :: mchrg_coulomb_type, mchrg_coulomb_cache, structure_type, wp
         class(mchrg_coulomb_type), intent(in) :: self
         type(structure_type), intent(in) :: mol
         type(mchrg_coulomb_cache), intent(in) :: cache
         real(wp), intent(inout) :: amat(:, :)
      end subroutine get_coulomb_matrix

      subroutine get_coulomb_derivs(self, mol, cache, qvec, dadr, dadL, atrace)
         import :: mchrg_coulomb_type, mchrg_coulomb_cache, structure_type, wp
         class(mchrg_coulomb_type), intent(in) :: self
         type(structure_type), intent(in) :: mol
         type(mchrg_coulomb_cache), intent(in) :: cache
         real(wp), intent(in) :: qvec(:)
         real(wp), intent(inout) :: dadr(:, :, :)
         real(wp), intent(inout) :: dadL(:, :, :)
         real(wp), intent(inout) :: atrace(:, :)
      end subroutine get_coulomb_derivs
   end interface


contains


subroutine get_energy(self, mol, cache, qvec, energy)
   class(mchrg_coulomb_type), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(mchrg_coulomb_cache), intent(in) :: cache
   real(wp), contiguous, intent(in) :: qvec(:)
   real(wp), intent(inout) :: energy

   real(wp), allocatable :: amat(:, :), vtmp(:)

   allocate(amat(mol%nat, mol%nat), vtmp(mol%nat))
   amat(:, :) = 0.0_wp

   call self%get_coulomb_matrix(mol, cache, amat)

   call symv(amat, qvec, vtmp)
   energy = energy + 0.5_wp * dot_product(qvec, vtmp)
end subroutine get_energy


subroutine get_gradient(self, mol, cache, qvec, gradient, sigma)
   class(mchrg_coulomb_type), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(mchrg_coulomb_cache), intent(in) :: cache
   real(wp), contiguous, intent(in) :: qvec(:)
   real(wp), contiguous, intent(inout) :: gradient(:, :)
   real(wp), contiguous, intent(inout) :: sigma(:, :)

   real(wp), allocatable :: dadr(:, :, :), dadL(:, :, :), atrace(:, :)

   allocate(dadr(3, mol%nat, mol%nat), dadL(3, 3, mol%nat), atrace(3, mol%nat))
   dadr(:, :, :) = 0.0_wp
   dadL(:, :, :) = 0.0_wp
   atrace(:, :) = 0.0_wp

   call self%get_coulomb_derivs(mol, cache, qvec, dadr, dadL, atrace)

   call gemv(dadr, qvec, gradient, beta=1.0_wp)
   call gemv(dadL, qvec, sigma, beta=1.0_wp, alpha=0.5_wp)
end subroutine get_gradient

end module multicharge_coulomb_type
