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

module multicharge_param
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use multicharge_model, only : mchrg_model_type, new_mchrg_model
   use multicharge_coulomb, only : mchrg_coulomb_type
   use multicharge_coulomb_gaussian, only : mchrg_gaussian_coulomb, new_gaussian_coulomb
   use multicharge_param_eeq2019, only : get_eeq_chi, get_eeq_eta, &
      & get_eeq_rad, get_eeq_kcn
   implicit none
   private

   public :: new_eeq2019_model

contains

subroutine new_eeq2019_model(mol, model)
   !> Molecular structure data
   type(structure_type), intent(in) :: mol
   !> Electronegativity equilibration model
   type(mchrg_model_type), intent(out) :: model

   real(wp), allocatable :: chi(:), eta(:), kcn(:), rad(:)
   type(mchrg_gaussian_coulomb), allocatable :: coulomb

   chi = get_eeq_chi(mol%num)
   eta = get_eeq_eta(mol%num)
   kcn = get_eeq_kcn(mol%num)
   rad = get_eeq_rad(mol%num)

   call new_mchrg_model(model, chi=chi, eta=eta, kcn=kcn)

   allocate(coulomb)
   call new_gaussian_coulomb(coulomb, rad=rad)
   call move_alloc(coulomb, model%coulomb)

end subroutine new_eeq2019_model

end module multicharge_param
