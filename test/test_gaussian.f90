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

module test_gaussian
   use mctc_env, only : wp
   use mctc_env_testing, only : new_unittest, unittest_type, error_type, check, test_failed
   use mctc_io_structure, only : structure_type
   use mstore, only : get_structure
   use multicharge_blas, only : gemv, symv
   use multicharge_coulomb, only : mchrg_coulomb_cache
   use multicharge_coulomb_gaussian, only : mchrg_gaussian_coulomb, new_gaussian_coulomb
   implicit none
   private

   public :: collect_gaussian

   real(wp), parameter :: thr = 100*epsilon(1.0_wp)
   real(wp), parameter :: thr2 = sqrt(epsilon(1.0_wp))

   abstract interface
      subroutine coulomb_maker(coulomb, mol)
         import :: mchrg_gaussian_coulomb, structure_type
         type(mchrg_gaussian_coulomb), intent(out) :: coulomb
         type(structure_type), intent(in) :: mol
      end subroutine coulomb_maker
   end interface

contains


!> Collect all exported unit tests
subroutine collect_gaussian(testsuite)

   !> Collection of tests
   type(unittest_type), allocatable, intent(out) :: testsuite(:)

   testsuite = [ &
      new_unittest("energy-1", test_e_effective_m01), &
      new_unittest("energy-2", test_e_effective_m02), &
      new_unittest("energy-pbc", test_e_effective_oxacb), &
      new_unittest("gradient-1", test_g_effective_m03), &
      new_unittest("gradient-2", test_g_effective_m04), &
      new_unittest("gradient-pbc", test_g_effective_co2), &
      new_unittest("sigma-1", test_s_effective_m05), &
      new_unittest("sigma-2", test_s_effective_m06), &
      new_unittest("sigma-pbc", test_s_effective_ammonia) &
      ]

end subroutine collect_gaussian


!> Factory to create electrostatic objects based on DFT-D4 values
subroutine make_coulomb4(coulomb, mol)
   use multicharge_param_eeq2019, only : get_eeq_rad

   !> New electrostatic object
   type(mchrg_gaussian_coulomb), intent(out) :: coulomb

   !> Molecular structure data
   type(structure_type), intent(in) :: mol

   real(wp), allocatable :: rad(:)

   rad = get_eeq_rad(mol%num)
   call new_gaussian_coulomb(coulomb, rad)

end subroutine make_coulomb4


subroutine test_generic(error, mol, charges, make_coulomb, ref)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Atomic partial charges for this structure
   real(wp), contiguous, intent(in) :: charges(:)

   !> Factory to create new electrostatic objects
   procedure(coulomb_maker) :: make_coulomb

   !> Reference value to check against
   real(wp), intent(in) :: ref

   type(mchrg_gaussian_coulomb) :: coulomb
   type(mchrg_coulomb_cache) :: cache
   real(wp) :: energy

   energy = 0.0_wp
   call make_coulomb(coulomb, mol)
   call cache%update(mol)
   call coulomb%get_energy(mol, cache, charges, energy)

   call check(error, energy, ref, thr=thr)
   if (allocated(error)) then
      print*,energy
   end if

end subroutine test_generic


subroutine test_numgrad(error, mol, charges, make_coulomb)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Atomic partial charges for this structure
   real(wp), contiguous, intent(in) :: charges(:)

   !> Factory to create new electrostatic objects
   procedure(coulomb_maker) :: make_coulomb

   integer :: iat, ic
   type(mchrg_gaussian_coulomb) :: coulomb
   type(mchrg_coulomb_cache) :: cache
   real(wp) :: er, el, sigma(3, 3)
   real(wp), allocatable :: gradient(:, :), numgrad(:, :)
   real(wp), parameter :: step = 1.0e-6_wp

   allocate(gradient(3, mol%nat), numgrad(3, mol%nat))
   gradient(:, :) = 0.0_wp
   sigma(:, :) = 0.0_wp
   call make_coulomb(coulomb, mol)

   do iat = 1, mol%nat
      do ic = 1, 3
         er = 0.0_wp
         el = 0.0_wp
         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         call cache%update(mol)
         call coulomb%get_energy(mol, cache, charges, er)

         mol%xyz(ic, iat) = mol%xyz(ic, iat) - 2*step
         call cache%update(mol)
         call coulomb%get_energy(mol, cache, charges, el)

         mol%xyz(ic, iat) = mol%xyz(ic, iat) + step
         numgrad(ic, iat) = 0.5_wp*(er - el)/step
      end do
   end do

   call cache%update(mol)
   call coulomb%get_gradient(mol, cache, charges, gradient, sigma)

   if (any(abs(gradient - numgrad) > thr2)) then
      call test_failed(error, "Gradient of electrostatic energy does not match")
      print'(3es21.14)', gradient-numgrad
   end if

end subroutine test_numgrad


subroutine test_numsigma(error, mol, charges, make_coulomb)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   !> Molecular structure data
   type(structure_type), intent(inout) :: mol

   !> Atomic partial charges for this structure
   real(wp), contiguous, intent(in) :: charges(:)

   !> Factory to create new electrostatic objects
   procedure(coulomb_maker) :: make_coulomb

   integer :: ic, jc
   type(mchrg_gaussian_coulomb) :: coulomb
   type(mchrg_coulomb_cache) :: cache
   real(wp) :: er, el, sigma(3, 3), eps(3, 3), numsigma(3, 3)
   real(wp), allocatable :: gradient(:, :), xyz(:, :), lattice(:, :)
   real(wp), parameter :: unity(3, 3) = reshape(&
      & [1, 0, 0, 0, 1, 0, 0, 0, 1], shape(unity))
   real(wp), parameter :: step = 1.0e-6_wp

   allocate(gradient(3, mol%nat), xyz(3, mol%nat))
   gradient(:, :) = 0.0_wp
   sigma(:, :) = 0.0_wp
   call make_coulomb(coulomb, mol)

   eps(:, :) = unity
   xyz(:, :) = mol%xyz
   if (any(mol%periodic)) lattice = mol%lattice
   do ic = 1, 3
      do jc = 1, 3
         er = 0.0_wp
         el = 0.0_wp
         eps(jc, ic) = eps(jc, ic) + step
         mol%xyz(:, :) = matmul(eps, xyz)
         if (allocated(lattice)) mol%lattice(:, :) = matmul(eps, lattice)
         call cache%update(mol)
         call coulomb%get_energy(mol, cache, charges, er)
         eps(jc, ic) = eps(jc, ic) - 2*step
         mol%xyz(:, :) = matmul(eps, xyz)
         if (allocated(lattice)) mol%lattice(:, :) = matmul(eps, lattice)
         call cache%update(mol)
         call coulomb%get_energy(mol, cache, charges, el)
         eps(jc, ic) = eps(jc, ic) + step
         mol%xyz(:, :) = xyz
         if (allocated(lattice)) mol%lattice = lattice
         numsigma(jc, ic) = 0.5_wp*(er - el)/step
      end do
   end do

   call cache%update(mol)
   call coulomb%get_gradient(mol, cache, charges, gradient, sigma)

   if (any(abs(sigma - numsigma) > thr2)) then
      call test_failed(error, "Strain derivatives do not match")
      print'(3es21.14)', sigma-numsigma
   end if

end subroutine test_numsigma


subroutine test_e_effective_m01(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: charges(*) = [&
      & 7.73347900345264E-1_wp, 1.07626888948184E-1_wp,-3.66999593831010E-1_wp,&
      & 4.92833325937897E-2_wp,-1.83332156197733E-1_wp, 2.33302086605469E-1_wp,&
      & 6.61837152062315E-2_wp,-5.43944165050002E-1_wp,-2.70264356583716E-1_wp,&
      & 2.66618968841682E-1_wp, 2.62725033202480E-1_wp,-7.15315510172571E-2_wp,&
      &-3.73300777019193E-1_wp, 3.84585237785621E-2_wp,-5.05851088366940E-1_wp,&
      & 5.17677238544189E-1_wp]

   call get_structure(mol, "MB16-43", "01")
   call test_generic(error, mol, charges, make_coulomb4, 0.29610787220648987_wp)

end subroutine test_e_effective_m01


subroutine test_e_effective_m02(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: charges(*) = [&
      & 7.38394711236234E-2_wp,-1.68354976558608E-1_wp,-3.47642833746823E-1_wp,&
      &-7.05489267186003E-1_wp, 7.73548301641266E-1_wp, 2.30207581365386E-1_wp,&
      & 1.02748501676354E-1_wp, 9.47818107467040E-2_wp, 2.44260351729187E-2_wp,&
      & 2.34984927037408E-1_wp,-3.17839896393030E-1_wp, 6.67112994818879E-1_wp,&
      &-4.78119977010488E-1_wp, 6.57536027459275E-2_wp, 1.08259054549882E-1_wp,&
      &-3.58215329983396E-1_wp]

   call get_structure(mol, "MB16-43", "02")
   call test_generic(error, mol, charges, make_coulomb4, 0.32591359741119319_wp)

end subroutine test_e_effective_m02


subroutine test_e_effective_oxacb(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: charges(*) = [&
      & 3.41731844312030E-1_wp, 3.41716020106239E-1_wp, 3.41730526585671E-1_wp,&
      & 3.41714427217954E-1_wp, 3.80996046757999E-1_wp, 3.80989821246195E-1_wp,&
      & 3.81000747720282E-1_wp, 3.80990494183703E-1_wp,-3.70406587264474E-1_wp,&
      &-3.70407565207006E-1_wp,-3.70417590212352E-1_wp,-3.70399716470705E-1_wp,&
      &-3.52322260586075E-1_wp,-3.52304269439196E-1_wp,-3.52313440903261E-1_wp,&
      &-3.52298498047004E-1_wp]

   call get_structure(mol, "X23", "oxacb")
   call test_generic(error, mol, charges, make_coulomb4, 0.27372827485350670_wp)

end subroutine test_e_effective_oxacb


subroutine test_g_effective_m03(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: charges(*) = [&
      &-1.77788256288236E-1_wp,-8.22943267808161E-1_wp, 4.04578389873281E-2_wp,&
      & 5.79710531992282E-1_wp, 6.99601887637659E-1_wp, 6.84309612639107E-2_wp,&
      &-3.42971414989811E-1_wp, 4.64954031865410E-2_wp, 6.77012204116428E-2_wp,&
      & 8.49931225363225E-2_wp,-5.22285304699699E-1_wp,-2.92515001764712E-1_wp,&
      &-3.98375452377043E-1_wp, 2.09769668297792E-1_wp, 7.23140464830357E-1_wp,&
      & 3.65775987838250E-2_wp]

   call get_structure(mol, "MB16-43", "03")
   call test_numgrad(error, mol, charges, make_coulomb4)

end subroutine test_g_effective_m03


subroutine test_g_effective_m04(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: charges(*) = [&
      & 9.33596160193497E-2_wp,-3.41088061922851E-1_wp, 7.32474961830646E-2_wp,&
      &-2.21649975471802E-1_wp, 6.24413528413759E-3_wp, 1.07366683260668E-1_wp,&
      & 1.25982547197317E-1_wp, 9.65935501843890E-2_wp, 1.02704543049803E-1_wp,&
      & 1.45380937882263E-1_wp,-1.55978251071729E-1_wp, 3.42948437914661E-1_wp,&
      & 5.65504846503244E-2_wp,-3.37789986050220E-1_wp, 1.13510089629769E-1_wp,&
      &-2.07382246739143E-1_wp]

   call get_structure(mol, "MB16-43", "04")
   call test_numgrad(error, mol, charges, make_coulomb4)

end subroutine test_g_effective_m04


subroutine test_g_effective_co2(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: charges(*) = [&
      & 4.56275672862067E-1_wp, 4.56284770386671E-1_wp, 4.56284770386671E-1_wp,&
      & 4.56284770386671E-1_wp,-2.28127680925611E-1_wp,-2.28138283131909E-1_wp,&
      &-2.28145770512561E-1_wp,-2.28145770512561E-1_wp,-2.28150142163058E-1_wp,&
      &-2.28145770512561E-1_wp,-2.28138283131909E-1_wp,-2.28138283131909E-1_wp]

   call get_structure(mol, "X23", "CO2")
   call test_numgrad(error, mol, charges, make_coulomb4)

end subroutine test_g_effective_co2


subroutine test_s_effective_m05(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: charges(*) = [&
      &-2.01138111283277E-1_wp, 1.30358706339300E-1_wp, 9.38825924720944E-2_wp,&
      & 8.92795900801844E-2_wp, 5.13625440660610E-2_wp,-2.65500121876709E-2_wp,&
      & 9.26496972837658E-2_wp,-9.61095258223972E-2_wp,-4.92845009674246E-1_wp,&
      & 2.66730531684206E-1_wp, 3.37256104303071E-2_wp, 1.63170419985976E-1_wp,&
      & 6.91343155032824E-2_wp, 1.04287482572171E-1_wp, 6.09307909835941E-2_wp,&
      &-3.38869622433350E-1_wp]

   call get_structure(mol, "MB16-43", "05")
   call test_numsigma(error, mol, charges, make_coulomb4)

end subroutine test_s_effective_m05


subroutine test_s_effective_m06(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: charges(*) = [&
      &-2.13983049532933E-1_wp,-5.10521279217923E-1_wp, 7.70190120699491E-2_wp,&
      &-3.68835155548212E-1_wp,-4.08747874260092E-1_wp,-4.09471309598929E-2_wp,&
      & 2.94164204769172E-1_wp, 9.76819709672870E-2_wp,-7.84337476935767E-3_wp,&
      & 7.07702520795024E-1_wp, 2.38774840136381E-1_wp, 1.08934666297455E-1_wp,&
      & 1.10156911889136E-1_wp, 9.25098455002779E-2_wp,-1.96776817442259E-1_wp,&
      & 2.07107093059868E-2_wp]

   call get_structure(mol, "MB16-43", "06")
   call test_numsigma(error, mol, charges, make_coulomb4)

end subroutine test_s_effective_m06


subroutine test_s_effective_ammonia(error)

   !> Error handling
   type(error_type), allocatable, intent(out) :: error

   type(structure_type) :: mol
   real(wp), parameter :: charges(*) = [&
      & 2.95376975876519E-1_wp, 2.95376975876519E-1_wp, 2.95376975876519E-1_wp,&
      & 2.95329109335847E-1_wp, 2.95332441468412E-1_wp, 2.95347202855778E-1_wp,&
      & 2.95347202855779E-1_wp, 2.95329109335848E-1_wp, 2.95332441468411E-1_wp,&
      & 2.95347202855777E-1_wp, 2.95329109335847E-1_wp, 2.95332441468412E-1_wp,&
      &-8.86118742099358E-1_wp,-8.86012815503436E-1_wp,-8.86012815503437E-1_wp,&
      &-8.86012815503434E-1_wp]

   call get_structure(mol, "X23", "ammonia")
   call test_numsigma(error, mol, charges, make_coulomb4)

end subroutine test_s_effective_ammonia


end module test_gaussian
