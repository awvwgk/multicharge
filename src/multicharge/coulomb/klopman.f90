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

module multicharge_coulomb_klopman
   use mctc_env, only : wp
   use mctc_io, only : structure_type
   use mctc_io_constants, only : pi
   use mctc_io_math, only : matdet_3x3, matinv_3x3
   use multicharge_coulomb_cache, only : mchrg_coulomb_cache
   use multicharge_coulomb_type, only : mchrg_coulomb_type
   use multicharge_cutoff, only : get_lattice_points
   use multicharge_wignerseitz, only : wignerseitz_cell_type
   implicit none
   private

   public :: mchrg_klopman_coulomb, new_klopman_coulomb
   public :: average_interface
   public :: harmonic_average, arithmetic_average, geometric_average

   type, extends(mchrg_coulomb_type) :: mchrg_klopman_coulomb
      real(wp) :: gexp
      real(wp), allocatable :: hardness(:, :)
   contains
      procedure :: get_coulomb_matrix
      procedure :: get_coulomb_derivs
   end type mchrg_klopman_coulomb


   abstract interface
      pure function average_interface(gi, gj) result(gij)
         import :: wp
         !> Hardness of shell i
         real(wp), intent(in) :: gi
         !> Hardness of shell j
         real(wp), intent(in) :: gj
         !> Averaged hardness
         real(wp) :: gij
      end function average_interface
   end interface

   real(wp), parameter :: twopi = 2 * pi
   real(wp), parameter :: sqrtpi = sqrt(pi)
   real(wp), parameter :: eps = sqrt(epsilon(0.0_wp))

contains

subroutine new_klopman_coulomb(self, gexp, hardness, average)
   type(mchrg_klopman_coulomb), intent(out) :: self
   real(wp), intent(in) :: gexp
   procedure(average_interface) :: average
   real(wp), intent(in) :: hardness(:)

   integer :: isp, jsp

   self%gexp = gexp

   allocate(self%hardness(size(hardness), size(hardness)))
   do isp = 1, size(hardness)
      do jsp = 1, size(hardness)
         self%hardness(jsp, isp) = average(hardness(jsp), hardness(isp))
      end do
   end do
end subroutine new_klopman_coulomb


!> Harmonic averaging functions for hardnesses in GFN1-xTB
pure function harmonic_average(gi, gj) result(gij)
   !> Hardness of shell i
   real(wp), intent(in) :: gi
   !> Hardness of shell j
   real(wp), intent(in) :: gj
   !> Averaged hardness
   real(wp) :: gij

   gij = 2.0_wp/(1.0_wp/gi+1.0_wp/gj)

end function harmonic_average

!> Arithmetic averaging functions for hardnesses in GFN2-xTB
pure function arithmetic_average(gi, gj) result(gij)
   !> Hardness of shell i
   real(wp), intent(in) :: gi
   !> Hardness of shell j
   real(wp), intent(in) :: gj
   !> Averaged hardness
   real(wp) :: gij

   gij = 0.5_wp*(gi+gj)

end function arithmetic_average

!> Geometric averaging functions for hardnesses
pure function geometric_average(gi, gj) result(gij)
   !> Hardness of shell i
   real(wp), intent(in) :: gi
   !> Hardness of shell j
   real(wp), intent(in) :: gj
   !> Averaged hardness
   real(wp) :: gij

   gij = sqrt(gi*gj)

end function geometric_average


subroutine get_coulomb_matrix(self, mol, cache, amat)
   class(mchrg_klopman_coulomb), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(mchrg_coulomb_cache), intent(in) :: cache
   real(wp), intent(inout) :: amat(:, :)

   if (any(mol%periodic)) then
      call get_amat_3d(mol, cache%wsc, cache%alpha, self%gexp, self%hardness, amat)
   else
      call get_amat_0d(mol, self%gexp, self%hardness, amat)
   end if
end subroutine get_coulomb_matrix

subroutine get_dir_trans(lattice, trans)
   real(wp), intent(in) :: lattice(:, :)
   real(wp), allocatable, intent(out) :: trans(:, :)
   integer, parameter :: rep(3) = 2

   call get_lattice_points(lattice, rep, .true., trans)

end subroutine get_dir_trans

subroutine get_rec_trans(lattice, trans)
   real(wp), intent(in) :: lattice(:, :)
   real(wp), allocatable, intent(out) :: trans(:, :)
   integer, parameter :: rep(3) = 2
   real(wp) :: rec_lat(3, 3)

   rec_lat = twopi*transpose(matinv_3x3(lattice))
   call get_lattice_points(rec_lat, rep, .false., trans)

end subroutine get_rec_trans


subroutine get_amat_0d(mol, gexp, hardness, amat)
   type(structure_type), intent(in) :: mol
   real(wp), intent(in) :: gexp
   real(wp), intent(in) :: hardness(:, :)
   real(wp), intent(inout) :: amat(:, :)

   integer :: iat, jat, izp, jzp
   real(wp) :: vec(3), r1, r1g, gam, tmp

   !$omp parallel do default(none) schedule(runtime) &
   !$omp reduction(+:amat) shared(mol, gexp, hardness) &
   !$omp private(iat, izp, jat, jzp, gam, vec, r1, r1g, tmp)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      do jat = 1, iat-1
         jzp = mol%id(jat)
         vec = mol%xyz(:, jat) - mol%xyz(:, iat)
         r1 = norm2(vec)
         r1g = r1**gexp
         gam = hardness(jzp, izp)
         tmp = 1.0_wp/(r1g + gam**(-gexp))**(1.0_wp/gexp)
         amat(jat, iat) = amat(jat, iat) + tmp
         amat(iat, jat) = amat(iat, jat) + tmp
      end do
      tmp = hardness(izp, izp)
      amat(iat, iat) = amat(iat, iat) + tmp
   end do

end subroutine get_amat_0d

subroutine get_amat_3d(mol, wsc, alpha, gexp, hardness, amat)
   type(structure_type), intent(in) :: mol
   type(wignerseitz_cell_type), intent(in) :: wsc
   real(wp), intent(in) :: alpha
   real(wp), intent(in) :: gexp
   real(wp), intent(in) :: hardness(:, :)
   real(wp), intent(inout) :: amat(:, :)

   integer :: iat, jat, izp, jzp, img
   real(wp) :: vec(3), gam, wsw, dtmp, rtmp, vol
   real(wp), allocatable :: dtrans(:, :), rtrans(:, :)

   vol = abs(matdet_3x3(mol%lattice))
   call get_dir_trans(mol%lattice, dtrans)
   call get_rec_trans(mol%lattice, rtrans)

   !$omp parallel do default(none) schedule(runtime) &
   !$omp reduction(+:amat) shared(mol, gexp, hardness, wsc, dtrans, rtrans, alpha, vol) &
   !$omp private(iat, izp, jat, jzp, gam, wsw, vec, dtmp, rtmp)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      do jat = 1, iat-1
         jzp = mol%id(jat)
         gam = hardness(jzp, izp)
         wsw = 1.0_wp / real(wsc%nimg(jat, iat), wp)
         do img = 1, wsc%nimg(jat, iat)
            vec = mol%xyz(:, iat) - mol%xyz(:, jat) - wsc%trans(:, wsc%tridx(img, jat, iat))
            call get_amat_dir_3d(vec, gexp, gam, alpha, dtrans, dtmp)
            call get_amat_rec_3d(vec, vol, alpha, rtrans, rtmp)
            amat(jat, iat) = amat(jat, iat) + (dtmp + rtmp) * wsw
            amat(iat, jat) = amat(iat, jat) + (dtmp + rtmp) * wsw
         end do
      end do

      gam = hardness(izp, izp)
      wsw = 1.0_wp / real(wsc%nimg(iat, iat), wp)
      do img = 1, wsc%nimg(iat, iat)
         vec = wsc%trans(:, wsc%tridx(img, iat, iat))
         call get_amat_dir_3d(vec, gexp, gam, alpha, dtrans, dtmp)
         call get_amat_rec_3d(vec, vol, alpha, rtrans, rtmp)
         amat(iat, iat) = amat(iat, iat) + (dtmp + rtmp) * wsw
      end do

      dtmp = gam - 2 * alpha / sqrtpi
      amat(iat, iat) = amat(iat, iat) + dtmp
   end do

end subroutine get_amat_3d

subroutine get_amat_dir_3d(rij, gam, gexp, alp, trans, amat)
   real(wp), intent(in) :: rij(3)
   real(wp), intent(in) :: gam
   real(wp), intent(in) :: gexp
   real(wp), intent(in) :: alp
   real(wp), intent(in) :: trans(:, :)
   real(wp), intent(out) :: amat

   integer :: itr
   real(wp) :: vec(3), r1, tmp

   amat = 0.0_wp

   do itr = 1, size(trans, 2)
      vec(:) = rij + trans(:, itr)
      r1 = norm2(vec)
      if (r1 < eps) cycle
      tmp = 1.0_wp/(r1**gexp + gam**(-gexp))**(1.0_wp/gexp) - erf(alp*r1)/r1
      amat = amat + tmp
   end do

end subroutine get_amat_dir_3d

subroutine get_amat_rec_3d(rij, vol, alp, trans, amat)
   real(wp), intent(in) :: rij(3)
   real(wp), intent(in) :: vol
   real(wp), intent(in) :: alp
   real(wp), intent(in) :: trans(:, :)
   real(wp), intent(out) :: amat

   integer :: itr
   real(wp) :: fac, vec(3), g2, tmp

   amat = 0.0_wp
   fac = 4*pi/vol

   do itr = 1, size(trans, 2)
      vec(:) = trans(:, itr)
      g2 = dot_product(vec, vec)
      if (g2 < eps) cycle
      tmp = cos(dot_product(rij, vec)) * fac * exp(-0.25_wp*g2/(alp*alp))/g2
      amat = amat + tmp
   end do

end subroutine get_amat_rec_3d

subroutine get_coulomb_derivs(self, mol, cache, qvec, dadr, dadL, atrace)
   class(mchrg_klopman_coulomb), intent(in) :: self
   type(structure_type), intent(in) :: mol
   type(mchrg_coulomb_cache), intent(in) :: cache
   real(wp), intent(in) :: qvec(:)
   real(wp), intent(inout) :: dadr(:, :, :)
   real(wp), intent(inout) :: dadL(:, :, :)
   real(wp), intent(inout) :: atrace(:, :)

   if (any(mol%periodic)) then
      call get_damat_3d(mol, cache%wsc, cache%alpha, self%gexp, self%hardness, qvec, &
         & dadr, dadL, atrace)
   else
      call get_damat_0d(mol, self%gexp, self%hardness, qvec, dadr, dadL, atrace)
   end if
end subroutine get_coulomb_derivs

subroutine get_damat_0d(mol, gexp, hardness, qvec, dadr, dadL, atrace)
   type(structure_type), intent(in) :: mol
   real(wp), intent(in) :: gexp
   real(wp), intent(in) :: hardness(:, :)
   real(wp), intent(in) :: qvec(:)
   real(wp), intent(inout) :: dadr(:, :, :)
   real(wp), intent(inout) :: dadL(:, :, :)
   real(wp), intent(inout) :: atrace(:, :)

   integer :: iat, jat, izp, jzp
   real(wp) :: vec(3), r1, gam, arg, dtmp, dG(3), dS(3, 3)

   !$omp parallel do default(none) schedule(runtime) &
   !$omp reduction(+:atrace, dadr, dadL) shared(mol, gexp, hardness, qvec) &
   !$omp private(iat, izp, jat, jzp, gam, r1, vec, dG, dS, dtmp, arg)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      do jat = 1, iat-1
         jzp = mol%id(jat)
         vec = mol%xyz(:, iat) - mol%xyz(:, jat)
         r1 = norm2(vec)
         gam = hardness(jzp, izp)
         dtmp = 1.0_wp / (r1**gexp + gam**(-gexp))
         dtmp = -r1**(gExp-2.0_wp) * dtmp * dtmp**(1.0_wp/gExp)
         dG = dtmp*vec
         dS = spread(dG, 1, 3) * spread(vec, 2, 3)
         atrace(:, iat) = +dG*qvec(jat) + atrace(:, iat)
         atrace(:, jat) = -dG*qvec(iat) + atrace(:, jat)
         dadr(:, iat, jat) = +dG*qvec(iat) + dadr(:, iat, jat)
         dadr(:, jat, iat) = -dG*qvec(jat) + dadr(:, jat, iat)
         dadL(:, :, jat) = +dS*qvec(iat) + dadL(:, :, jat)
         dadL(:, :, iat) = +dS*qvec(jat) + dadL(:, :, iat)
      end do
   end do

end subroutine get_damat_0d

subroutine get_damat_3d(mol, wsc, alpha, gexp, hardness, qvec, dadr, dadL, atrace)
   type(structure_type), intent(in) :: mol
   type(wignerseitz_cell_type), intent(in) :: wsc
   real(wp), intent(in) :: alpha
   real(wp), intent(in) :: gexp
   real(wp), intent(in) :: hardness(:, :)
   real(wp), intent(in) :: qvec(:)
   real(wp), intent(inout) :: dadr(:, :, :)
   real(wp), intent(inout) :: dadL(:, :, :)
   real(wp), intent(inout) :: atrace(:, :)

   integer :: iat, jat, izp, jzp, img
   real(wp) :: vol, gam, wsw, vec(3), dG(3), dS(3, 3)
   real(wp) :: dGd(3), dSd(3, 3), dGr(3), dSr(3, 3)
   real(wp), allocatable :: dtrans(:, :), rtrans(:, :)

   vol = abs(matdet_3x3(mol%lattice))
   call get_dir_trans(mol%lattice, dtrans)
   call get_rec_trans(mol%lattice, rtrans)

   !$omp parallel do default(none) schedule(runtime) &
   !$omp reduction(+:atrace, dadr, dadL) &
   !$omp shared(mol, gexp, hardness, wsc, alpha, vol, dtrans, rtrans, qvec) &
   !$omp private(iat, izp, jat, jzp, img, gam, wsw, vec, dG, dS, &
   !$omp& dGr, dSr, dGd, dSd)
   do iat = 1, mol%nat
      izp = mol%id(iat)
      do jat = 1, iat-1
         jzp = mol%id(jat)
         dG(:) = 0.0_wp
         dS(:, :) = 0.0_wp
         gam = hardness(jzp, izp)
         wsw = 1.0_wp / real(wsc%nimg(jat, iat), wp)
         do img = 1, wsc%nimg(jat, iat)
            vec = mol%xyz(:, iat) - mol%xyz(:, jat) - wsc%trans(:, wsc%tridx(img, jat, iat))
            call get_damat_dir_3d(vec, gexp, gam, alpha, dtrans, dGd, dSd)
            call get_damat_rec_3d(vec, vol, alpha, rtrans, dGr, dSr)
            dG = dG + (dGd + dGr) * wsw
            dS = dS + (dSd + dSr) * wsw
         end do
         atrace(:, iat) = +dG*qvec(jat) + atrace(:, iat)
         atrace(:, jat) = -dG*qvec(iat) + atrace(:, jat)
         dadr(:, iat, jat) = +dG*qvec(iat) + dadr(:, iat, jat)
         dadr(:, jat, iat) = -dG*qvec(jat) + dadr(:, jat, iat)
         dadL(:, :, jat) = +dS*qvec(iat) + dadL(:, :, jat)
         dadL(:, :, iat) = +dS*qvec(jat) + dadL(:, :, iat)
      end do

      dS(:, :) = 0.0_wp
      gam = hardness(izp, izp)
      wsw = 1.0_wp / real(wsc%nimg(iat, iat), wp)
      do img = 1, wsc%nimg(iat, iat)
         vec = wsc%trans(:, wsc%tridx(img, iat, iat))
         call get_damat_dir_3d(vec, gexp, gam, alpha, dtrans, dGd, dSd)
         call get_damat_rec_3d(vec, vol, alpha, rtrans, dGr, dSr)
         dS = dS + (dSd + dSr) * wsw
      end do
      dadL(:, :, iat) = +dS*qvec(iat) + dadL(:, :, iat)
   end do

end subroutine get_damat_3d

subroutine get_damat_dir_3d(rij, gam, gexp, alp, trans, dg, ds)
   real(wp), intent(in) :: rij(3)
   real(wp), intent(in) :: gam
   real(wp), intent(in) :: gexp
   real(wp), intent(in) :: alp
   real(wp), intent(in) :: trans(:, :)
   real(wp), intent(out) :: dg(3)
   real(wp), intent(out) :: ds(3, 3)

   integer :: itr
   real(wp) :: vec(3), r1, r2, gtmp, atmp, gam2, alp2

   dg(:) = 0.0_wp
   ds(:, :) = 0.0_wp

   alp2 = alp*alp

   do itr = 1, size(trans, 2)
      vec(:) = rij + trans(:, itr)
      r1 = norm2(vec)
      if (r1 < eps) cycle
      r2 = r1*r1
      gtmp = 1.0_wp / (r1**gexp + gam**(-gexp))
      gtmp = -r1**(gexp-2.0_wp) * gtmp * gtmp**(1.0_wp/gexp)
      atmp = -2*alp*exp(-r2*alp2)/(sqrtpi*r2) + erf(r1*alp)/(r2*r1)
      dg(:) = dg + (gtmp + atmp) * vec
      ds(:, :) = ds + (gtmp + atmp) * spread(vec, 1, 3) * spread(vec, 2, 3)
   end do

end subroutine get_damat_dir_3d

subroutine get_damat_rec_3d(rij, vol, alp, trans, dg, ds)
   real(wp), intent(in) :: rij(3)
   real(wp), intent(in) :: vol
   real(wp), intent(in) :: alp
   real(wp), intent(in) :: trans(:, :)
   real(wp), intent(out) :: dg(3)
   real(wp), intent(out) :: ds(3, 3)

   integer :: itr
   real(wp) :: fac, vec(3), g2, gv, etmp, dtmp, alp2
   real(wp), parameter :: unity(3, 3) = reshape(&
      & [1, 0, 0, 0, 1, 0, 0, 0, 1], shape(unity))

   dg(:) = 0.0_wp
   ds(:, :) = 0.0_wp
   fac = 4*pi/vol
   alp2 = alp*alp

   do itr = 1, size(trans, 2)
      vec(:) = trans(:, itr)
      g2 = dot_product(vec, vec)
      if (g2 < eps) cycle
      gv = dot_product(rij, vec)
      etmp = fac * exp(-0.25_wp*g2/alp2)/g2
      dtmp = -sin(gv) * etmp
      dg(:) = dg + dtmp * vec
      ds(:, :) = ds + etmp * cos(gv) &
         & * ((2.0_wp/g2 + 0.5_wp/alp2) * spread(vec, 1, 3)*spread(vec, 2, 3) - unity)
   end do

end subroutine get_damat_rec_3d


end module multicharge_coulomb_klopman
