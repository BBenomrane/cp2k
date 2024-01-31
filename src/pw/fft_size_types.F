!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \par History
!>      JGH (30-Nov-2000): ESSL FFT Library added
!>      JGH (05-Jan-2001): Added SGI library FFT
!>      JGH (14-Jan-2001): Added parallel 3d FFT
!>      JGH (10-Feb-2006): New interface type
!>      JGH (31-Mar-2008): Remove local allocates and reshapes (performance)
!>                         Possible problems can be related with setting arrays
!>                         not to zero
!>                         Some interfaces could be further simplified by avoiding
!>                         an initial copy. However, this assumes contiguous arrays
!>      IAB (15-Oct-2008): Moved mp_cart_sub calls out of cube_tranpose_* and into
!>                         fft_scratch type, reducing number of calls dramatically
!>      IAB (05-Dec-2008): Moved all other non-essential MPI calls into scratch type
!>      IAB (09-Jan-2009): Added fft_plan_type to store FFT data, including cached FFTW plans
!>      IAB (13-Feb-2009): Extended plan caching to serial 3D FFT (fft3d_s)
!>      IAB (09-Oct-2009): Added OpenMP directives to parallel 3D FFT
!>                         (c) The Numerical Algorithms Group (NAG) Ltd, 2008-2009 on behalf of the HECToR project
!> \author JGH
! **************************************************************************************************
MODULE fft_size_types
   USE message_passing,                 ONLY: mp_cart_type,&
                                              mp_comm_type

!$ USE OMP_LIB, ONLY: omp_get_max_threads, omp_get_thread_num, omp_get_num_threads

#include "../base/base_uses.f90"

   IMPLICIT NONE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'fft_size_types'

   ! Types for the pool of scratch data needed in FFT routines
   ! keep the subroutine "is_equal" up-to-date
   ! needs a default initialization
   TYPE fft_size_type
      INTEGER                              :: nx = 0, ny = 0, nz = 0
      INTEGER                              :: lmax = 0, mmax = 0, nmax = 0
      INTEGER                              :: mx1 = 0, mx2 = 0, mx3 = 0
      INTEGER                              :: my1 = 0, my2 = 0, my3 = 0
      INTEGER                              :: mz1 = 0, mz2 = 0, mz3 = 0
      INTEGER                              :: mcz1 = 0, mcz2 = 0, mcy3 = 0, mcx2 = 0
      INTEGER                              :: lg = 0, mg = 0
      INTEGER                              :: nbx = 0, nbz = 0
      INTEGER                              :: nmray = 0, nyzray = 0
      TYPE(mp_comm_type)                   :: gs_group = mp_comm_type()
      TYPE(mp_cart_type)                   :: rs_group = mp_cart_type()
      INTEGER, DIMENSION(2)                :: g_pos = 0, r_pos = 0, r_dim = 0
      INTEGER                              :: numtask = 0
   CONTAINS
      PROCEDURE, PRIVATE, PASS(this), NON_OVERRIDABLE :: fft_size_is_equal
      GENERIC, PUBLIC :: OPERATOR(==) => fft_size_is_equal
   END TYPE fft_size_type

   PRIVATE
   PUBLIC :: fft_size_type

CONTAINS

! **************************************************************************************************
!> \brief  test data structures for equality. It is assumed that if they are
!>         different for one mpi task they are different for all (??)
!> \param this ...
!> \param other ...
!> \return ...
! **************************************************************************************************
   ELEMENTAL FUNCTION fft_size_is_equal(this, other) RESULT(equal)
      CLASS(fft_size_type), INTENT(IN)                    :: this, other
      LOGICAL                                            :: equal

      equal = .TRUE.

      equal = equal .AND. this%nx == other%nx
      equal = equal .AND. this%ny == other%ny
      equal = equal .AND. this%nz == other%nz

      equal = equal .AND. this%lmax == other%lmax
      equal = equal .AND. this%mmax == other%mmax
      equal = equal .AND. this%nmax == other%nmax

      equal = equal .AND. this%mx1 == other%mx1
      equal = equal .AND. this%mx2 == other%mx2
      equal = equal .AND. this%mx3 == other%mx3

      equal = equal .AND. this%my1 == other%my1
      equal = equal .AND. this%my2 == other%my2
      equal = equal .AND. this%my3 == other%my3

      equal = equal .AND. this%mcz1 == other%mcz1
      equal = equal .AND. this%mcx2 == other%mcx2
      equal = equal .AND. this%mcz2 == other%mcz2
      equal = equal .AND. this%mcy3 == other%mcy3

      equal = equal .AND. this%lg == other%lg
      equal = equal .AND. this%mg == other%mg

      equal = equal .AND. this%nbx == other%nbx
      equal = equal .AND. this%nbz == other%nbz

      equal = equal .AND. this%nmray == other%nmray
      equal = equal .AND. this%nyzray == other%nyzray

      equal = equal .AND. this%gs_group == other%gs_group
      equal = equal .AND. this%rs_group == other%rs_group

      equal = equal .AND. ALL(this%g_pos == other%g_pos)
      equal = equal .AND. ALL(this%r_pos == other%r_pos)
      equal = equal .AND. ALL(this%r_dim == other%r_dim)

      equal = equal .AND. this%numtask == other%numtask

   END FUNCTION fft_size_is_equal

END MODULE fft_size_types
