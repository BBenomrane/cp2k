!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

PROGRAM pw_types_unittest
   USE kinds, ONLY: dp
   USE machine, ONLY: &
      default_output_unit
   USE message_passing, ONLY: mp_comm_type, &
                              mp_para_env_release, &
                              mp_para_env_type, &
                              mp_world_finalize, &
                              mp_world_init
   USE pw_grid_types, ONLY: FULLSPACE, &
                            pw_grid_type
   USE pw_grids, ONLY: pw_grid_create, &
                       pw_grid_release, &
                       pw_grid_setup
   USE pw_types, ONLY: REALSPACE, RECIPROCALSPACE, &
                       pw_type, pw_c1d_type, pw_r1d_type, pw_r3d_type, pw_c3d_type, &
                       REALDATA3D, COMPLEXDATA3D, REALDATA1D, COMPLEXDATA1D
   USE pw_methods, ONLY: pw_copy, pw_copy_to_array, pw_copy_from_array, pw_set, pw_gather, pw_scatter, &
                         pw_transfer, pw_integrate_function
   USE fft_tools, ONLY: fft3d, &
                        finalize_fft, &
                        init_fft

   #:include "pw_types.fypp"

   #:mute
      #:def kind2comp(kind)
         #{if kind[1]=="1"}#c${kind[0]}$#{else}#c${kind}$#{endif}#
      #:enddef
   #:endmute

   IMPLICIT NONE

   INTEGER                          :: i_arg, nargs, stat, npts
   LOGICAL                          :: ionode
   REAL(KIND=dp)                    :: cell_parameter, hmat(3, 3)
   TYPE(pw_grid_type), POINTER      :: grid
   TYPE(mp_comm_type) :: mpi_comm
   TYPE(mp_para_env_type), POINTER :: para_env
   CHARACTER(len=32)                :: arg
   CHARACTER(LEN=*), PARAMETER                      :: wisdom_file = "/etc/fftw/wisdom"
   #:for kind in pw_kinds
      TYPE(pw_type) :: pw_${kind}$_real_old_1, pw_${kind}$_real_old_2
      TYPE(pw_${kind}$_type) :: pw_${kind}$_real_new_1, pw_${kind}$_real_new_2
      TYPE(pw_type) :: pw_${kind}$_rec_old_1, pw_${kind}$_rec_old_2
      TYPE(pw_${kind}$_type) :: pw_${kind}$_rec_new_1, pw_${kind}$_rec_new_2
   #:endfor

   npts = 20
   cell_parameter = 5.0
   nargs = command_argument_count()

   IF (nargs > 1) &
      ERROR STOP "Usage pw_types_unittest"

   IF (.FALSE.) THEN
      i_arg = 0
      arg_loop: DO WHILE (i_arg < nargs)
         i_arg = i_arg + 1
         CALL GET_COMMAND_ARGUMENT(i_arg, arg, status=stat)
         IF (stat /= 0) &
            ERROR STOP "Usage pw_types_unittest [--npts <int>]"
         SELECT CASE (arg)
         CASE ("--npts")
            i_arg = i_arg + 1
            CALL GET_COMMAND_ARGUMENT(i_arg, arg, status=stat)
            READ (arg, *, iostat=stat) npts
            IF (stat /= 0) &
               ERROR STOP "Usage pw_types_unittest [--npts <int>]"
         END SELECT
      END DO arg_loop
   END IF

   CALL mp_world_init(mpi_comm)
   ionode = mpi_comm%is_source()

   CALL init_fft_library()

   CALL RANDOM_SEED(mpi_comm%mepos)

   hmat = RESHAPE((/cell_parameter, 0.0_dp, 0.0_dp, 0.0_dp, cell_parameter, 0.0_dp, 0.0_dp, 0.0_dp, cell_parameter/), (/3, 3/))

! pw_grid_types require reference countable communicators, so we have to convert it first
   ALLOCATE (para_env)
   para_env = mpi_comm

! Use the para_env to create the grid
   NULLIFY (grid)
   CALL pw_grid_create(grid, para_env)
   CALL pw_grid_setup(hmat, grid, grid_span=FULLSPACE, npts=(/npts, npts, npts/), fft_usage=.TRUE.) !, iounit=iw)

! Create some grids
   CALL pw_c3d_real_old_1%create(grid, COMPLEXDATA3D, REALSPACE)
   CALL pw_c3d_real_old_2%create(grid, COMPLEXDATA3D, REALSPACE)
   CALL pw_c3d_real_new_1%create(grid, REALSPACE)
   CALL pw_c3d_real_new_2%create(grid, REALSPACE)
   CALL pw_r3d_real_old_1%create(grid, REALDATA3D, REALSPACE)
   CALL pw_r3d_real_old_2%create(grid, REALDATA3D, REALSPACE)
   CALL pw_r3d_real_new_1%create(grid, REALSPACE)
   CALL pw_r3d_real_new_2%create(grid, REALSPACE)
   CALL pw_c1d_real_old_1%create(grid, COMPLEXDATA1D, REALSPACE)
   CALL pw_c1d_real_old_2%create(grid, COMPLEXDATA1D, REALSPACE)
   CALL pw_c1d_real_new_1%create(grid, REALSPACE)
   CALL pw_c1d_real_new_2%create(grid, REALSPACE)
   CALL pw_r1d_real_old_1%create(grid, REALDATA1D, REALSPACE)
   CALL pw_r1d_real_old_2%create(grid, REALDATA1D, REALSPACE)
   CALL pw_r1d_real_new_1%create(grid, REALSPACE)
   CALL pw_r1d_real_new_2%create(grid, REALSPACE)
   CALL pw_c3d_rec_old_1%create(grid, COMPLEXDATA3D, RECIPROCALSPACE)
   CALL pw_c3d_rec_old_2%create(grid, COMPLEXDATA3D, RECIPROCALSPACE)
   CALL pw_c3d_rec_new_1%create(grid, RECIPROCALSPACE)
   CALL pw_c3d_rec_new_2%create(grid, RECIPROCALSPACE)
   CALL pw_r3d_rec_old_1%create(grid, REALDATA3D, RECIPROCALSPACE)
   CALL pw_r3d_rec_old_2%create(grid, REALDATA3D, RECIPROCALSPACE)
   CALL pw_r3d_rec_new_1%create(grid, RECIPROCALSPACE)
   CALL pw_r3d_rec_new_2%create(grid, RECIPROCALSPACE)
   CALL pw_c1d_rec_old_1%create(grid, COMPLEXDATA1D, RECIPROCALSPACE)
   CALL pw_c1d_rec_old_2%create(grid, COMPLEXDATA1D, RECIPROCALSPACE)
   CALL pw_c1d_rec_new_1%create(grid, RECIPROCALSPACE)
   CALL pw_c1d_rec_new_2%create(grid, RECIPROCALSPACE)
   CALL pw_r1d_rec_old_1%create(grid, REALDATA1D, RECIPROCALSPACE)
   CALL pw_r1d_rec_old_2%create(grid, REALDATA1D, RECIPROCALSPACE)
   CALL pw_r1d_rec_new_1%create(grid, RECIPROCALSPACE)
   CALL pw_r1d_rec_new_2%create(grid, RECIPROCALSPACE)

   CALL test_pw_copy()
   CALL test_pw_copy_array()
   CALL test_pw_gather_scatter()
   CALL test_pw_transfer()

   #:for kind in pw_kinds
      #:for space in ["real", "rec"]
         #:for age in ["old", "new"]
            CALL pw_${kind}$_${space}$_${age}$_1%release()
            CALL pw_${kind}$_${space}$_${age}$_2%release()
         #:endfor
      #:endfor
   #:endfor

   CALL pw_grid_release(grid)

   CALL finalize_fft_library()

   CALL mp_para_env_release(para_env)
   CALL mp_world_finalize()

CONTAINS

   SUBROUTINE test_pw_copy()

      REAL(KIND=dp) :: error

      #:for space in ["real", "rec"]
         #:for kind, kind2 in pw_kinds2_sameD
            CALL setup_pw_${kind}$ (pw_${kind}$_${space}$_new_1)
            CALL pw_copy(pw_${kind}$_${space}$_new_1, pw_${kind}$_${space}$_new_2)

            error = MAXVAL(ABS(pw_${kind}$_${space}$_new_2%array - pw_${kind}$_${space}$_new_1%array))
            CALL mpi_comm%max(error)
            IF (error > 0.0_dp) &
               ERROR STOP "Copying new ${kind}$ to new ${kind}$ ${space}$ grid does not work!"

            CALL pw_copy(pw_${kind}$_${space}$_new_2, pw_${kind2}$_${space}$_old_1)

           error = MAXVAL(ABS(REAL(pw_${kind}$_${space}$_new_2%array - pw_${kind2}$_${space}$_old_1%${kind2comp(kind2)}$, KIND=dp)))
            CALL mpi_comm%max(error)
            IF (error > 0.0_dp) &
               ERROR STOP "Copying new ${kind}$ to old ${kind2}$ ${space}$ grid does not work!"

            CALL pw_copy(pw_${kind2}$_${space}$_old_1, pw_${kind2}$_${space}$_old_2)

            error = MAXVAL(ABS(pw_${kind2}$_${space}$_old_1%${kind2comp(kind2)}$-pw_${kind2}$_${space}$_old_2%${kind2comp(kind2)}$))
            CALL mpi_comm%max(error)
            IF (error > 0.0_dp) &
               ERROR STOP "Copying old ${kind2}$ to old ${kind2}$ ${space}$ grid does not work!"

            CALL pw_copy(pw_${kind2}$_${space}$_old_2, pw_${kind}$_${space}$_new_1)

            error = MAXVAL(ABS(REAL(pw_${kind2}$_${space}$_old_2%${kind2comp(kind2)}$-pw_${kind}$_${space}$_new_1%array, KIND=dp)))
            CALL mpi_comm%max(error)
            IF (error > 0.0_dp) &
               ERROR STOP "Copying old ${kind2}$ to new ${kind}$ ${space}$ grid does not work!"
         #:endfor
      #:endfor

      WRITE (default_output_unit, "(T3,A)") "Testing pw_copy successful!"

   END SUBROUTINE test_pw_copy

   SUBROUTINE test_pw_copy_array()

      REAL(KIND=dp) :: error

      ! We use these routines to copy between realspace and reciprocal space grid
      ! This is usually not their purpose but here a mere test

      #:for kind, kind2 in pw_kinds2_sameD
         CALL setup_pw_${kind}$ (pw_${kind}$_${space}$_new_1)
         CALL pw_copy_to_array(pw_${kind}$_real_new_1, pw_${kind2}$_real_new_2%array)

         error = MAXVAL(ABS(REAL(pw_${kind}$_real_new_1%array - pw_${kind2}$_real_new_2%array, KIND=dp)))
         CALL mpi_comm%max(error)
         IF (error > 0.0_dp) &
            ERROR STOP "Copying new ${kind}$ to new array ${kind2}$ grid does not work!"

         CALL pw_copy_from_array(pw_${kind}$_rec_new_1, pw_${kind2}$_real_new_2%array)

         error = MAXVAL(ABS(pw_${kind}$_rec_new_1%array - pw_${kind2}$_real_new_2%array))
         CALL mpi_comm%max(error)
         IF (error > 0.0_dp) &
            ERROR STOP "Copying array ${kind}$ to new ${kind}$ grid does not work!"

         CALL setup_pw_pw(pw_${kind}$_real_old_1)
         CALL pw_copy_to_array(pw_${kind}$_real_old_1, pw_${kind2}$_real_old_2%${kind2comp(kind2)}$)

         error = MAXVAL(ABS(REAL(pw_${kind}$_real_old_1%${kind2comp(kind)}$-pw_${kind2}$_real_old_2%${kind2comp(kind2)}$, KIND=dp)))
         CALL mpi_comm%max(error)
         IF (error > 0.0_dp) &
            ERROR STOP "Copying old ${kind}$ to array ${kind2}$ grid does not work!"

         CALL pw_copy_from_array(pw_${kind}$_rec_old_2, pw_${kind2}$_real_old_2%${kind2comp(kind2)}$)

         error = MAXVAL(ABS(REAL(pw_${kind}$_rec_old_2%${kind2comp(kind)}$-pw_${kind2}$_real_old_2%${kind2comp(kind2)}$, KIND=dp)))
         CALL mpi_comm%max(error)
         IF (error > 0.0_dp) &
            ERROR STOP "Copying array ${kind}$ to old ${kind}$ grid does not work!"
      #:endfor

      WRITE (default_output_unit, "(T3,A)") "Testing pw_copy_to_array and pw_copy_from_array successful!"

   END SUBROUTINE test_pw_copy_array

   SUBROUTINE test_pw_gather_scatter()
      REAL(KIND=dp) :: error

      #:for kind, kind2 in pw_kinds2
         #:if kind[1]=="1" and kind2[1]=="3"
            CALL setup_pw_${kind}$ (pw_${kind}$_rec_new_2)
            CALL pw_gather(pw_${kind}$_rec_new_2, pw_${kind2}$_rec_new_1%array)
            CALL pw_scatter(pw_${kind}$_rec_new_2, pw_${kind2}$_real_new_2%array)

            error = MAXVAL(ABS(REAL(pw_${kind2}$_rec_new_1%array - pw_${kind2}$_real_new_2%array, KIND=dp)))
            CALL mpi_comm%max(error)
            IF (error > 0.0_dp) &
               ERROR STOP "Scatter between ${kind2}$ and ${kind}$ does not work!"
         #:endif
      #:endfor

      WRITE (default_output_unit, "(T3,A)") "Testing pw_gather and pw_scatter successful!"

   END SUBROUTINE test_pw_gather_scatter

   SUBROUTINE test_pw_transfer()

      REAL(KIND=dp) :: error, tmp

! We try to transform and check the integral
      CALL setup_pw_r3d(pw_r3d_real_new_1)
      CALL pw_transfer(pw_r3d_real_new_1, pw_c1d_rec_new_1)
      tmp = pw_integrate_function(pw_r3d_real_new_1)
      error = ABS(tmp - pw_integrate_function(pw_c1d_rec_new_1))
      IF (error > 1.0E-12) THEN
         WRITE (default_output_unit, "(A,T64,ES18.10)") "FFT from new r3d real to new c1d reciprocal does not work: ", error
         ERROR STOP
      END IF

      CALL pw_copy(pw_r3d_real_new_1, pw_r3d_real_old_1)
! We try to transform and check the integral
      CALL pw_transfer(pw_r3d_real_old_1, pw_c1d_rec_old_1)
      tmp = pw_integrate_function(pw_r3d_real_old_1)
      error = ABS(tmp - pw_integrate_function(pw_c1d_rec_old_1))
      IF (error > 1.0E-12) THEN
         WRITE (default_output_unit, "(A,T64,ES18.10)") "FFT from new r3d real to old c1d reciprocal does not work: ", error
         ERROR STOP
      END IF

      error = MAXVAL(ABS(pw_c1d_rec_old_1%cc - pw_c1d_rec_new_1%array))
      IF (error > 1.0E-12) THEN
         WRITE (default_output_unit, "(A,T64,ES18.10)") "FFT from new r3d real to old c1d reciprocal does not work (2): ", error
         ERROR STOP
      END IF

      CALL pw_transfer(pw_c1d_rec_new_1, pw_r3d_real_new_2)
      tmp = pw_integrate_function(pw_c1d_rec_new_1)
      error = ABS(tmp - pw_integrate_function(pw_r3d_real_new_2))
      IF (error > 1.0E-12) THEN
         WRITE (default_output_unit, "(A,T28,F18.10)") "FFT from new c1d real to new r3d reciprocal does not work: ", error
         ERROR STOP
      END IF
      error = MAXVAL(ABS(pw_r3d_real_new_1%array - pw_r3d_real_new_2%array))
      IF (error > 1.0E-12) THEN
         WRITE (default_output_unit, "(A,T64,ES18.10)") "FFT from new c1d real to new r3d reciprocal does not work: ", error
         ERROR STOP
      END IF

      CALL pw_transfer(pw_c1d_rec_new_1, pw_r3d_real_old_2)
      tmp = pw_integrate_function(pw_c1d_rec_new_1)
      error = ABS(tmp - pw_integrate_function(pw_r3d_real_old_2))
      IF (error > 1.0E-12) THEN
         WRITE (default_output_unit, "(A,T28,F18.10)") "FFT from new c1d real to new r3d reciprocal does not work: ", error
         ERROR STOP
      END IF
      error = MAXVAL(ABS(pw_r3d_real_old_2%cr3d - pw_r3d_real_new_2%array))
      IF (error > 1.0E-12) THEN
         WRITE (default_output_unit, "(A,T64,ES18.10)") "FFT from new c1d real to new r3d reciprocal does not work: ", error
         ERROR STOP
      END IF

      CALL pw_transfer(pw_c1d_rec_old_1, pw_r3d_real_old_2)
      tmp = pw_integrate_function(pw_c1d_rec_old_1)
      error = ABS(tmp - pw_integrate_function(pw_r3d_real_old_2))
      IF (error > 1.0E-12) THEN
         WRITE (default_output_unit, "(A,T28,F18.10)") "FFT from old c1d real to old r3d reciprocal does not work: ", error
         ERROR STOP
      END IF
      error = MAXVAL(ABS(pw_r3d_real_new_1%array - pw_r3d_real_old_2%cr3d))
      IF (error > 1.0E-12) THEN
         WRITE (default_output_unit, "(A,T64,ES18.10)") "FFT from r3d real to c1d reciprocal does not work (2): ", error
         ERROR STOP
      END IF

      CALL pw_transfer(pw_c1d_rec_old_1, pw_r3d_real_new_2)
      tmp = pw_integrate_function(pw_c1d_rec_old_1)
      error = ABS(tmp - pw_integrate_function(pw_r3d_real_new_2))
      IF (error > 1.0E-12) THEN
         WRITE (default_output_unit, "(A,T28,F18.10)") "FFT from old c1d reciprocal to new r3d real does not work: ", error
         ERROR STOP
      END IF
      error = MAXVAL(ABS(pw_r3d_real_old_1%cr3d - pw_r3d_real_new_2%array))
      IF (error > 1.0E-12) THEN
         WRITE (default_output_unit, "(A,T64,ES18.10)") "FFT from r3d real to c1d reciprocal does not work (2): ", error
         ERROR STOP
      END IF

      WRITE (default_output_unit, "(T3,A)") "Testing pw_transfer successful!"

   END SUBROUTINE test_pw_transfer

   SUBROUTINE init_fft_library()
      ! This parameter is given by fftw_plan_estimate in input_constants which is not allowed to be referenced from here
      INTEGER, PARAMETER                                            :: fftw_plan_type = 1
      INTEGER                                            :: n(3), stat
      COMPLEX(KIND=dp), DIMENSION(4, 4, 4)               :: zz

      CALL init_fft("FFTW3", alltoall=.FALSE., fftsg_sizes=.TRUE., &
                    pool_limit=10, &
                    wisdom_file=wisdom_file, &
                    plan_style=fftw_plan_type)
      n = 4
      zz = 0.0_dp
      CALL fft3d(1, n, zz, status=stat)
      IF (stat /= 0) THEN
         CALL init_fft("FFTSG", alltoall=.FALSE., fftsg_sizes=.TRUE., &
                       pool_limit=10, &
                       wisdom_file=wisdom_file, &
                       plan_style=fftw_plan_type)
         CALL fft3d(1, n, zz, status=stat)
         IF (stat /= 0) &
            ERROR STOP "No FFT library available"
      END IF

   END SUBROUTINE init_fft_library

   SUBROUTINE finalize_fft_library()
      CALL finalize_fft(para_env, wisdom_file)
   END SUBROUTINE finalize_fft_library

   #:for kind in pw_kinds
      SUBROUTINE setup_pw_${kind}$ (pw)
         TYPE(pw_${kind}$_type), INTENT(INOUT) :: pw
         REAL(KIND=dp), DIMENSION(SIZE(pw%array)) :: array
         #:if kind[0]=="c"
            REAL(KIND=dp), DIMENSION(SIZE(pw%array)) :: array2
         #:endif

         CALL RANDOM_NUMBER(array)
         #:if kind[0]=="r"
            pw%array = RESHAPE(array, SHAPE(pw%array))
         #:else
            CALL RANDOM_NUMBER(array2)
            pw%array = CMPLX(RESHAPE(array, SHAPE(pw%array)), RESHAPE(array2, SHAPE(pw%array)), KIND=dp)
         #:endif

      END SUBROUTINE setup_pw_${kind}$
   #:endfor

   SUBROUTINE setup_pw_pw(pw)
      TYPE(pw_type), INTENT(INOUT) :: pw
      REAL(KIND=dp), DIMENSION(:), ALLOCATABLE :: array_1d, array_1d_2
      REAL(KIND=dp), DIMENSION(:, :, :), ALLOCATABLE :: array_3d, array_3d_2

      SELECT CASE (pw%in_use)
      CASE (REALDATA3D)
         CALL RANDOM_NUMBER(pw%cr3d)
      CASE (COMPLEXDATA3D)
         ALLOCATE (array_3d(SIZE(pw%cc3d, 1), SIZE(pw%cc3d, 2), SIZE(pw%cc3d, 3)), &
                   array_3d_2(SIZE(pw%cc3d, 1), SIZE(pw%cc3d, 2), SIZE(pw%cc3d, 3)))
         CALL RANDOM_NUMBER(array_3d)
         CALL RANDOM_NUMBER(array_3d_2)
         pw%cc3d = CMPLX(array_3d, array_3d_2, KIND=dp)
      CASE (REALDATA1D)
         CALL RANDOM_NUMBER(pw%cr)
      CASE (COMPLEXDATA1D)
         ALLOCATE (array_1d(SIZE(pw%cc, 1)), array_1d_2(SIZE(pw%cc, 1)))
         CALL RANDOM_NUMBER(array_1d)
         CALL RANDOM_NUMBER(array_1d_2)
         pw%cc = CMPLX(array_1d, array_1d_2, KIND=dp)
      END SELECT

   END SUBROUTINE setup_pw_pw

END PROGRAM pw_types_unitTEST
! vim: set ts=3 sw=3 tw=132 :
