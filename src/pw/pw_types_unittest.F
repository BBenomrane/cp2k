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
                       REALDATA3D, COMPLEXDATA3D
   USE pw_methods, ONLY: pw_copy, pw_copy_to_array, pw_copy_from_array, pw_set

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
   #:for kind in pw_kinds
      #{if kind[1]=="3"}# TYPE(pw_type) :: pw_${kind}$_real_old_1, pw_${kind}$_real_old_2#{endif}#
      TYPE(pw_${kind}$_type) :: pw_${kind}$_real_new_1, pw_${kind}$_real_new_2
      #{if kind[1]=="3"}# TYPE(pw_type) :: pw_${kind}$_rec_old_1, pw_${kind}$_rec_old_2#{endif}#
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
   CALL pw_c1d_real_new_1%create(grid, REALSPACE)
   CALL pw_c1d_real_new_2%create(grid, REALSPACE)
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
   CALL pw_c1d_rec_new_1%create(grid, RECIPROCALSPACE)
   CALL pw_c1d_rec_new_2%create(grid, RECIPROCALSPACE)
   CALL pw_r1d_rec_new_1%create(grid, RECIPROCALSPACE)
   CALL pw_r1d_rec_new_2%create(grid, RECIPROCALSPACE)

   CALL test_pw_copy()
   CALL test_pw_copy_array()

   CALL pw_c3d_real_old_1%release()
   CALL pw_c3d_real_old_2%release()
   CALL pw_c3d_real_new_1%release()
   CALL pw_c3d_real_new_2%release()
   CALL pw_r3d_real_old_1%release()
   CALL pw_r3d_real_old_2%release()
   CALL pw_r3d_real_new_1%release()
   CALL pw_r3d_real_new_2%release()
   CALL pw_c1d_real_new_1%release()
   CALL pw_c1d_real_new_2%release()
   CALL pw_r1d_real_new_1%release()
   CALL pw_r1d_real_new_2%release()
   CALL pw_c3d_rec_old_1%release()
   CALL pw_c3d_rec_old_2%release()
   CALL pw_c3d_rec_new_1%release()
   CALL pw_c3d_rec_new_2%release()
   CALL pw_r3d_rec_old_1%release()
   CALL pw_r3d_rec_old_2%release()
   CALL pw_r3d_rec_new_1%release()
   CALL pw_r3d_rec_new_2%release()
   CALL pw_c1d_rec_new_1%release()
   CALL pw_c1d_rec_new_2%release()
   CALL pw_r1d_rec_new_1%release()
   CALL pw_r1d_rec_new_2%release()

   CALL pw_grid_release(grid)

   CALL mp_para_env_release(para_env)
   CALL mp_world_finalize()

CONTAINS

   SUBROUTINE test_pw_copy()

      REAL(KIND=dp) :: error

      #:for space in ["real", "rec"]
         #:for kind, kind2 in pw_kinds2_sameD
            CALL pw_set(pw_${kind}$_${space}$_new_1, 1.0_dp)
            CALL pw_copy(pw_${kind}$_${space}$_new_1, pw_${kind}$_${space}$_new_2)

            error = MAXVAL(ABS(pw_${kind}$_${space}$_new_2%array - pw_${kind}$_${space}$_new_1%array))
            CALL mpi_comm%max(error)
            IF (error > 0.0_dp) &
               ERROR STOP "Copying new ${kind}$ to new ${kind}$ ${space}$ grid does not work!"

            #:if kind[1]=="3"
               CALL pw_copy(pw_${kind}$_${space}$_new_2, pw_${kind2}$_${space}$_old_1)

               error = MAXVAL(ABS(pw_${kind}$_${space}$_new_2%array - pw_${kind2}$_${space}$_old_1%${kind2comp(kind2)}$))
               CALL mpi_comm%max(error)
               IF (error > 0.0_dp) &
                  ERROR STOP "Copying new ${kind}$ to old ${kind2}$ ${space}$ grid does not work!"

               CALL pw_copy(pw_${kind2}$_${space}$_old_1, pw_${kind2}$_${space}$_old_2)

            error = MAXVAL(ABS(pw_${kind2}$_${space}$_old_1%${kind2comp(kind2)}$-pw_${kind2}$_${space}$_old_2%${kind2comp(kind2)}$))
               CALL mpi_comm%max(error)
               IF (error > 0.0_dp) &
                  ERROR STOP "Copying old ${kind2}$ to old ${kind2}$ ${space}$ grid does not work!"

               CALL pw_copy(pw_${kind2}$_${space}$_old_2, pw_${kind}$_${space}$_new_1)

               error = MAXVAL(ABS(pw_${kind2}$_${space}$_old_2%${kind2comp(kind2)}$-pw_${kind}$_${space}$_new_1%array))
               CALL mpi_comm%max(error)
               IF (error > 0.0_dp) &
                  ERROR STOP "Copying old ${kind2}$ to new ${kind}$ ${space}$ grid does not work!"
            #:endif
         #:endfor
      #:endfor

      WRITE (default_output_unit, "(T3,A)") "Testing pw_copy successful!"

   END SUBROUTINE test_pw_copy

   SUBROUTINE test_pw_copy_array()

      REAL(KIND=dp) :: error

      ! We use these routines to copy between realspace and reciprocal space grid
      ! This is usually not their purpose but here a mere test

      #:for kind, kind2 in pw_kinds2_sameD
         CALL pw_set(pw_${kind}$_${space}$_new_1, 2.0_dp)
         WRITE (*, *) "${kind}$", "${kind2}$", "new", SHAPE(pw_${kind}$_real_new_2%array), SHAPE(pw_${kind2}$_real_new_1%array)
         CALL pw_copy_to_array(pw_${kind}$_real_new_1, pw_${kind2}$_real_new_2%array)

         error = MAXVAL(ABS(pw_${kind}$_real_new_2%array - pw_${kind2}$_real_new_1%array))
         CALL mpi_comm%max(error)
         IF (error > 0.0_dp) &
            ERROR STOP "Copying new ${kind}$ to new array ${kind}$ grid does not work!"

         WRITE (*, *) "${kind}$", "${kind2}$", "new", SHAPE(pw_${kind}$_rec_new_1%array), SHAPE(pw_${kind2}$_real_new_2%array)
         CALL pw_copy_from_array(pw_${kind}$_rec_new_1, pw_${kind2}$_real_new_2%array)

         error = MAXVAL(ABS(pw_${kind}$_rec_new_1%array - pw_${kind2}$_real_new_2%array))
         CALL mpi_comm%max(error)
         IF (error > 0.0_dp) &
            ERROR STOP "Copying array ${kind}$ to new ${kind}$ grid does not work!"

         #:if kind[1]=="3"
            CALL pw_set(pw_${kind}$_real_old_1, 3.0_dp)
            WRITE (*, *) "${kind}$", "${kind2}$", "old", SHAPE(pw_${kind}$_real_old_2%${kind2comp(kind)}$), SHAPE(pw_${kind2}$_real_old_1%${kind2comp(kind2)}$)
            CALL pw_copy_to_array(pw_${kind}$_real_old_1, pw_${kind2}$_real_old_2%${kind2comp(kind2)}$)

            error = MAXVAL(ABS(pw_${kind}$_real_old_1%${kind2comp(kind)}$-pw_${kind2}$_real_old_2%${kind2comp(kind2)}$))
            CALL mpi_comm%max(error)
            IF (error > 0.0_dp) &
               ERROR STOP "Copying old ${kind}$ to array ${kind2}$ grid does not work!"

            WRITE (*, *) "${kind}$", "${kind2}$", "old", SHAPE(pw_${kind}$_rec_old_2%${kind2comp(kind)}$), SHAPE(pw_${kind2}$_real_old_2%${kind2comp(kind2)}$)
            CALL pw_copy_from_array(pw_${kind}$_rec_old_2, pw_${kind2}$_real_old_2%${kind2comp(kind2)}$)

            error = MAXVAL(ABS(pw_${kind}$_rec_old_2%${kind2comp(kind)}$-pw_${kind2}$_real_old_2%${kind2comp(kind2)}$))
            CALL mpi_comm%max(error)
            IF (error > 0.0_dp) &
               ERROR STOP "Copying array ${kind}$ to old ${kind}$ grid does not work!"
         #:endif
      #:endfor

      WRITE (default_output_unit, "(T3,A)") "Testing pw_copy_to_array and pw_copy_from_array successful!"

   END SUBROUTINE test_pw_copy_array

END PROGRAM pw_types_unitTEST
! vim: set ts=3 sw=3 tw=132 :
