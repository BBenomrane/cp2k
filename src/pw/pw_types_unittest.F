!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

PROGRAM pw_types_unittest
   USE kinds,                           ONLY: dp
   USE message_passing,                 ONLY: mp_comm_type,&
                                              mp_para_env_release,&
                                              mp_para_env_type,&
                                              mp_world_finalize,&
                                              mp_world_init
   USE pw_grid_types,                   ONLY: FULLSPACE,&
                                              pw_grid_type
   USE pw_grids,                        ONLY: pw_grid_create,&
                                              pw_grid_release,&
                                              pw_grid_setup

   IMPLICIT NONE

   INTEGER                          :: i_arg, nargs, stat, npts
   LOGICAL                          :: ionode
   REAL(KIND=dp)                    :: cell_parameter, hmat(3, 3)
   TYPE(pw_grid_type), POINTER      :: grid
   TYPE(mp_comm_type) :: mpi_comm
   TYPE(mp_para_env_type), POINTER :: para_env
   CHARACTER(len=32)                :: arg

   npts = 200
   cell_parameter = 5.0
   nargs = command_argument_count()

   WRITE (*, *) nargs
   DO i_arg = 1, nargs
      CALL GET_COMMAND_ARGUMENT(i_arg, arg)
      WRITE (*, *) i_arg, ">", arg, "<"
   END DO

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

   CALL pw_grid_release(grid)

   CALL mp_para_env_release(para_env)
   CALL mp_world_finalize()

END PROGRAM pw_types_unitTEST
! vim: set ts=3 sw=3 tw=132 :
