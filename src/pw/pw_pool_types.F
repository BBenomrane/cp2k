!--------------------------------------------------------------------------------------------------!
!   CP2K: A general program to perform molecular dynamics simulations                              !
!   Copyright 2000-2024 CP2K developers group <https://cp2k.org>                                   !
!                                                                                                  !
!   SPDX-License-Identifier: GPL-2.0-or-later                                                      !
!--------------------------------------------------------------------------------------------------!

! **************************************************************************************************
!> \brief Manages a pool of grids (to be used for example as tmp objects),
!>      but can also be used to instantiate grids that are never given back.
!>
!>      Multigrid pools are just an array of pw_pools
!> \note
!>      The pool could also work without pointers (doing = each time),
!>      but I find it *very* ugly.
!>
!>      The pool could be integrated into pw_grid_type, I don't know if
!>      it would be a good or bad idea (but would add a circular dependence
!>      between pw and pw_grid types).
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
MODULE pw_pool_types
   USE cp_linked_list_pw, ONLY: &
      cp_sll_3d_c_dealloc, cp_sll_3d_c_get_first_el, cp_sll_3d_c_get_length, &
      cp_sll_3d_c_insert_el, cp_sll_3d_c_next, cp_sll_3d_c_rm_first_el, cp_sll_3d_c_type, &
      cp_sll_3d_r_dealloc, cp_sll_3d_r_get_first_el, cp_sll_3d_r_get_length, &
      cp_sll_3d_r_insert_el, cp_sll_3d_r_next, cp_sll_3d_r_rm_first_el, cp_sll_3d_r_type, &
      cp_sll_pw_dealloc, cp_sll_pw_get_first_el, cp_sll_pw_get_length, cp_sll_pw_insert_el, &
      cp_sll_pw_next, cp_sll_pw_rm_first_el, cp_sll_pw_type, &
      cp_sll_1d_c_dealloc, cp_sll_1d_c_get_first_el, cp_sll_1d_c_get_length, &
      cp_sll_1d_c_insert_el, cp_sll_1d_c_next, cp_sll_1d_c_rm_first_el, cp_sll_1d_c_type, &
      cp_sll_1d_r_dealloc, cp_sll_1d_r_get_first_el, cp_sll_1d_r_get_length, &
      cp_sll_1d_r_insert_el, cp_sll_1d_r_next, cp_sll_1d_r_rm_first_el, cp_sll_1d_r_type
   USE kinds, ONLY: dp
   USE pw_grid_types, ONLY: pw_grid_type
   USE pw_grids, ONLY: pw_grid_compare, &
                       pw_grid_release, &
                       pw_grid_retain
   USE pw_types, ONLY: COMPLEXDATA1D, &
                       COMPLEXDATA3D, &
                       REALDATA1D, &
                       REALDATA3D, &
                       REALSPACE, &
                       RECIPROCALSPACE, &
                       pw_type
#include "../base/base_uses.f90"

   #:include 'pw_types.fypp'

   IMPLICIT NONE
   PRIVATE

   CHARACTER(len=*), PARAMETER, PRIVATE :: moduleN = 'pw_pool_types'
   INTEGER, PARAMETER :: default_max_cache = 75, max_max_cache = 150

   PUBLIC :: pw_pool_type, pw_pool_p_type
   PUBLIC :: pw_pool_create, pw_pool_release
   PUBLIC :: pw_pools_copy, pw_pools_dealloc, &
             pw_pools_create_pws, pw_pools_give_back_pws

! **************************************************************************************************
!> \brief Manages a pool of grids (to be used for example as tmp objects),
!>      but can also be used to instantiate grids that are never given back.
!> \param ref_count reference count (see /cp2k/doc/ReferenceCounting.html)
!> \param real 1d_array, c1d_array, complex3d_array: liked list with
!>        the cached grids of the corresponding type
!> \note
!>      As of now I would like replace the linked lists by arrays
!>      (no annoying list elements that are allocated would show up when
!>      tracking leaks) [fawzi]
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
   TYPE pw_pool_type
      INTEGER :: ref_count = 0, max_cache = 0
      TYPE(pw_grid_type), POINTER :: pw_grid => NULL()
      TYPE(cp_sll_3d_r_type), POINTER :: r3d_array => NULL()
      TYPE(cp_sll_3d_c_type), POINTER :: c3d_array => NULL()
      TYPE(cp_sll_1d_r_type), POINTER :: r1d_array => NULL()
      TYPE(cp_sll_1d_c_type), POINTER :: c1d_array => NULL()
   CONTAINS
      PROCEDURE, PUBLIC, NON_OVERRIDABLE :: retain => pw_pool_retain
      PROCEDURE, PUBLIC, NON_OVERRIDABLE :: create_pw => pw_pool_create_pw
      PROCEDURE, PUBLIC, NON_OVERRIDABLE :: give_back_pw => pw_pool_give_back_pw
      PROCEDURE, PUBLIC, NON_OVERRIDABLE :: create_cr3d => pw_pool_create_cr3d
      PROCEDURE, PUBLIC, NON_OVERRIDABLE :: give_back_cr3d => pw_pool_give_back_cr3d
   END TYPE pw_pool_type

! **************************************************************************************************
!> \brief to create arrays of pools
!> \param pool the pool
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
   TYPE pw_pool_p_type
      TYPE(pw_pool_type), POINTER :: pool => NULL()
   END TYPE pw_pool_p_type

CONTAINS

! **************************************************************************************************
!> \brief creates a pool for pw
!> \param pool the pool to create
!> \param pw_grid the grid that is used to create the pw
!> \param max_cache ...
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
   SUBROUTINE pw_pool_create(pool, pw_grid, max_cache)
      TYPE(pw_pool_type), POINTER                        :: pool
      TYPE(pw_grid_type), POINTER                        :: pw_grid
      INTEGER, OPTIONAL                                  :: max_cache

      ALLOCATE (pool)
      pool%pw_grid => pw_grid
      CALL pw_grid_retain(pw_grid)
      pool%ref_count = 1
      pool%max_cache = default_max_cache
      IF (PRESENT(max_cache)) pool%max_cache = max_cache
      pool%max_cache = MIN(max_max_cache, pool%max_cache)
   END SUBROUTINE pw_pool_create

! **************************************************************************************************
!> \brief retains the pool (see cp2k/doc/ReferenceCounting.html)
!> \param pool the pool to retain
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
   SUBROUTINE pw_pool_retain(pool)
      CLASS(pw_pool_type), INTENT(INOUT)                  :: pool

      CPASSERT(pool%ref_count > 0)

      pool%ref_count = pool%ref_count + 1
   END SUBROUTINE pw_pool_retain

! **************************************************************************************************
!> \brief deallocates all the cached grids
!> \param pool the pool to flush
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
   SUBROUTINE pw_pool_flush_cache(pool)
      TYPE(pw_pool_type), INTENT(INOUT)                  :: pool

      #:for kind, type in zip(pw_kinds, pw_types)
         ${type}$, CONTIGUOUS, POINTER                      :: ${kind}$_att
         TYPE(cp_sll_${kind[1:]}$_${kind[0]}$_type), POINTER   :: ${kind}$_iterator
      #:endfor

      #:for kind in pw_kinds
         NULLIFY (${kind}$_iterator, ${kind}$_att)
         ${kind}$_iterator => pool%${kind}$_array
         DO
            IF (.NOT. cp_sll_${kind[1:]}$_${kind[0]}$_next(${kind}$_iterator, el_att=${kind}$_att)) EXIT
            DEALLOCATE (${kind}$_att)
         END DO
         CALL cp_sll_${kind[1:]}$_${kind[0]}$_dealloc(pool%${kind}$_array)
      #:endfor

   END SUBROUTINE pw_pool_flush_cache

! **************************************************************************************************
!> \brief releases the given pool (see cp2k/doc/ReferenceCounting.html)
!> \param pool the pool to release
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
   SUBROUTINE pw_pool_release(pool)
      TYPE(pw_pool_type), POINTER                        :: pool

      IF (ASSOCIATED(pool)) THEN
         CPASSERT(pool%ref_count > 0)
         pool%ref_count = pool%ref_count - 1
         IF (pool%ref_count == 0) THEN
            CALL pw_pool_flush_cache(pool)
            CALL pw_grid_release(pool%pw_grid)

            DEALLOCATE (pool)
         END IF
      END IF
      NULLIFY (pool)
   END SUBROUTINE pw_pool_release

   #:for kind, type in zip(pw_kinds, pw_types)
! **************************************************************************************************
!> \brief tries to pop an element from the given list (no error on failure)
!> \param list the list to pop
!> \return ...
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
!> \note
!>      private function
! **************************************************************************************************
      FUNCTION try_pop_${kind}$ (list) RESULT(res)
         TYPE(cp_sll_${kind[1:]}$_${kind[0]}$_type), POINTER                    :: list
         ${type}$, CONTIGUOUS, POINTER                                         :: res

         IF (ASSOCIATED(list)) THEN
            res => cp_sll_${kind[1:]}$_${kind[0]}$_get_first_el(list)
            CALL cp_sll_${kind[1:]}$_${kind[0]}$_rm_first_el(list)
         ELSE
            NULLIFY (res)
         END IF
      END FUNCTION try_pop_${kind}$
   #:endfor

! **************************************************************************************************
!> \brief returns a pw, allocating it if none is in the pool
!> \param pool the pool from where you get the pw
!> \param pw will contain the new pw
!> \param use_data which data it uses: REALDATA1D, COMPLEXDATA1D,
!>        REALDATA3D, COMPLEXDATA3D
!> \param in_space in which space it is: REALSPACE, RECIPROCALSPACE
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
   SUBROUTINE pw_pool_create_pw(pool, pw, use_data, in_space)
      CLASS(pw_pool_type), INTENT(IN)                     :: pool
      TYPE(pw_type), INTENT(OUT)                         :: pw
      INTEGER, INTENT(in)                                :: use_data, in_space

      CHARACTER(len=*), PARAMETER                        :: routineN = 'pw_pool_create_pw'

      INTEGER                                            :: handle
      #:for kind, type in zip(pw_kinds, pw_types)
         ${type}$, CONTIGUOUS, POINTER                      :: ${kind}$_ptr
      #:endfor

      CALL timeset(routineN, handle)
      #:for kind in pw_kinds
         NULLIFY (${kind}$_ptr)
      #:endfor

      SELECT CASE (use_data)
      CASE (REALDATA1D)
         r1d_ptr => try_pop_r1d(pool%r1d_array)
         CALL pw%create(pool%pw_grid, use_data=use_data, &
                        in_space=in_space, r1d_ptr=r1d_ptr)
      CASE (REALDATA3D)
         r3d_ptr => try_pop_r3d(pool%r3d_array)
         CALL pw%create(pool%pw_grid, use_data=use_data, &
                        in_space=in_space, r3d_ptr=r3d_ptr)
      CASE (COMPLEXDATA1D)
         c1d_ptr => try_pop_c1d(pool%c1d_array)
         CALL pw%create(pool%pw_grid, use_data=use_data, &
                        in_space=in_space, c1d_ptr=c1d_ptr)
      CASE (COMPLEXDATA3D)
         c3d_ptr => try_pop_c3d(pool%c3d_array)
         CALL pw%create(pool%pw_grid, use_data=use_data, &
                        in_space=in_space, c3d_ptr=c3d_ptr)
      CASE default
! unknown use_data
         CPABORT("")
      END SELECT

      CALL timestop(handle)

   END SUBROUTINE pw_pool_create_pw

! **************************************************************************************************
!> \brief returns the pw to the pool
!> \param pool the pool where to reintegrate the pw
!> \param pw the pw to give back
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
   SUBROUTINE pw_pool_give_back_pw(pool, pw)
      CLASS(pw_pool_type), INTENT(IN)                     :: pool
      TYPE(pw_type), INTENT(INOUT)                       :: pw

      CHARACTER(len=*), PARAMETER :: routineN = 'pw_pool_give_back_pw'

      INTEGER                                            :: handle

      CALL timeset(routineN, handle)
      IF (ASSOCIATED(pw%pw_grid)) THEN
         IF (pw_grid_compare(pw%pw_grid, pool%pw_grid)) THEN
            SELECT CASE (pw%in_use)
            CASE (REALDATA1D)
               IF (ASSOCIATED(pw%cr)) THEN
                  IF (cp_sll_1d_r_get_length(pool%r1d_array) < pool%max_cache) THEN
                     CALL cp_sll_1d_r_insert_el(pool%r1d_array, el=pw%cr)
                     NULLIFY (pw%cr)
                  ELSE IF (max_max_cache >= 0) THEN
                     CPWARN("hit max_cache")
                  END IF
               END IF
            CASE (REALDATA3D)
               IF (ASSOCIATED(pw%cr3d)) THEN
                  IF (cp_sll_3d_r_get_length(pool%r3d_array) < pool%max_cache) THEN
                     CALL cp_sll_3d_r_insert_el(pool%r3d_array, el=pw%cr3d)
                     NULLIFY (pw%cr3d)
                  ELSE IF (max_max_cache >= 0) THEN
                     CPWARN("hit max_cache")
                  END IF
               END IF
            CASE (COMPLEXDATA1D)
               IF (ASSOCIATED(pw%cc)) THEN
                  IF (cp_sll_1d_c_get_length(pool%c1d_array) < pool%max_cache) THEN
                     CALL cp_sll_1d_c_insert_el(pool%c1d_array, el=pw%cc)
                     NULLIFY (pw%cc)
                  ELSE IF (max_max_cache >= 0) THEN
                     CPWARN("hit max_cache")
                  END IF
               END IF
            CASE (COMPLEXDATA3D)
               IF (ASSOCIATED(pw%cc3d)) THEN
                  IF (cp_sll_3d_c_get_length(pool%c3d_array) < pool%max_cache) THEN
                     CALL cp_sll_3d_c_insert_el(pool%c3d_array, el=pw%cc3d)
                     NULLIFY (pw%cc3d)
                  ELSE IF (max_max_cache >= 0) THEN
                     CPWARN("hit max_cache")
                  END IF
               END IF
            CASE default
               ! unknown in_use
               CPABORT("")
            END SELECT
         END IF
      END IF
      CALL pw%release()
      CALL timestop(handle)
   END SUBROUTINE pw_pool_give_back_pw

! **************************************************************************************************
!> \brief returns a 3d real array of coefficients as the one used by pw with
!>      REALDATA3D, allocating it if none is present in the pool
!> \param pw_pool the pool that caches the cr3d
!> \param cr3d the pointer that will contain the array
!> \par History
!>      11.2003 created [fawzi]
!> \author fawzi
! **************************************************************************************************
   SUBROUTINE pw_pool_create_cr3d(pw_pool, cr3d)
      CLASS(pw_pool_type), INTENT(IN)                     :: pw_pool
      REAL(kind=dp), DIMENSION(:, :, :), POINTER         :: cr3d

      TYPE(pw_type), POINTER                             :: pw

      NULLIFY (pw)
      IF (ASSOCIATED(pw_pool%r3d_array)) THEN
         cr3d => cp_sll_3d_r_get_first_el(pw_pool%r3d_array)
         CALL cp_sll_3d_r_rm_first_el(pw_pool%r3d_array)
      END IF
      IF (.NOT. ASSOCIATED(cr3d)) THEN
         ALLOCATE (cr3d(pw_pool%pw_grid%bounds_local(1, 1):pw_pool%pw_grid%bounds_local(2, 1), &
                        pw_pool%pw_grid%bounds_local(1, 2):pw_pool%pw_grid%bounds_local(2, 2), &
                        pw_pool%pw_grid%bounds_local(1, 3):pw_pool%pw_grid%bounds_local(2, 3)))
      END IF
   END SUBROUTINE pw_pool_create_cr3d

! **************************************************************************************************
!> \brief returns a 3d real array of coefficients as the one used by pw with
!>      REALDATA3D, allocating it if none is present in the pool
!> \param pw_pool the pool that caches the cr3d
!> \param cr3d the pointer that will contain the array
!> \par History
!>      11.2003 created [fawzi]
!> \author fawzi
! **************************************************************************************************
   SUBROUTINE pw_pool_give_back_cr3d(pw_pool, cr3d)
      CLASS(pw_pool_type), INTENT(IN)                     :: pw_pool
      REAL(kind=dp), CONTIGUOUS, DIMENSION(:, :, :), &
         POINTER                                         :: cr3d

      LOGICAL                                            :: compatible

      IF (ASSOCIATED(cr3d)) THEN
         compatible = ALL(MERGE(pw_pool%pw_grid%bounds_local(1, :) == LBOUND(cr3d) .AND. &
                                pw_pool%pw_grid%bounds_local(2, :) == UBOUND(cr3d), &
                                pw_pool%pw_grid%bounds_local(2, :) < pw_pool%pw_grid%bounds_local(1, :), &
                                UBOUND(cr3d) >= LBOUND(cr3d)))
         IF (compatible) THEN
            IF (cp_sll_3d_r_get_length(pw_pool%r3d_array) < pw_pool%max_cache) THEN
               CALL cp_sll_3d_r_insert_el(pw_pool%r3d_array, el=cr3d)
            ELSE
               IF (max_max_cache >= 0) &
                  CPWARN("hit max_cache")
               DEALLOCATE (cr3d)
            END IF
         ELSE
            DEALLOCATE (cr3d)
         END IF
      END IF
      NULLIFY (cr3d)
   END SUBROUTINE pw_pool_give_back_cr3d

! **************************************************************************************************
!> \brief creates a multigrid structure
!> \param pools the multigrid pool (i.e. an array of pw_pool)
!> \param pws the multigrid of coefficent you want to initialize
!> \param use_data which data it uses: REALDATA1D, COMPLEXDATA1D,
!>        REALDATA3D, COMPLEXDATA3D
!> \param in_space ...
!> \par History
!>      07.2004 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
   SUBROUTINE pw_pools_create_pws(pools, pws, use_data, in_space)
      TYPE(pw_pool_p_type), DIMENSION(:), INTENT(IN)     :: pools
      TYPE(pw_type), ALLOCATABLE, DIMENSION(:), &
         INTENT(OUT)                                     :: pws
      INTEGER, INTENT(in)                                :: use_data
      INTEGER, INTENT(in), OPTIONAL                      :: in_space

      INTEGER                                            :: i

      ALLOCATE (pws(SIZE(pools)))
      DO i = 1, SIZE(pools)
         CALL pools(i)%pool%create_pw(pws(i), use_data, &
                                      in_space=in_space)
      END DO
   END SUBROUTINE pw_pools_create_pws

! **************************************************************************************************
!> \brief returns the pw part of the coefficients into the pools
!> \param pools the pools that will cache the pws %pw
!> \param pws the coefficients to give back
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
   SUBROUTINE pw_pools_give_back_pws(pools, pws)
      TYPE(pw_pool_p_type), DIMENSION(:), INTENT(IN)     :: pools
      TYPE(pw_type), ALLOCATABLE, DIMENSION(:), &
         INTENT(INOUT)                                   :: pws

      INTEGER                                            :: i

      CPASSERT(SIZE(pws) == SIZE(pools))
      DO i = 1, SIZE(pools)
         CALL pools(i)%pool%give_back_pw(pws(i))
      END DO
      DEALLOCATE (pws)
   END SUBROUTINE pw_pools_give_back_pws

! **************************************************************************************************
!> \brief copies a multigrid pool, the underlying pools are shared
!> \param source_pools the pools to copy
!> \param target_pools will hold the copy of the pools
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
   SUBROUTINE pw_pools_copy(source_pools, target_pools)
      TYPE(pw_pool_p_type), DIMENSION(:), INTENT(IN)     :: source_pools
      TYPE(pw_pool_p_type), DIMENSION(:), POINTER        :: target_pools

      INTEGER                                            :: i

      ALLOCATE (target_pools(SIZE(source_pools)))
      DO i = 1, SIZE(source_pools)
         target_pools(i)%pool => source_pools(i)%pool
         CALL source_pools(i)%pool%retain()
      END DO
   END SUBROUTINE pw_pools_copy

! **************************************************************************************************
!> \brief deallocates the given pools (releasing each of the underlying
!>      pools)
!> \param pools the pols to deallocate
!> \par History
!>      08.2002 created [fawzi]
!> \author Fawzi Mohamed
! **************************************************************************************************
   SUBROUTINE pw_pools_dealloc(pools)
      TYPE(pw_pool_p_type), DIMENSION(:), POINTER        :: pools

      INTEGER                                            :: i

      IF (ASSOCIATED(pools)) THEN
         DO i = 1, SIZE(pools)
            CALL pw_pool_release(pools(i)%pool)
         END DO
         DEALLOCATE (pools)
      END IF
      NULLIFY (pools)
   END SUBROUTINE pw_pools_dealloc

END MODULE pw_pool_types
