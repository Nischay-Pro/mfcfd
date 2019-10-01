!        Generated by TAPENADE     (INRIA, Ecuador team)
!  Tapenade 3.14 (r7079) -  5 Oct 2018 09:56
!
!        Generated by TAPENADE     (INRIA, Ecuador team)
!  Tapenade 3.14 (r7079) -  5 Oct 2018 09:56
!
MODULE Q_LSKUM_MOD_DIFF_DIFF
  USE DATA_STRUCTURE_MOD_DIFF_DIFF
  USE POINT_NORMALS_MOD_DIFF_DIFF
  USE GENERATE_CONNECTIVITY_MOD_DIFF_DIFF
  USE FPI_SOLVER_MOD_DIFF_DIFF
  USE INITIAL_CONDITIONS_MOD_DIFF_DIFF
  IMPLICIT NONE

CONTAINS
!  Differentiation of q_lskum_d in reverse (adjoint) mode (with options fixinterface):
!   gradient     of useful results: cmd cdd cld
!   with respect to varying inputs: point.x point.y
!   RW status of diff variables: pointd.prim:(loc) pointd.prim_old:(loc)
!                pointd.flux_res:(loc) pointd.q:(loc) pointd.dq:(loc)
!                pointd.qm:(loc) pointd.delta:(loc) cmd:in-killed
!                cdd:in-killed point.x:out point.y:out point.nx:(loc)
!                point.ny:(loc) point.prim:(loc) point.prim_old:(loc)
!                point.flux_res:(loc) point.q:(loc) point.dq:(loc)
!                point.qm:(loc) point.delta:(loc) cld:in-killed
!  Differentiation of q_lskum in forward (tangent) mode (with options fixinterface):
!   variations   of useful results: cd cl cm
!   with respect to varying inputs: aoa
!   RW status of diff variables: q_inf:(loc) theta:(loc) aoa:in
!                q_init:(loc) cd:out cl:out cm:out point.prim:(loc)
!                point.prim_old:(loc) point.flux_res:(loc) point.q:(loc)
!                point.dq:(loc) point.qm:(loc) point.delta:(loc)
  SUBROUTINE Q_LSKUM_D_B()
    IMPLICIT NONE
    INTEGER :: i
    CALL INITIAL_CONDITIONS_D()
    CALL COMPUTE_NORMALS()
    CALL GENERATE_CONNECTIVITY()
    IF (restart .EQ. 0) THEN
      itr = 0
      pointd%prim_old = 0.0_8
      pointd%flux_res = 0.0_8
      pointd%q = 0.0_8
      pointd%dq = 0.0_8
      pointd%qm = 0.0_8
      pointd%delta = 0.0_8
    ELSE
      pointd%prim_old = 0.0_8
      pointd%flux_res = 0.0_8
      pointd%q = 0.0_8
      pointd%dq = 0.0_8
      pointd%qm = 0.0_8
      pointd%delta = 0.0_8
    END IF
    DO it=itr+1,itr+max_iters
      CALL PUSHREAL8ARRAY(point%delta, max_points)
      CALL PUSHREAL8ARRAY(point%qm, 2*4*max_points)
      CALL PUSHREAL8ARRAY(point%dq, 2*4*max_points)
      CALL PUSHREAL8ARRAY(point%q, 4*max_points)
      CALL PUSHREAL8ARRAY(point%flux_res, 4*max_points)
      CALL PUSHREAL8ARRAY(point%prim_old, 4*max_points)
      CALL PUSHREAL8ARRAY(point%prim, 4*max_points)
      CALL PUSHREAL8ARRAY(pointd%delta, max_points)
      CALL PUSHREAL8ARRAY(pointd%qm, 2*4*max_points)
      CALL PUSHREAL8ARRAY(pointd%dq, 2*4*max_points)
      CALL PUSHREAL8ARRAY(pointd%q, 4*max_points)
      CALL PUSHREAL8ARRAY(pointd%flux_res, 4*max_points)
      CALL PUSHREAL8ARRAY(pointd%prim_old, 4*max_points)
      CALL PUSHREAL8ARRAY(pointd%prim, 4*max_points)
      CALL FPI_SOLVER_D(it)
    END DO
    pointdb%prim = 0.0_8
    pointdb%prim_old = 0.0_8
    pointdb%flux_res = 0.0_8
    pointdb%q = 0.0_8
    pointdb%dq = 0.0_8
    pointdb%qm = 0.0_8
    pointdb%delta = 0.0_8
    pointb%x = 0.0_8
    pointb%y = 0.0_8
    pointb%nx = 0.0_8
    pointb%ny = 0.0_8
    pointb%prim = 0.0_8
    pointb%prim_old = 0.0_8
    pointb%flux_res = 0.0_8
    pointb%q = 0.0_8
    pointb%dq = 0.0_8
    pointb%qm = 0.0_8
    pointb%delta = 0.0_8
    DO it=itr+max_iters,itr+1,-1
      CALL POPREAL8ARRAY(pointd%prim, 4*max_points)
      CALL POPREAL8ARRAY(pointd%prim_old, 4*max_points)
      CALL POPREAL8ARRAY(pointd%flux_res, 4*max_points)
      CALL POPREAL8ARRAY(pointd%q, 4*max_points)
      CALL POPREAL8ARRAY(pointd%dq, 2*4*max_points)
      CALL POPREAL8ARRAY(pointd%qm, 2*4*max_points)
      CALL POPREAL8ARRAY(pointd%delta, max_points)
      CALL POPREAL8ARRAY(point%prim, 4*max_points)
      CALL POPREAL8ARRAY(point%prim_old, 4*max_points)
      CALL POPREAL8ARRAY(point%flux_res, 4*max_points)
      CALL POPREAL8ARRAY(point%q, 4*max_points)
      CALL POPREAL8ARRAY(point%dq, 2*4*max_points)
      CALL POPREAL8ARRAY(point%qm, 2*4*max_points)
      CALL POPREAL8ARRAY(point%delta, max_points)
      CALL FPI_SOLVER_D_B(it)
      cmdb = 0.0_8
      cddb = 0.0_8
      cldb = 0.0_8
    END DO
    CALL COMPUTE_NORMALS_B()
  END SUBROUTINE Q_LSKUM_D_B

!  Differentiation of q_lskum in forward (tangent) mode (with options fixinterface):
!   variations   of useful results: cd cl cm
!   with respect to varying inputs: aoa
!   RW status of diff variables: q_inf:(loc) theta:(loc) aoa:in
!                q_init:(loc) cd:out cl:out cm:out point.prim:(loc)
!                point.prim_old:(loc) point.flux_res:(loc) point.q:(loc)
!                point.dq:(loc) point.qm:(loc) point.delta:(loc)
  SUBROUTINE Q_LSKUM_D()
    IMPLICIT NONE
    INTEGER :: i
    CALL INITIAL_CONDITIONS_D()
    CALL COMPUTE_NORMALS()
    CALL GENERATE_CONNECTIVITY()
! Set U_old to U for first iteration
    DO i=1,local_points
      point%u_old(1, i) = point%prim(1, i)
      point%u_old(2, i) = point%prim(1, i)*point%prim(2, i)
      point%u_old(3, i) = point%prim(1, i)*point%prim(3, i)
      point%u_old(4, i) = 2.5d0*point%prim(4, i) + 0.5d0*point%prim(1, i&
&       )*(point%prim(2, i)*point%prim(2, i)+point%prim(3, i)*point%prim&
&       (3, i))
    END DO
    IF (restart .EQ. 0) THEN
      itr = 0
      cdd = 0.0_8
      cld = 0.0_8
      cmd = 0.0_8
      pointd%prim_old = 0.0_8
      pointd%flux_res = 0.0_8
      pointd%q = 0.0_8
      pointd%dq = 0.0_8
      pointd%qm = 0.0_8
      pointd%delta = 0.0_8
    ELSE
      cdd = 0.0_8
      cld = 0.0_8
      cmd = 0.0_8
      pointd%prim_old = 0.0_8
      pointd%flux_res = 0.0_8
      pointd%q = 0.0_8
      pointd%dq = 0.0_8
      pointd%qm = 0.0_8
      pointd%delta = 0.0_8
    END IF
    DO it=itr+1,itr+max_iters
      CALL FPI_SOLVER_D(it)
    END DO
  END SUBROUTINE Q_LSKUM_D

  SUBROUTINE Q_LSKUM()
    IMPLICIT NONE
    INTEGER :: i
    CALL INITIAL_CONDITIONS()
    CALL COMPUTE_NORMALS()
    CALL GENERATE_CONNECTIVITY()
! Set U_old to U for first iteration
    DO i=1,local_points
      point%u_old(1, i) = point%prim(1, i)
      point%u_old(2, i) = point%prim(1, i)*point%prim(2, i)
      point%u_old(3, i) = point%prim(1, i)*point%prim(3, i)
      point%u_old(4, i) = 2.5d0*point%prim(4, i) + 0.5d0*point%prim(1, i&
&       )*(point%prim(2, i)*point%prim(2, i)+point%prim(3, i)*point%prim&
&       (3, i))
    END DO
    IF (restart .EQ. 0) itr = 0
    DO it=itr+1,itr+max_iters
      CALL FPI_SOLVER(it)
    END DO
  END SUBROUTINE Q_LSKUM

END MODULE Q_LSKUM_MOD_DIFF_DIFF

