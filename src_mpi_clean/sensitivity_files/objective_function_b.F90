!        Generated by TAPENADE     (INRIA, Ecuador team)
!  Tapenade 3.16 (master) -  9 Oct 2020 17:47
!
MODULE OBJECTIVE_FUNCTION_MOD_DIFF
  USE DATA_STRUCTURE_MOD_DIFF
  USE COMPUTE_FORCE_COEFFS_MOD_DIFF
  USE COMPUTE_ENTROPY_MOD_DIFF
  USE COMPUTE_ENSTROPHY_MOD_DIFF
  USE STAGNATION_VALUES_MOD_DIFF
  IMPLICIT NONE

CONTAINS
!  Differentiation of objective_function in reverse (adjoint) mode (with options fixinterface):
!   gradient     of useful results: *cd *vector_cost_func *(point.x)
!                *(point.y) *(point.nx) *(point.ny) *(point.prim)
!   with respect to varying inputs: *cd *vector_cost_func *(point.x)
!                *(point.y) *(point.nx) *(point.ny) *(point.prim)
!   Plus diff mem management of: cd:in vector_cost_func:in point.x:in
!                point.y:in point.nx:in point.ny:in point.prim:in
  SUBROUTINE OBJECTIVE_FUNCTION_B()
    IMPLICIT NONE
    cdb = cdb + vector_cost_funcb
    vector_cost_funcb = 0.0_8
    CALL COMPUTE_CL_CD_CM_B()
  END SUBROUTINE OBJECTIVE_FUNCTION_B

  SUBROUTINE OBJECTIVE_FUNCTION()
    IMPLICIT NONE
    CALL COMPUTE_CL_CD_CM()
! call compute_entropy()
! call compute_enstrophy()
! call objective_function_J()
    vector_cost_func = cd
  END SUBROUTINE OBJECTIVE_FUNCTION

END MODULE OBJECTIVE_FUNCTION_MOD_DIFF

