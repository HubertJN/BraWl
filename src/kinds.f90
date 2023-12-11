!----------------------------------------------------------------------!
! kinds.f90                                                            !
!                                                                      !
! Module for defining f2008 KINDS.                                     !
!                                                                      !
! C. D. Woodgate,  Warwick                                        2023 !
!----------------------------------------------------------------------!
module kinds

  implicit none

  integer, parameter :: int8 = selected_int_kind(2)
  integer, parameter :: int16 = selected_int_kind(4)
  integer, parameter :: int32 = selected_int_kind(9)
  integer, parameter :: int64 = selected_int_kind(15)
  integer, parameter :: real32 = selected_real_kind(6, 37)
  integer, parameter :: real64 = selected_real_kind(15, 307)
  integer, parameter :: real128 = selected_real_kind(33, 4931)

end module
