!#########################################################
! RAISON D'ÊTRE:
! - Implementation of those routines needed to define a CRP PES 
! UPDATES:
! - Created: 18/12/2013 ---> Alberto Muzas 
! FUNCTIONALITY:
! - Type definition 
! IDEAS FOR THE FUTURE:
! - None 
!##########################################################
MODULE CRP_MOD
   USE PES_MOD
! Initial declarations
IMPLICIT NONE
TYPE, EXTENDS(PES) :: CRP
END TYPE CRP
END MODULE CRP_MOD
