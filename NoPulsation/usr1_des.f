!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvv!
!                                                                      !
!  Module name: URS1_DES                                               !
!                                                                      !
!  Purpose: This routine is called within the discrete phase time loop !
!  after the source terms have been calculated but before they are     !
!  applied. The user may insert code in this routine or call user      !
!  defined subroutines.                                                !
!                                                                      !
!  This routine is called from the time loop, but no indices (fluid    !
!  cell or particle) are defined.                                      !
!                                                                      !
!  Author: J.Musser                                   Date: 06-Nov-12  !
!                                                                      !
!  Comments:                                                           !
!                                                                      !
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^!
      SUBROUTINE USR1_DES
      
      USE constant, only : gravity, gravity_x, gravity_y, gravity_z, Pi
      USE run, only: time, dt, tstop
      USE discretelement, only : grav, FC, TOW, MAX_PIP, PIJK, DES_VEL_NEW, DES_USR_VAR, F_drag
      USE bc, only: BC_V_g
      USE set_bc0_flow_mod, only: set_bc0_flow
      IMPLICIT NONE

      INTEGER :: NP, IJK

      GRAVITY_X=0
      GRAVITY_Z=0
      GRAVITY_Y = -9.81
      GRAV(1) = GRAVITY_X
      GRAV(2) = GRAVITY_Y
      GRAV(3) = GRAVITY_Z
     
      ! Left inlet velocity
      bc_v_g(3) = 0.32 
      ! Right inlet velocity
      bc_v_g(5) = 0.32 

      ! Spout inlet velocity
      bc_v_g(4) = 11.91 

      call set_bc0_flow

      DO NP = 1,MAX_PIP
            DES_USR_VAR(1,NP)=F_drag(NP,1)
            DES_USR_VAR(2,NP)=F_drag(NP,2)
            DES_USR_VAR(3,NP)=F_drag(NP,3)				  
      ENDDO	


      RETURN
      END SUBROUTINE USR1_DES
