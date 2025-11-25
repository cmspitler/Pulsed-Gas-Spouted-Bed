MODULE CFNEWVALUES_MOD
      USE SQ_OBB_MOD
      USE SQ_CONTACT_NEWTON_DPmethod_MOD
      USE SQ_ROTATION_MOD
      USE SQ_PROPERTIES_MOD
      USE Sq_math_mod
      USE SQ_CONTACT_WALL

CONTAINS
!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: CFNEWVALUES                                            C
!
!  Purpose: DES - Calculate the new values of particle velocity,
!           position, angular velocity etc
!
!                                                                      C
!  Comments: Implements Eqns 1, 2, 3, 4 & 5  from the following paper:
!    Tsuji Y., Kawaguchi T., and Tanak T., "Lagrangian numerical
!    simulation of plug glow of cohesionless particles in a
!    horizontal pipe", Powder technology, 71, 239-250, 1992
!
!  pradeep : changes for parallel processing
!          1. periodic boundaries might lie in different proc. so adjust
!             particle position for periodic removed
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
   SUBROUTINE CFNEWVALUES

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE constant
      USE des_bc
      USE discretelement
      USE fldvar
      USE mfix_pic
      USE mpi_utility
      USE parallel
      USE run
      USE param
      USE param1
      USE physprop
      use geometry, only: DO_K, NO_K

!**************************************qg***************************************************
      use discretelement, only: F_drag
      use randomno	 
!**************************************qg***************************************************             

      IMPLICIT NONE
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
      INTEGER :: L
      LOGICAL, SAVE :: FIRST_PASS = .TRUE.
      DOUBLE PRECISION :: OMEGA_MAG,OMEGA_UNIT(3),ROT_ANGLE

!-----------------------------------------------

! Adams-Bashforth defaults to Euler for the first time step.
      IF(FIRST_PASS .AND. INTG_ADAMS_BASHFORTH) THEN
         DO L =1, MAX_PIP
            IF(IS_NONEXISTENT(L)) CYCLE                       ! Only real particles
            IF(IS_ENTERING(L).or.IS_ENTERING_GHOST(L)) CYCLE  ! Only non-entering
            IF(IS_GHOST(L)) CYCLE                             ! Skip ghost particles
            DES_ACC_OLD(L,:) = FC(L,:)/PMASS(L) + GRAV(:)
            ROT_ACC_OLD(L,:) = TOW(L,:)
         ENDDO
      ENDIF


!!$omp parallel default(none)                    &
!!$omp shared(MAX_PIP,INTG_EULER,INTG_ADAMS_BASHFORTH,fc,tow,do_nsearch,   &
!!$omp       omega_new,omega_old,pmass,grav,des_vel_new,des_pos_new,       &
!!$omp       des_vel_old,des_pos_old,dtsolid,omoi,des_acc_old,rot_acc_old, &
!!$omp       ppos,neighbor_search_rad_ratio,des_radius,DO_OLD, iGlobal_ID, &
!!$omp       residence_time, &
!!$omp       particle_orientation,orientation,particle_state) &
!!$omp private(l,rot_angle,omega_mag,omega_unit,aabb)

! If a particle is classified as new, then forces are ignored.
! Classification from new to existing is performed in routine
! des_check_new_particle.f

! Advance particle position, velocity
! first-order method
      IF (INTG_EULER) THEN
!!$omp sections
!!$omp section
         WHERE(PARTICLE_STATE == NORMAL_PARTICLE) DES_VEL_NEW(:,1) =   &
            DES_VEL_NEW(:,1) + DTSOLID*(FC(:,1)/PMASS(:) + GRAV(1))
         WHERE(PARTICLE_STATE < NORMAL_GHOST) DES_POS_NEW(:,1) =       &
            DES_POS_NEW(:,1) + DES_VEL_NEW(:,1)*DTSOLID
         FC(:,1) = ZERO

!!$omp section
         WHERE(PARTICLE_STATE == NORMAL_PARTICLE) DES_VEL_NEW(:,2) =   &
            DES_VEL_NEW(:,2) + DTSOLID*(FC(:,2)/PMASS(:) + GRAV(2))
         WHERE(PARTICLE_STATE < NORMAL_GHOST) DES_POS_NEW(:,2) =       &
            DES_POS_NEW(:,2) + DES_VEL_NEW(:,2)*DTSOLID
         FC(:,2) = ZERO

!!$omp section
         WHERE(PARTICLE_STATE == NORMAL_PARTICLE) DES_VEL_NEW(:,3) =   &
            DES_VEL_NEW(:,3) + DTSOLID*(FC(:,3)/PMASS(:) + GRAV(3))
         WHERE(PARTICLE_STATE < NORMAL_GHOST) DES_POS_NEW(:,3) =       &
            DES_POS_NEW(:,3) + DES_VEL_NEW(:,3)*DTSOLID
         FC(:,3) = ZERO

!!$omp section
         WHERE(PARTICLE_STATE == NORMAL_PARTICLE) OMEGA_NEW(:,1) =     &
            OMEGA_NEW(:,1) + TOW(:,1)*OMOI(:)*DTSOLID
         TOW(:,1) = ZERO

!!$omp section
         WHERE(PARTICLE_STATE == NORMAL_PARTICLE)OMEGA_NEW(:,2) =      &
            OMEGA_NEW(:,2) + TOW(:,2)*OMOI(:)*DTSOLID
         TOW(:,2) = ZERO

!!$omp section
         WHERE(PARTICLE_STATE == NORMAL_PARTICLE)OMEGA_NEW(:,3) =      &
            OMEGA_NEW(:,3) + TOW(:,3)*OMOI(:)*DTSOLID
         TOW(:,3) = ZERO
!!$omp end sections

! Second-order Adams-Bashforth/Trapezoidal scheme
      ELSEIF (INTG_ADAMS_BASHFORTH) THEN

!!$omp sections
!!$omp section
         WHERE(PARTICLE_STATE == NORMAL_PARTICLE) &
            FC(:MAX_PIP,1) = FC(:MAX_PIP,1)/PMASS(:MAX_PIP) + GRAV(1)
         WHERE(PARTICLE_STATE == NORMAL_PARTICLE) &
            DES_VEL_NEW(:MAX_PIP,1) = DES_VEL_OLD(:MAX_PIP,1) + 0.5d0* &
               (3.d0*FC(:MAX_PIP,1)-DES_ACC_OLD(:MAX_PIP,1) )*DTSOLID
            DES_ACC_OLD(:MAX_PIP,1) = FC(:MAX_PIP,1)

         WHERE(PARTICLE_STATE < NORMAL_GHOST)                          &
            DES_POS_NEW(:MAX_PIP,1) = DES_POS_OLD(:MAX_PIP,1) + 0.5d0* &
               (DES_VEL_OLD(:MAX_PIP,1)+DES_VEL_NEW(:MAX_PIP,1))*DTSOLID
         FC(:MAX_PIP,1) = ZERO

!!$omp section
         WHERE(PARTICLE_STATE == NORMAL_PARTICLE) &
            FC(:MAX_PIP,2) = FC(:MAX_PIP,2)/PMASS(:MAX_PIP) + GRAV(2)
         WHERE(PARTICLE_STATE == NORMAL_PARTICLE) &
            DES_VEL_NEW(:MAX_PIP,2) = DES_VEL_OLD(:MAX_PIP,2) + 0.5d0* &
               (3.d0*FC(:MAX_PIP,2)-DES_ACC_OLD(:MAX_PIP,2) )*DTSOLID
            DES_ACC_OLD(:MAX_PIP,2) = FC(:MAX_PIP,2)

         WHERE(PARTICLE_STATE < NORMAL_GHOST)                          &
            DES_POS_NEW(:MAX_PIP,2) = DES_POS_OLD(:MAX_PIP,2) + 0.5d0* &
               (DES_VEL_OLD(:MAX_PIP,2)+DES_VEL_NEW(:MAX_PIP,2))*DTSOLID
         FC(:MAX_PIP,2) = ZERO

!!$omp section
         WHERE(PARTICLE_STATE == NORMAL_PARTICLE) &
            FC(:MAX_PIP,3) = FC(:MAX_PIP,3)/PMASS(:MAX_PIP) + GRAV(3)
         WHERE(PARTICLE_STATE == NORMAL_PARTICLE) &
            DES_VEL_NEW(:MAX_PIP,3) = DES_VEL_OLD(:MAX_PIP,3) + 0.5d0* &
                 (3.d0*FC(:MAX_PIP,3)-DES_ACC_OLD(:MAX_PIP,3) )*DTSOLID
            DES_ACC_OLD(:MAX_PIP,3) = FC(:MAX_PIP,3)

         WHERE(PARTICLE_STATE < NORMAL_GHOST)                          &
            DES_POS_NEW(:MAX_PIP,3) = DES_POS_OLD(:MAX_PIP,3) + 0.5d0* &
               (DES_VEL_OLD(:MAX_PIP,3)+DES_VEL_NEW(:MAX_PIP,3))*DTSOLID
         FC(:MAX_PIP,3) = ZERO

!!$omp section
        WHERE(PARTICLE_STATE(:MAX_PIP) == NORMAL_PARTICLE) &
           OMEGA_NEW(:MAX_PIP,1) = OMEGA_OLD(:MAX_PIP,1) + 0.5d0*     &
               (3.d0*TOW(:MAX_PIP,1)*OMOI(:MAX_PIP) -                  &
               ROT_ACC_OLD(:MAX_PIP,1))*DTSOLID
        WHERE(PARTICLE_STATE(:MAX_PIP) == NORMAL_PARTICLE) &
            ROT_ACC_OLD(:MAX_PIP,1) = TOW(:MAX_PIP,1)*OMOI(:MAX_PIP)
         TOW(:MAX_PIP,1) = ZERO

!!$omp section
         WHERE(PARTICLE_STATE == NORMAL_PARTICLE) &
            OMEGA_NEW(:MAX_PIP,2) = OMEGA_OLD(:MAX_PIP,2) + 0.5d0*     &
                 (3.d0*TOW(:MAX_PIP,2)*OMOI(:MAX_PIP)-                 &
                 ROT_ACC_OLD(:MAX_PIP,2) )*DTSOLID
         WHERE(PARTICLE_STATE == NORMAL_PARTICLE) &
            ROT_ACC_OLD(:MAX_PIP,2) = TOW(:MAX_PIP,2)*OMOI(:MAX_PIP)
         TOW(:MAX_PIP,2) = ZERO

!!$omp section
        WHERE(PARTICLE_STATE == NORMAL_PARTICLE) &
            OMEGA_NEW(:MAX_PIP,3) = OMEGA_OLD(:MAX_PIP,3) + 0.5d0*     &
               (3.d0*TOW(:MAX_PIP,3)*OMOI(:MAX_PIP)-                   &
               ROT_ACC_OLD(:MAX_PIP,3) )*DTSOLID
        WHERE(PARTICLE_STATE == NORMAL_PARTICLE) &
            ROT_ACC_OLD(:MAX_PIP,3) = TOW(:MAX_PIP,3)*OMOI(:MAX_PIP)
         TOW(:MAX_PIP,3) = ZERO

!!$omp end sections
      ENDIF

! Update particle orientation - Always first order
! When omega is non-zero, compute the rotation angle, and apply the
! Rodrigues' rotation formula

      IF(PARTICLE_ORIENTATION) THEN
         DO L = 1, MAX_PIP
            OMEGA_MAG = OMEGA_NEW(L,1)**2 + OMEGA_NEW(L,2)**2 + OMEGA_NEW(L,3)**2

            IF(OMEGA_MAG>ZERO) THEN
               OMEGA_MAG=DSQRT(OMEGA_MAG)
               OMEGA_UNIT(:) = OMEGA_NEW(L,:)/OMEGA_MAG
               ROT_ANGLE = OMEGA_MAG * DTSOLID

               ORIENTATION(L,:) = ORIENTATION(L,:)*DCOS(ROT_ANGLE) + &
                  DES_CROSSPRDCT(OMEGA_UNIT,ORIENTATION(L,:))*DSIN(ROT_ANGLE) + &
                  OMEGA_UNIT(:)*DOT_PRODUCT(OMEGA_UNIT,ORIENTATION(L,:))*&
                  (ONE-DCOS(ROT_ANGLE))
            ENDIF
         ENDDO
      ENDIF

! Residence time
!!$omp section
         WHERE(PARTICLE_STATE == NORMAL_PARTICLE) RESIDENCE_TIME(:) = RESIDENCE_TIME(:) + DTSOLID
!!$omp end sections

! Check if the particle has moved a distance greater than or equal to
! its radius since the last time a neighbor search was called. if so,
! make sure that neighbor is called in des_time_march
      IF(.NOT.DO_NSEARCH) THEN
!!$omp do reduction (.or.:do_nsearch)
         DO L = 1, MAX_PIP
            DO_NSEARCH = DO_NSEARCH .or. &
               (DES_POS_NEW(L,1) - PPOS(L,1))**2+              &
               (DES_POS_NEW(L,2) - PPOS(L,2))**2+              &
               (DES_POS_NEW(L,3) - PPOS(L,3))**2  >=           &
               (NEIGHBOR_SEARCH_RAD_RATIO*DES_RADIUS(L))**2
         ENDDO
      ENDIF

!!$omp end parallel

      FIRST_PASS = .FALSE.

      RETURN

      contains

      include 'functions.inc'

   END SUBROUTINE CFNEWVALUES

!vvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvvC
!                                                                      C
!  Module name: SuperDEM_CFNEWVALUES                                   C
!                                                                      C
!  Author: Xi Gao                                 Date: 25-Feb-2019    C
!                                                                      C
!  Purpose: SuperDEM - Calculate the new values of particle velocity,  C
!           position, angular velocity etc                             C
!                                                                      C
!                                                                      C
!^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^C
   SUBROUTINE SuperDEM_CFNEWVALUES

!-----------------------------------------------
! Modules
!-----------------------------------------------
      USE constant
      USE des_bc
      USE discretelement
      USE fldvar
      USE mfix_pic
      USE mpi_utility
      USE parallel
      USE run
      USE param
      USE param1
      USE physprop
      use geometry, only: DO_K, NO_K

      IMPLICIT NONE
!-----------------------------------------------
! Local Variables
!-----------------------------------------------
      INTEGER :: L, ixg
      LOGICAL, SAVE :: FIRST_PASS = .TRUE.
      DOUBLE PRECISION :: OMEGA_MAG,OMEGA_UNIT(3),ROT_ANGLE
      DOUBLE PRECISION :: Qc(4),Qc_new(4),DW(3),DW_local(3)
      DOUBLE PRECISION :: OMEGA_NEW_global(3),OMEGA_NEW_LOCAL(3)
      DOUBLE PRECISION :: tow_global(3),tow_local(3)
      LOGICAL, PARAMETER :: ldebug = .false.

!-----------------------------------------------

! Adams-Bashforth defaults to Euler for the first time step.
      IF(FIRST_PASS .AND. INTG_ADAMS_BASHFORTH) THEN
         DO L =1, MAX_PIP
            IF(IS_NONEXISTENT(L)) CYCLE                       ! Only real particles
            IF(IS_ENTERING(L).or.IS_ENTERING_GHOST(L)) CYCLE  ! Only non-entering
            IF(IS_GHOST(L)) CYCLE                             ! Skip ghost particles
            DES_ACC_OLD(L,:) = FC(L,:)/PMASS(L) + GRAV(:)
            ROT_ACC_OLD(L,:) = TOW(L,:)
         ENDDO
      ENDIF

! If a particle is classified as new, then forces are ignored.
! Classification from new to existing is performed in routine
! des_check_new_particle.f

! Advance particle position, velocity

         WHERE(PARTICLE_STATE == NORMAL_PARTICLE) DES_VEL_NEW(:,1) =   &
            DES_VEL_NEW(:,1) + DTSOLID*(FC(:,1)/PMASS(:) + GRAV(1))
         WHERE(PARTICLE_STATE < NORMAL_GHOST) DES_POS_NEW(:,1) =       &
            DES_POS_NEW(:,1) + DES_VEL_NEW(:,1)*DTSOLID
         FC(:,1) = ZERO

!!$omp section
         WHERE(PARTICLE_STATE == NORMAL_PARTICLE) DES_VEL_NEW(:,2) =   &
            DES_VEL_NEW(:,2) + DTSOLID*(FC(:,2)/PMASS(:) + GRAV(2))
         WHERE(PARTICLE_STATE < NORMAL_GHOST) DES_POS_NEW(:,2) =       &
            DES_POS_NEW(:,2) + DES_VEL_NEW(:,2)*DTSOLID
         FC(:,2) = ZERO

!!$omp section
         WHERE(PARTICLE_STATE == NORMAL_PARTICLE) DES_VEL_NEW(:,3) =   &
            DES_VEL_NEW(:,3) + DTSOLID*(FC(:,3)/PMASS(:) + GRAV(3))
         WHERE(PARTICLE_STATE < NORMAL_GHOST) DES_POS_NEW(:,3) =       &
            DES_POS_NEW(:,3) + DES_VEL_NEW(:,3)*DTSOLID
         FC(:,3) = ZERO


! Update of angular velocity
      IF (SuperDEM) THEN
        do L = 1, MAX_PIP
! Normal_particle
           if (PARTICLE_STATE(L) == normal_particle) then
           Qc(1)= super_q(L,1)
           Qc(2)= super_q(L,2)
           Qc(3)= super_q(L,3)
           Qc(4)= super_q(L,4)
           OMEGA_NEW_global(1)=OMEGA_NEW(L,1)
           OMEGA_NEW_global(2)=OMEGA_NEW(L,2)
           OMEGA_NEW_global(3)=OMEGA_NEW(L,3)
! Rotate angular velocity from global to local
           ixg=1
           CALL QROTATE(Qc,OMEGA_NEW_global,OMEGA_NEW_local,Ixg)
           TOW_global(1)=TOW(L,1)
           TOW_global(2)=TOW(L,2)
           TOW_global(3)=TOW(L,3)
           ixg=1
           CALL QROTATE(Qc,TOW_global,TOW_local,Ixg)

            OMEGA_NEW_LOCAL(1) = OMEGA_NEW_LOCAL(1) + &
                        ((1.0/OMOI3(L,2)-1.0/OMOI3(L,3))*OMEGA_NEW_LOCAL(2)*&
                        OMEGA_NEW_LOCAL(3)+TOW_LOCAL(1))*DTSOLID*OMOI3(L,1)

            OMEGA_NEW_LOCAL(2) = OMEGA_NEW_LOCAL(2) + &
                        ((1.0/OMOI3(L,3)-1.0/OMOI3(L,1))*OMEGA_NEW_LOCAL(3)*&
                        OMEGA_NEW_LOCAL(1)+TOW_LOCAL(2))*DTSOLID*OMOI3(L,2)

            OMEGA_NEW_LOCAL(3) = OMEGA_NEW_LOCAL(3) + &
                        ((1.0/OMOI3(L,1)-1.0/OMOI3(L,2))*OMEGA_NEW_LOCAL(1)*&
                        OMEGA_NEW_LOCAL(2)+TOW_LOCAL(3))*DTSOLID*OMOI3(L,3)
! Rotate angular velocity from local to global
            IXG=2
            CALL QROTATE(Qc,OMEGA_NEW_global,OMEGA_NEW_local,Ixg)
            OMEGA_NEW(L,1)= OMEGA_NEW_global(1)
            OMEGA_NEW(L,2)= OMEGA_NEW_global(2)
            OMEGA_NEW(L,3)= OMEGA_NEW_global(3)
            tow(L,1)=0
            tow(L,2)=0
            tow(L,3)=0
            endif
        enddo
      ENDIF !SuperDEM

!Update particle quaternion
     if (SuperDEM) then
        do L = 1, MAX_PIP
           if (PARTICLE_STATE(L) == normal_particle) then
              Qc(1)= super_q(L,1)
              Qc(2)= super_q(L,2)
              Qc(3)= super_q(L,3)
              Qc(4)= super_q(L,4)
! Incremental of angular velocity in the global frame
              DW(1)= OMEGA_NEW(L,1)*DTSOLID
              DW(2)= OMEGA_NEW(L,2)*DTSOLID
              DW(3)= OMEGA_NEW(L,3)*DTSOLID
! Apply the incremental rotation to the orientation quaternion
              ixg=1
              CALL QROTATE(Qc,DW,DW_local,Ixg)
              CALL QINCROTATE2(Qc,Qc_new,Dw_local)
              super_q(L,1)= Qc_new(1)
              super_q(L,2)= Qc_new(2)
              super_q(L,3)= Qc_new(3)
              super_q(L,4)= Qc_new(4)

              IF(PARTICLE_ORIENTATION) THEN
                 ! Convert quaternion to orientation vector
                 ! Input: quaternion Qc_new
                 !        initial orientation vector INIT_ORIENTATION
                 ! Output: orientation vector (global frame) ORIENTATION
                 ! Note the last argument is 2 to apply the rotation on the
                 ! initial orientation
                 CALL QROTATE(Qc_new,ORIENTATION(L,:),INIT_ORIENTATION,2)
              ENDIF
           endif
        ENDDO
     endif !SuperDEM

! Residence time
!!$omp section
         WHERE(PARTICLE_STATE == NORMAL_PARTICLE) RESIDENCE_TIME(:) = RESIDENCE_TIME(:) + DTSOLID
!!$omp end sections

! Check if the particle has moved a distance greater than or equal to
! its radius since the last time a neighbor search was called. if so,
! make sure that neighbor is called in des_time_march
      IF(.NOT.DO_NSEARCH) THEN
!!$omp do reduction (.or.:do_nsearch)
         DO L = 1, MAX_PIP
            DO_NSEARCH = DO_NSEARCH .or. &
               (DES_POS_NEW(L,1) - PPOS(L,1))**2+              &
               (DES_POS_NEW(L,2) - PPOS(L,2))**2+              &
               (DES_POS_NEW(L,3) - PPOS(L,3))**2  >=           &
               (NEIGHBOR_SEARCH_RAD_RATIO*DES_RADIUS(L))**2
         ENDDO
      ENDIF

      FIRST_PASS = .FALSE.

      RETURN

      contains

      include 'functions.inc'

   END SUBROUTINE SuperDEM_CFNEWVALUES

END MODULE CFNEWVALUES_MOD
