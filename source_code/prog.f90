!+--------------------------------------------------------+
!|                      PADDI-Code                        |
!|                        =====                           |
!|          (Pa)rallel (D)oube (Di)ffusion Code           |
!|                                                        |
!| Spectral Code to solve the double diffusive Boussinesq |
!| equations in a triply periodic cube.                   |
!+--------------------------------------------------------+
!| Written by:    Stephan Stellmach         (2008-2011)   |
!|                (stellma@uni-muenster.de)               |
!|                                                        |
!| To prevent overlap, please always contact me before    |
!| using this code on a particular problem!               |
!| For personal use only. Do not distribute!              |
!+--------------------------------------------------------+
PROGRAM double_diffusion
  USE defprecision_module
  USE message_passing_module, ONLY: start_mpi,stop_mpi,myid
  USE main_module, ONLY: init,read_parameter,free_allocated_memory, &
      &                  timestep_AB_BDF3,timestep_RK2_CN, allocate_drag, &
      &                  compute_drag, deallocate_drag 
  USE state_module, ONLY: velocity, buoyancy
  USE parameter_module, ONLY: nsteps
  USE IO_module, ONLY: open_files,write_compressed_file,write_diagnostics_file,close_files, &
      &                write_output_files
  USE message_passing_module, ONLY : myid
#ifdef MPI_MODULE
  USE MPI
  IMPLICIT NONE
#else
  IMPLICIT NONE
  INCLUDE "mpif.h"
#endif
  CHARACTER :: restarted
  TYPE(velocity) :: u,up
  TYPE(buoyancy) :: Temp,Chem,Part
  INTEGER(kind=ki) :: istep,startstep
  REAL(kind=kr) :: t,dt,dt1,dt2
  REAL(kind=krd)     :: time_begin=0.,time_end=0.
  REAL(kind=kr),POINTER :: drag(:,:,:,:) 
  
! Initialize MPI
  CALL start_mpi


! Read parameter values for the run
  CALL read_parameter(restarted)


! perform necessary initializations
  CALL init(u,Temp,Chem,up,Part,t,dt,dt1,dt2,startstep,restarted)
#ifdef PARTICLE_FIELD
  CALL allocate_drag(drag) 
#endif

!  write(*,*) myid,M_init

! open output files 
  CALL open_files
  CALL write_diagnostics_file(u,Temp,Chem,up,Part,0,t,dt)
  CALL write_compressed_file(u,Temp,Chem,up,Part,0,t,dt)
   
! start time stepping

#ifdef AB_BDF3
  ! Third order Adams-Bashforth / Backward-Differencing multi-step method
  time_stepping_loop : DO istep = startstep+1,startstep+nsteps
     IF (myid.EQ.0) WRITE(*,'(a,I7,a,F8.3,a)') "step :",istep,"  ",time_end-time_begin," s"
     time_begin=MPI_WTIME()
     IF (istep.LT.startstep+3) THEN
        ! use second order Runge-Kutta / Crank Nicholson as starting scheme
#ifdef PARTICLE_FIELD
        CALL compute_drag(u,up,drag)
#else  
        drag = 0._kr
#endif
        CALL timestep_RK2_CN(u,Temp,Chem,up,Part,drag,istep,t,dt,dt1,dt2)
     ELSE 
#ifdef PARTICLE_FIELD
        CALL compute_drag(u,up,drag)
#else  
        drag = 0._kr
#endif
        CALL timestep_AB_BDF3(u,Temp,Chem,up,Part,drag,istep,t,dt,dt1,dt2)
     ENDIF
     time_end=MPI_WTIME()
     CALL write_output_files(u,Temp,Chem,up,Part,t,dt,istep)
  END DO time_stepping_loop
#endif

#ifdef PARTICLE_FIELD
  CALL deallocate_drag(drag)
#endif

! close output files
  CALL close_files

! free allocated memory
  CALL free_allocated_memory(u,Temp,Chem,up,Part)

! Stop MPI
  CALL stop_mpi

END PROGRAM double_diffusion
