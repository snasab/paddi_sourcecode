! This subroutine performs a time step by the second order RK2/CN method.
! The one-step algorithm is used as a starting scheme for the multistep AB/BDI3 method.
! Theory tells us that a second order starting scheme is sufficient to guarantee that the 
! overall method is third order accurate. 
! The RK2/CN algorithm exists in various variants. We use a combination of the Crank 
! Nicholson method with the "modified Euler version" of the 2nd order Runge Kutta method 
! here.
! For an ODE of type 
!
!   d phi / dt = L(phi) + N(phi)
!
! where L is linear and N nonlinear, the method is defined as
!
!   phi1 := phi^n + dt ( N(phi^n) + L(phi^n) )
!   phi^(n+1) - 0.5 * dt * L(phi^(n+1)) = 0.5 * (phi1 + phi^n) + 0.5 * dt N(phi^n)
!
! where the superscripts denote time levels. 
! 
SUBROUTINE timestep_RK2_CN(u,Temp,Chem,up,Part,drag,istep,t,dt,dt1,dt2)  
  USE defprecision_module
  USE state_module, ONLY: velocity,buoyancy,shift_time_pointers,ltime0,ltime1,rtime1
  USE parameter_module, ONLY : Nx,D_visc,D_therm,D_comp,S_therm,S_comp,T_part,G_part,&
                               & D_part,S_part,Nmax,Dv_part
  USE IO_module, ONLY: write_compressed_file,write_diagnostics_file,n_comp_diag
  USE mpi_transf_module, ONLY: mysy_phys,myey_phys,mysz_phys,myez_phys, &
      &                        mysx_spec,myex_spec,mysy_spec,myey_spec
  IMPLICIT NONE
  TYPE(velocity) :: u,up
  TYPE(buoyancy) :: Temp,Chem,Part
  INTEGER(kind=ki) :: istep
  REAL(kind=kr) :: t                ! dimensionless time
  REAL(kind=kr) :: dt,dt1,dt2       ! last three time step sizes 
  COMPLEX(kind=kr),POINTER :: work_scalar(:,:,:),work_vector(:,:,:,:)
  REAL(kind=kr),POINTER :: drag(:,:,:,:) ! RHS 

  ! shift time pointers (t_{n-1} -> t_{n-2}, t_n -> t_{n-1})
  CALL shift_time_pointers
  dt2=dt1
  dt1=dt

  ! time step size control
  dt = dt !we do not change the time step during startup
  t=t+dt

  ! Compute RHS at timelevel n-1 and store for future use in the
  ! third order Adams-Bashforth part
#ifdef TEMPERATURE_FIELD
  CALL crhs_buoyancy(Temp%rhs(:,:,:,rtime1),Temp,u,ltime1,S_therm)
#endif
#ifdef CHEMICAL_FIELD
  CALL crhs_buoyancy(Chem%rhs(:,:,:,rtime1),Chem,u,ltime1,S_comp)
#endif
#ifdef PARTICLE_FIELD
  CALL crhs_particle(Part%rhs(:,:,:,rtime1),Part,up,ltime1)
  CALL crhs_part_velocity(up%rhs(:,:,:,:,rtime1),u,up,Part,drag)
#endif
  CALL crhs_velocity(u%rhs(:,:,:,:,rtime1),u,Temp,Chem,up,Part,drag)
  

  ! perform first step of the RK2/CN method
  ! save the intermediate values phi1 in the array location reserved for the values 
  ! at $t_n$ (new time level ltime0)
#ifdef TEMPERATURE_FIELD
  CALL step1_RK2_CN(Temp%spec(:,:,:,ltime0),Temp%spec(:,:,:,ltime1), &
       &           Temp%rhs(:,:,:,rtime1),dt,D_therm)
#endif
#ifdef CHEMICAL_FIELD
  CALL step1_RK2_CN(Chem%spec(:,:,:,ltime0),Chem%spec(:,:,:,ltime1), &
       &           Chem%rhs(:,:,:,rtime1),dt,D_comp)
#endif
#ifdef PARTICLE_FIELD
  CALL step1_RK2_CN(Part%spec(:,:,:,ltime0),Part%spec(:,:,:,ltime1), &
       &           Part%rhs(:,:,:,rtime1),dt,D_part)
  CALL step1_RK2_CN(up%spec(:,:,:,vec_x,ltime0),up%spec(:,:,:,vec_x,ltime1),       &
	   &           up%rhs(:,:,:,vec_x,rtime1),dt,Dv_part)
#ifndef TWO_DIMENSIONAL
  CALL step1_RK2_CN(up%spec(:,:,:,vec_y,ltime0),up%spec(:,:,:,vec_y,ltime1),       &
	   &           up%rhs(:,:,:,vec_y,rtime1),dt,Dv_part)
#endif
  CALL step1_RK2_CN(up%spec(:,:,:,vec_z,ltime0),up%spec(:,:,:,vec_z,ltime1),       &
	   &           up%rhs(:,:,:,vec_z,rtime1),dt,Dv_part)
#endif
  CALL step1_RK2_CN(u%spec(:,:,:,vec_x,ltime0),u%spec(:,:,:,vec_x,ltime1),         &
                  &                             u%rhs(:,:,:,vec_x,rtime1),dt,D_visc)
#ifndef TWO_DIMENSIONAL
  CALL step1_RK2_CN(u%spec(:,:,:,vec_y,ltime0),u%spec(:,:,:,vec_y,ltime1),         &
                  &                             u%rhs(:,:,:,vec_y,rtime1),dt,D_visc)
#endif
  CALL step1_RK2_CN(u%spec(:,:,:,vec_z,ltime0),u%spec(:,:,:,vec_z,ltime1),         &
                  &                             u%rhs(:,:,:,vec_z,rtime1),dt,D_visc)
  
  ! convert intermediate Temp,Chem,u fields to physical space
  CALL compute_phys_space_vars(u,Temp,Chem,up,Part)

  ! compute RHS at intermediate time level and solve for new state
  ! a new array is needed to hold this information
#ifdef TEMPERATURE_FIELD
  ALLOCATE(work_scalar(0:2*Nmax-1,mysx_spec:myex_spec,mysy_spec:myey_spec))
  CALL crhs_buoyancy(work_scalar,Temp,u,ltime0,S_therm)
  CALL step2_RK2_CN(Temp%spec(:,:,:,ltime0),Temp%spec(:,:,:,ltime1), &
       &           work_scalar,dt,D_therm)
  DEALLOCATE(work_scalar)
#endif
#ifdef CHEMICAL_FIELD
  ALLOCATE(work_scalar(0:2*Nmax-1,mysx_spec:myex_spec,mysy_spec:myey_spec))
  CALL crhs_buoyancy(work_scalar,Chem,u,ltime0,S_comp)
  CALL step2_RK2_CN(Chem%spec(:,:,:,ltime0),Chem%spec(:,:,:,ltime1), &
       &           work_scalar,dt,D_comp)
  DEALLOCATE(work_scalar)  
#endif
#ifdef PARTICLE_FIELD
  ALLOCATE(work_scalar(0:2*Nmax-1,mysx_spec:myex_spec,mysy_spec:myey_spec))
  CALL crhs_particle(work_scalar,Part,up,ltime0)
  CALL step2_RK2_CN(Part%spec(:,:,:,ltime0),Part%spec(:,:,:,ltime1), &
       &           work_scalar,dt,D_part)
  DEALLOCATE(work_scalar)  
  ALLOCATE(work_vector(0:2*Nmax-1,mysx_spec:myex_spec,mysy_spec:myey_spec,3)) 
  CALL crhs_part_velocity(work_vector,u,up,Part,drag)
  CALL step2_RK2_CN(up%spec(:,:,:,vec_x,ltime0),up%spec(:,:,:,vec_x,ltime1),work_vector(:,:,:,vec_x),dt,Dv_part)
#ifndef TWO_DIMENSIONAL
  CALL step2_RK2_CN(up%spec(:,:,:,vec_y,ltime0),up%spec(:,:,:,vec_y,ltime1),work_vector(:,:,:,vec_y),dt,Dv_part)
#endif
  CALL step2_RK2_CN(up%spec(:,:,:,vec_z,ltime0),up%spec(:,:,:,vec_z,ltime1),work_vector(:,:,:,vec_z),dt,Dv_part)
  DEALLOCATE(work_vector)
#endif

  ALLOCATE(work_vector(0:2*Nmax-1,mysx_spec:myex_spec,mysy_spec:myey_spec,3)) 
  CALL crhs_velocity(work_vector,u,Temp,Chem,up,Part,drag)
  CALL step2_RK2_CN(u%spec(:,:,:,vec_x,ltime0),u%spec(:,:,:,vec_x,ltime1),work_vector(:,:,:,vec_x),dt,D_visc)
#ifndef TWO_DIMENSIONAL
  CALL step2_RK2_CN(u%spec(:,:,:,vec_y,ltime0),u%spec(:,:,:,vec_y,ltime1),work_vector(:,:,:,vec_y),dt,D_visc)
#endif
  CALL step2_RK2_CN(u%spec(:,:,:,vec_z,ltime0),u%spec(:,:,:,vec_z,ltime1),work_vector(:,:,:,vec_z),dt,D_visc)
  DEALLOCATE(work_vector)
  ! convert new Temp,Chem,u fields to physical space
!#ifdef PARTICLE_FIELD
! For now, just copy solutions at t_{n-1} into solution at t_{n}
 ! up%spec(:,:,:,vec_x,ltime0) = up%spec(:,:,:,vec_x,ltime1)
  !up%spec(:,:,:,vec_y,ltime0) = up%spec(:,:,:,vec_y,ltime1)
  !up%spec(:,:,:,vec_z,ltime0) = up%spec(:,:,:,vec_z,ltime1)
!#endif
  CALL compute_phys_space_vars(u,Temp,Chem,up,Part) 

END SUBROUTINE timestep_RK2_CN
