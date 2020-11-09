! This subroutine computes the drag term. 
SUBROUTINE compute_drag(u,up,drag)
  USE defprecision_module
  USE state_module, ONLY: velocity
  !USE state_module, ONLY: buoyancy,velocity,ltime0,ltime1,ltime2,rtime1,rtime2
  USE mpi_transf_module, ONLY: mysx_spec,myex_spec,mysy_spec,myey_spec, &
      &                        mysy_phys,myey_phys,mysz_phys,myez_phys
  USE parameter_module, ONLY : Nx,T_part
  IMPLICIT NONE
  !TYPE(buoyancy)   :: Part ! Temperature and Chemical field
  TYPE(velocity)   :: u,up ! velocity field

  REAL(kind=kr) :: drag(0:Nx-1,mysy_phys:myey_phys,mysz_phys:myez_phys,dim_vec)
 !Is the best way to set 0 in prog or here?  
  drag = (u%phys(:,:,:,:)- up%phys(:,:,:,:))/T_part
  

END SUBROUTINE compute_drag
