! This subroutine computes the drag term. 
SUBROUTINE allocate_drag(drag)
  USE defprecision_module
  USE defs_2D_3D_module
  use parameter_module, ONLY: Nx
  USE mpi_transf_module, ONLY: mysy_phys,myey_phys, &
      &                        mysz_phys,myez_phys
  IMPLICIT NONE
  REAL(kind=kr),POINTER :: drag(:,:,:,:)       
  
  ALLOCATE(drag(0:Nx-1,mysy_phys:myey_phys,mysz_phys:myez_phys,dim_vec)) 

END SUBROUTINE allocate_drag
