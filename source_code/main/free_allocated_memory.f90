SUBROUTINE free_allocated_memory(u,Temp,Chem,up,Part)
  USE defprecision_module
  USE state_module, ONLY: velocity,buoyancy,deallocate_uTC
  USE mpi_transf_module, ONLY: free_transforms
  USE parameter_module, ONLY: deallocate_fft_storage_scheme
  IMPLICIT NONE
  TYPE(velocity) :: u,up
  TYPE(buoyancy) :: Temp,Chem,Part

  CALL free_transforms
  CALL deallocate_fft_storage_scheme
  CALL deallocate_uTC(u,Temp,Chem,up,Part)

END SUBROUTINE free_allocated_memory
