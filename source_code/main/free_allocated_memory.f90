SUBROUTINE free_allocated_memory(u,Temp,Chem)
  USE defprecision_module
  USE state_module, ONLY: velocity,buoyancy,deallocate_uTC
  USE mpi_transf_module, ONLY: free_transforms
  USE parameter_module, ONLY: deallocate_fft_storage_scheme
  IMPLICIT NONE
  TYPE(velocity) :: u
  TYPE(buoyancy) :: Temp,Chem

  CALL free_transforms
  CALL deallocate_fft_storage_scheme
  CALL deallocate_uTC(u,Temp,Chem)

END SUBROUTINE free_allocated_memory
