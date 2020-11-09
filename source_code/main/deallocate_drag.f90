subroutine deallocate_drag(drag)
  USE mpi_transf_module, ONLY: mysx_spec,myex_spec,mysy_spec,myey_spec, &
      &                        mysy_phys,myey_phys,mysz_phys,myez_phys
  USE parameter_module, ONLY : Nx,T_part
  IMPLICIT NONE
  REAL(kind=kr), POINTER :: drag(:,:,:,:)
  !REAL(kind=kr) :: drag(0:Nx-1,mysy_phys:myey_phys,mysz_phys:myez_phys,dim_vec)

#ifdef PARTICLE_FIELD
  DEALLOCATE(drag)
#endif

end subroutine deallocate_drag
