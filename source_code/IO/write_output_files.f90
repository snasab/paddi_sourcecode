SUBROUTINE write_output_files(u,Temp,Chem,up,Part,t,dt,istep)
  USE defprecision_module
  USE state_module, ONLY: velocity,buoyancy
  USE pnetCDF_IO_module, ONLY: pn_write_dump,pn_write_step_simdat_file
  IMPLICIT NONE
  TYPE(velocity) :: u,up
  TYPE(buoyancy) :: Temp,Chem,Part
  REAL(kind=kr) :: t,dt
  INTEGER(kind=ki) :: istep

  IF (MOD(istep,n_comp_diag)==0)             &       !Write out diagnostics
       &   CALL write_diagnostics_file(u,Temp,Chem,up,Part,istep,t,dt)
  IF (MOD(istep,n_wrt_spec)==0)            &         !Write out horizontal spectra
       &   CALL  write_horizontal_spectra(u,Temp,Chem,up,Part,istep,t)
  IF (MOD(istep,n_wrt_spec)==0)            &          !Write out vertical spectra
       &   CALL  write_vertical_spectra(u,Temp,Chem,up,Part,istep,t)
  IF (MOD(istep,n_wrt_prof)==0)            &          !Write out vertical profiles
       &   CALL  write_z_profile(u,Temp,Chem,up,Part,istep,t)
  IF (MOD(istep,n_wrt_jc)==0)                &       !Write compressed file  
       &   CALL write_compressed_file(u,Temp,Chem,up,Part,istep,t,dt) 
  IF (MOD(istep,n_wrt_dump)==0)              &       !Write restart file (netCDF)
       &   CALL pn_write_dump(u,Temp,Chem,up,Part,t,dt,istep)
  IF (MOD(istep,n_wrt_netCDF)==0)            &       !Write simulation data file (netCDF)
       &   CALL pn_write_step_simdat_file(u,Temp,Chem,up,Part,istep,t,dt)


END SUBROUTINE write_output_files
