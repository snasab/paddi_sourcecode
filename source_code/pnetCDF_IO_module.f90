#include "defs_pnetCDF.h"
MODULE pnetCDF_IO_module
!+--------------------------------------------------------------------+
!| This module provides subroutines to perform the parallel IO        |
!| with parallel netCDF.                                              |
!| Currently, this is used for writing the simulation data to a disk  |
!| file for postprocession purposes (e.g. visualization) and to write |
!| restart data to a single file in a portable way. It is thus easy   |
!| to restart a computation on a different machine and with a         |
!| different number of processors.                                    |
!+--------------------------------------------------------------------+
  USE defprecision_module
  USE defs_2D_3D_module
  IMPLICIT NONE
  SAVE
!****************************************************
!****  Global Data for IO                        ****
!****************************************************

  LOGICAL  :: write_pnetCDF_sim_dat ! write simulation data in
                                    ! parallel to a netCDF file?
  CHARACTER (len = 100) :: netCDF_simdat_file_name  ! file name for simulation data
  CHARACTER (len = 100) :: netCDF_in_simdat_file_name  ! file name used in simulation data file is 
                                                       ! used to restart a simulation 
  CHARACTER (len = 100) :: netCDF_in_dump_file_name,netCDF_out_dump_file_name ! file names for restart files



  INTEGER  :: ncid_simdat            ! ncid_simdat holds the netCDF 
                                     ! id for the file containing
                                     ! (single precision) simulation data                                      
  INTEGER  :: ncid_dump              ! ncid_dump holds netCDF id for 
                                     ! restart (dump) file
  
  ! Stuff for the simulation data file 
  INTEGER  :: x_dimid_simdat,y_dimid_simdat,z_dimid_simdat, &
            & time_dimid_simdat
  INTEGER  :: x_varid_simdat,y_varid_simdat,z_varid_simdat, &
            & time_varid_simdat,timestep_varid_simdat
  INTEGER  :: B_therm_varid_simdat,B_comp_varid_simdat
  INTEGER  :: D_visc_varid_simdat,D_therm_varid_simdat,D_comp_varid_simdat
  INTEGER  :: S_therm_varid_simdat,S_comp_varid_simdat
  INTEGER  :: T_part_varid_simdat,G_part_varid_simdat,D_part_varid_simdat,S_part_varid_simdat
  INTEGER  :: R_part_varid_simdat, Dv_part_varid_simdat
  INTEGER  :: Temp_varid_simdat
  INTEGER  :: Chem_varid_simdat
  INTEGER  :: Part_varid_simdat
  INTEGER  :: ux_varid_simdat,uy_varid_simdat,uz_varid_simdat
  INTEGER  :: upx_varid_simdat,upy_varid_simdat,upz_varid_simdat
  INTEGER  :: Gammax_varid_simdat,Gammay_varid_simdat,Gammaz_varid_simdat
  INTEGER  :: CFL_varid_simdat,dt_varid_simdat,dt_max_varid_simdat,dt_initial_varid_simdat

  ! Stuff for the restart files 
  INTEGER  :: ri_dimid_dump,l_dimid_dump,m_dimid_dump,n_dimid_dump,time_dimid_dump
  INTEGER  :: xy_dimid_dump
  INTEGER  :: Gammax_varid_dump,Gammay_varid_dump,Gammaz_varid_dump
  INTEGER  :: B_therm_varid_dump,B_comp_varid_dump
  INTEGER  :: D_visc_varid_dump,D_therm_varid_dump,D_comp_varid_dump
  INTEGER  :: S_therm_varid_dump,S_comp_varid_dump
  INTEGER  :: T_part_varid_dump,G_part_varid_dump,D_part_varid_dump,S_part_varid_dump
  INTEGER  :: R_part_varid_dump, Dv_part_varid_dump
  INTEGER  :: istep_varid_dump,time_varid_dump
  INTEGER  :: dt_varid_dump
  INTEGER  :: ri_varid_dump,kx_varid_dump,ky_varid_dump,kz_varid_dump
  INTEGER  :: Temp_varid_dump
  INTEGER  :: Chem_varid_dump
  INTEGER  :: Part_varid_dump
  INTEGER  :: ux_varid_dump,uy_varid_dump,uz_varid_dump
  INTEGER  :: upx_varid_dump,upy_varid_dump,upz_varid_dump

CONTAINS

#include  "./pnetCDF_IO/pn_open_dump.f90"
#include  "./pnetCDF_IO/pn_write_dump.f90"
#include  "./pnetCDF_IO/pn_read_size_and_pa_from_dump.f90"
#include  "./pnetCDF_IO/pn_read_state_from_dump.f90"
#include  "./pnetCDF_IO/pn_read_size_and_pa_from_simdat.f90"
#include  "./pnetCDF_IO/pn_read_state_from_simdat.f90"
#include  "./pnetCDF_IO/pn_close_dump.f90"

#include  "./pnetCDF_IO/pn_open_simdat_file.f90"
#include  "./pnetCDF_IO/pn_write_step_simdat_file.f90"
#include  "./pnetCDF_IO/pn_close_simdat_file.f90"

#include  "./pnetCDF_IO/pn_check.f90"
  
END MODULE pnetCDF_IO_module
