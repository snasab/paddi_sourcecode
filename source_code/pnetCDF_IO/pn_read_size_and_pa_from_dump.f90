! This subroutine reads some header information from the restart netCDF file
! and determines if everything is consistent with the user given information
SUBROUTINE pn_read_size_and_pa_from_dump
  USE defprecision_module
  USE parameter_module,   ONLY : Lmax,Mmax,Nmax,Gammax,Gammay,Gammaz,    &
  &                              B_therm,B_comp,D_visc,D_therm,D_comp,   &
  &                              S_therm,S_comp,T_part,G_part,D_part,    &
  &                              S_part,R_part,Dv_part
  USE message_passing_module, ONLY: start_mpi,stop_mpi,myid
#ifdef MPI_MODULE
  USE MPI
  IMPLICIT NONE
#else 
  IMPLICIT NONE
  INCLUDE "mpif.h"
#endif
  INTEGER (kind=ki) :: Lmax_dump,Mmax_dump,Nmax_dump
  REAL(kind=kr)     :: Gammax_dump,Gammay_dump,Gammaz_dump
  REAL(kind=kr)     :: B_therm_dump,B_comp_dump, &
                   &   D_visc_dump,D_therm_dump,D_comp_dump,S_therm_dump,S_comp_dump, &
                   &   T_part_dump,G_part_dump,D_part_dump,S_part_dump,R_part_dump,   &
                   &   Dv_part_dump
  INTEGER (kind=MPI_OFFSET_KIND) :: idummy
#include 'pnetcdf.inc'
!
  ! Open the restart file - read only
  CALL pn_check( nfmpi_open(MPI_COMM_WORLD,TRIM(netCDF_in_dump_file_name), &
       &  NF_NOWRITE,MPI_INFO_NULL,ncid_dump                  ) )
  ! get the dimensions
  CALL pn_check (nfmpi_inq_dimid(ncid_dump,'ri',ri_dimid_dump)      )
  CALL pn_check (nfmpi_inq_dimid(ncid_dump,'n',n_dimid_dump)        )
  CALL pn_check (nfmpi_inq_dimlen(ncid_dump,n_dimid_dump,idummy)    )
  Nmax_dump = idummy / 2 
  CALL pn_check (nfmpi_inq_dimid(ncid_dump,'l',l_dimid_dump)        )
  CALL pn_check (nfmpi_inq_dimlen(ncid_dump,l_dimid_dump,idummy)    )
  Lmax_dump = idummy - 1 
  CALL pn_check (nfmpi_inq_dimid(ncid_dump,'m',m_dimid_dump)        )
  CALL pn_check (nfmpi_inq_dimlen(ncid_dump,m_dimid_dump,idummy)    )
  Mmax_dump = idummy / 2 
  !
  IF ( Lmax.NE.Lmax_dump  .OR. Mmax.NE.Mmax_dump .OR. Nmax.NE.Nmax_dump )  THEN
     WRITE(*,'(a,i4,a)') "Dump file for process ",myid, &
          &                   "is not valid for this problem size."
     WRITE(*,'(a)') "Aborting..."
     CALL stop_mpi
     STOP
  ENDIF
  
  ! get the varids
  CALL pn_check (nfmpi_inq_varid(ncid_dump,'Gammax',Gammax_varid_dump)       )
#ifndef TWO_DIMENSIONAL
  CALL pn_check (nfmpi_inq_varid(ncid_dump,'Gammay',Gammay_varid_dump)       )
#endif
  CALL pn_check (nfmpi_inq_varid(ncid_dump,'Gammaz',Gammaz_varid_dump)       )

  CALL pn_check (nfmpi_inq_varid(ncid_dump,'B_therm',B_therm_varid_dump)     )
  CALL pn_check (nfmpi_inq_varid(ncid_dump,'B_comp', B_comp_varid_dump)      )
  CALL pn_check (nfmpi_inq_varid(ncid_dump,'D_visc',D_visc_varid_dump)       )
  CALL pn_check (nfmpi_inq_varid(ncid_dump,'D_therm',D_therm_varid_dump)     )
  CALL pn_check (nfmpi_inq_varid(ncid_dump,'D_comp' ,D_comp_varid_dump)      )
  CALL pn_check (nfmpi_inq_varid(ncid_dump,'S_therm' ,S_therm_varid_dump)    )
  CALL pn_check (nfmpi_inq_varid(ncid_dump,'S_comp' ,S_comp_varid_dump)      )
  CALL pn_check (nfmpi_inq_varid(ncid_dump,'T_part' ,T_part_varid_dump)      )
  CALL pn_check (nfmpi_inq_varid(ncid_dump,'G_part' ,G_part_varid_dump)      )
  CALL pn_check (nfmpi_inq_varid(ncid_dump,'Dv_part',Dv_part_varid_dump)      )
  CALL pn_check (nfmpi_inq_varid(ncid_dump,'D_part' ,D_part_varid_dump)      )
  CALL pn_check (nfmpi_inq_varid(ncid_dump,'S_part' ,S_part_varid_dump)      )
  CALL pn_check (nfmpi_inq_varid(ncid_dump,'R_part' ,R_part_varid_dump)      )
  CALL pn_check (nfmpi_inq_varid(ncid_dump,'istep',istep_varid_dump)         )
  CALL pn_check (nfmpi_inq_varid(ncid_dump,'time',time_varid_dump)           )
  CALL pn_check (nfmpi_inq_varid(ncid_dump,'dt',dt_varid_dump)               )
  CALL pn_check (nfmpi_inq_varid(ncid_dump,'ri',ri_varid_dump)               )
  CALL pn_check (nfmpi_inq_varid(ncid_dump,'kx',kx_varid_dump)               )
#ifndef TWO_DIMENSIONAL
  CALL pn_check (nfmpi_inq_varid(ncid_dump,'ky',ky_varid_dump)               )
#endif
  CALL pn_check (nfmpi_inq_varid(ncid_dump,'kz',kz_varid_dump)               )
#ifdef TEMPERATURE_FIELD
  CALL pn_check (nfmpi_inq_varid(ncid_dump,'Temp',Temp_varid_dump)  )
#endif
#ifdef CHEMICAL_FIELD
  CALL pn_check (nfmpi_inq_varid(ncid_dump,'Chem',Chem_varid_dump)  )
#endif
  CALL pn_check (nfmpi_inq_varid(ncid_dump,'ux',ux_varid_dump)      )
#ifndef TWO_DIMENSIONAL
  CALL pn_check (nfmpi_inq_varid(ncid_dump,'uy',uy_varid_dump)      )
#endif
  CALL pn_check (nfmpi_inq_varid(ncid_dump,'uz',uz_varid_dump)      )
  
  
#ifdef PARTICLE_FIELD
  CALL pn_check (nfmpi_inq_varid(ncid_dump,'Part',Part_varid_dump)  )
  CALL pn_check (nfmpi_inq_varid(ncid_dump,'upx',upx_varid_dump)      )
#ifndef TWO_DIMENSIONAL
  CALL pn_check (nfmpi_inq_varid(ncid_dump,'upy',upy_varid_dump)      )
#endif
  CALL pn_check (nfmpi_inq_varid(ncid_dump,'upz',upz_varid_dump)      )
#endif

  
  ! read Parameter values and check if they have changed
  ! write warning to standart output if they have 
  CALL pn_check( PM_NFMPI_GET_VAR_FLOAT_ALL(ncid_dump,Gammax_varid_dump,Gammax_dump) )
#ifndef TWO_DIMENSIONAL
  CALL pn_check( PM_NFMPI_GET_VAR_FLOAT_ALL(ncid_dump,Gammay_varid_dump,Gammay_dump) )
#endif
  CALL pn_check( PM_NFMPI_GET_VAR_FLOAT_ALL(ncid_dump,Gammaz_varid_dump,Gammaz_dump) )
  CALL pn_check( PM_NFMPI_GET_VAR_FLOAT_ALL(ncid_dump,B_therm_varid_dump,B_therm_dump) )
  CALL pn_check( PM_NFMPI_GET_VAR_FLOAT_ALL(ncid_dump,B_comp_varid_dump,B_comp_dump) )
  CALL pn_check( PM_NFMPI_GET_VAR_FLOAT_ALL(ncid_dump,D_visc_varid_dump,D_visc_dump) )
  CALL pn_check( PM_NFMPI_GET_VAR_FLOAT_ALL(ncid_dump,D_therm_varid_dump,D_therm_dump) )
  CALL pn_check( PM_NFMPI_GET_VAR_FLOAT_ALL(ncid_dump,D_comp_varid_dump,D_comp_dump) )
  CALL pn_check( PM_NFMPI_GET_VAR_FLOAT_ALL(ncid_dump,S_therm_varid_dump,S_therm_dump) )
  CALL pn_check( PM_NFMPI_GET_VAR_FLOAT_ALL(ncid_dump,S_comp_varid_dump,S_comp_dump) )
  CALL pn_check( PM_NFMPI_GET_VAR_FLOAT_ALL(ncid_dump,T_part_varid_dump,T_part_dump) )
  CALL pn_check( PM_NFMPI_GET_VAR_FLOAT_ALL(ncid_dump,G_part_varid_dump,G_part_dump) )
  CALL pn_check( PM_NFMPI_GET_VAR_FLOAT_ALL(ncid_dump,Dv_part_varid_dump,Dv_part_dump) )
  CALL pn_check( PM_NFMPI_GET_VAR_FLOAT_ALL(ncid_dump,D_part_varid_dump,D_part_dump) )
  CALL pn_check( PM_NFMPI_GET_VAR_FLOAT_ALL(ncid_dump,S_part_varid_dump,S_part_dump) )
  CALL pn_check( PM_NFMPI_GET_VAR_FLOAT_ALL(ncid_dump,R_part_varid_dump,R_part_dump) )
  IF (Gammax_dump.NE.Gammax) WRITE(*,'(a,i4,a)') "Warning: Gammax on id ",myid, &
       & " is different from value in restart file."
#ifndef TWO_DIMENSIONAL
  IF (Gammay_dump.NE.Gammay) WRITE(*,'(a,i4,a)') "Warning: Gammay on id ",myid, &
       & " is different from value in restart file."
#endif
  IF (Gammaz_dump.NE.Gammaz) WRITE(*,'(a,i4,a)') "Warning: Gammaz on id ",myid, &
       & " is different from value in restart file."
  IF (B_therm_dump.NE.B_therm) WRITE(*,'(a,i4,a)') "Warning: B_therm on id ",myid, &
       & " is different from value in restart file."
  IF (B_comp_dump.NE.B_comp) WRITE(*,'(a,i4,a)') "Warning: B_comp on id ",myid, &
       & " is different from value in restart file."
  IF (D_visc_dump.NE.D_visc) WRITE(*,'(a,i4,a)') "Warning: D_visc on id ",myid, &
       & " is different from value in restart file."
  IF (D_therm_dump.NE.D_therm) WRITE(*,'(a,i4,a)') "Warning: D_therm on id ",myid, &
       & " is different from value in restart file."
  IF (D_comp_dump.NE.D_comp) WRITE(*,'(a,i4,a)') "Warning: D_comp on id ",myid, &
       & " is different from value in restart file."
  IF (S_therm_dump.NE.S_therm) WRITE(*,'(a,i4,a)') "Warning: S_therm on id ",myid, &
       & " is different from value in restart file."
  IF (S_comp_dump.NE.S_comp) WRITE(*,'(a,i4,a)') "Warning: S_comp on id ",myid, &
       & " is different from value in restart file."
  IF (T_part_dump.NE.T_part) WRITE(*,'(a,i4,a)') "Warning: T_part on id ",myid, &
  	   & " is different from value in restart file."
  IF (G_part_dump.NE.G_part) WRITE(*,'(a,i4,a)') "Warning: G_part on id ",myid, &
  	   & " is different from value in restart file."
  IF (Dv_part_dump.NE.Dv_part) WRITE(*,'(a,i4,a)') "Warning: Dv_part on id ",myid, &
       & " is different from value in restart file."
  IF (D_part_dump.NE.D_part) WRITE(*,'(a,i4,a)') "Warning: D_part on id ",myid, &
       & " is different from value in restart file."
  IF (S_part_dump.NE.S_part) WRITE(*,'(a,i4,a)') "Warning: S_part on id ",myid, &
       & " is different from value in restart file."
  IF (R_part_dump.NE.R_part) WRITE(*,'(a,i4,a)') "Warning: R_part on id ",myid, &
       & " is different from value in restart file."

END SUBROUTINE pn_read_size_and_pa_from_dump
