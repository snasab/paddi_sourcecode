! +-----------------------------------------------------------------------------+
! | The folowing subroutine writes the fields (Temp,Chem,Part,u,curl u)(x=0,y=0,z)   |
! | to a disk file.                                                             |
! +-----------------------------------------------------------------------------+
SUBROUTINE wrtout_single_profile(u,Temp,Chem,Part,istep,dt,t)
  USE defprecision_module
  USE state_module
  USE message_passing_module, ONLY : myid
  USE parameter_module, ONLY :Nz
  USE mpi_transf_module, ONLY : mysx_phys,mysy_phys,mysz_phys,myez_phys,  &
  &                             mysx_spec,mysy_spec
  IMPLICIT NONE
  TYPE(velocity)   :: u
  TYPE(buoyancy)   :: Temp,Chem,Part
  INTEGER(kind=ki) :: istep
  REAL(kind=kr)    :: dt,t
  INTEGER(kind=ki) :: nprofile,ncount
  REAL(kind=kr),ALLOCATABLE :: profile(:,:)
  REAL(kind=kr),ALLOCATABLE :: local_z(:,:)
  INTEGER(kind=ki) :: k

  nprofile=6
#if defined(TEMPERATURE_FIELD)
  nprofile = nprofile + 1
#endif
#if defined(CHEMICAL_FIELD)
  nprofile = nprofile + 1
#endif
#if defined(PARTICLE_FIELD)
  nprofile = nprofile + 7
#endif

  
  allocate(profile(nprofile,0:Nz))
  allocate(local_z(nprofile,mysz_phys:myez_phys))

  ! copy the fields(x=0,y=0,z) to local_z if you have those fields in your
  ! local memory
  IF (mysx_phys.EQ.0 .AND. mysy_phys.EQ.0) THEN
     local_z(1,:)=u%phys(0,0,:,1)
     local_z(2,:)=u%phys(0,0,:,2)
     local_z(3,:)=u%phys(0,0,:,3)
     local_z(4,:)=u%curl(0,0,:,1)
     local_z(5,:)=u%curl(0,0,:,2)
     local_z(6,:)=u%curl(0,0,:,3)
     ncount=6
#if defined(TEMPERATURE_FIELD)
     local_z(ncount+1,:)=Temp%phys(0,0,:)
     ncount=ncount+1
#endif
#if defined(CHEMICAL_FIELD)
     local_z(ncount+1,:)=Chem%phys(0,0,:)
#endif
#if defined(PARTICLE_FIELD)
     local_z(ncount+1,:)=Part%phys(0,0,:)
     local_z(ncount+2,:)=up%phys(0,0,:,1)
     local_z(ncount+3,:)=up%phys(0,0,:,2)
     local_z(ncount+4,:)=up%phys(0,0,:,3)
     local_z(ncount+5,:)=up%curl(0,0,:,1)
     local_z(ncount+6,:)=up%curl(0,0,:,2)
     local_z(ncount+7,:)=up%curl(0,0,:,3)
#endif
  ENDIF

  ! collect vertical profile at x=y=0 at node 0
  CALL collect_profile(local_z,profile,nprofile)

  deallocate(profile,local_z)

end SUBROUTINE wrtout_single_profile


SUBROUTINE collect_profile(local_z,profile,nfields)
  USE defprecision_module
  USE message_passing_module, ONLY : myid,numtasks,nprocs2
  USE mpi_transf_module, ONLY : mysx_phys,myex_phys,mysy_phys,myey_phys,         &
  &                             mysz_phys,myez_phys,pencil_transpose_info_xy_zx 
  USE parameter_module, ONLY : Nz
#ifdef MPI_MODULE
  USE MPI
  IMPLICIT NONE
#else
  IMPLICIT NONE
  INCLUDE "mpif.h"
#endif
  INTEGER(kind=ki) :: nfields
  REAL(kind=kr)    :: local_z(nfields,mysz_phys:myez_phys),profile(nfields,0:Nz-1)
  INTEGER(kind=ki),ALLOCATABLE :: sz(0:nprocs2-1),ez(0:nprocs2-1) !sz(0:numtasks-1),ez(0:numtasks-1)
  INTEGER(kind=kiMPI) :: recvcounts(0:numtasks-1),rdispls(0:numtasks-1)
  INTEGER(kind=kiMPI) :: ierr,comm
  
#ifdef TWO_DIMENSIONAL
  comm = 2d_transp_info%communicator
#else
  comm = pencil_transpose_info_xy_zx%info_2nd_transpose%communicator
#endif

  IF (mysx_phys.EQ.0 .AND. mysy_phys.EQ.0) THEN
     CALL MPI_GATHER(mysz_phys,1,MPI_INTEGER,sz,1,MPI_INTEGER,0,           &
          &   comm,ierr)
     CALL MPI_GATHER(myez_phys,1,MPI_INTEGER,ez,1,MPI_INTEGER,0,           &
          &   comm,ierr)
     recvcounts(:) = (ez(:) - sz(:) + 1) * nfields
     rdispls(:) = sz(:) * nfields
     CALL MPI_GATHERV(local_z,(myez_phys-mysz_phys+1)*nfields,                &
          &    PM_MPI_FLOAT_TYPE,                                          &
          &    profile,recvcounts,rdispls,PM_MPI_FLOAT_TYPE,0,             &
          &    comm,ierr)
  ENDIF

END SUBROUTINE collect_profile
