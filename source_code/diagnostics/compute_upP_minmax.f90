!+--------------------------------------------------------------------+
!| The following function computes the minumum and maximum of         |
!|         P,up_x,up_y,up_z,VORTp_x,VORTp_y,VORTp_z                       |
!| as well as the maxima of                                           |
!|                | up | and of | VORT |.                              |            
!| Caution: It is assumed that this routine is not called at every    |
!|          time step so that it does not need to be optimized.       |
!|          In fact it is likely to be highly inefficient...          |
!+--------------------------------------------------------------------+
SUBROUTINE compute_upP_minmax(Part_min,Part_max,up_min,up_max,       &
		   &                  VORTp_min,VORTp_max,                   &
		   &                  up_max_abs,VORTp_max_abs,              &
		   &                  up,Part                                )                 							 
		   
  USE defprecision_module
  USE state_module,ONLY : velocity,buoyancy
#ifdef MPI_MODULE
  USE MPI
  IMPLICIT NONE
#else
  IMPLICIT NONE
  INCLUDE "mpif.h"
#endif
  TYPE(velocity) :: up
  TYPE(buoyancy) :: Part
  REAL(kind=kr) :: Part_min,Part_max,                             &
	   &           up_min(3),up_max(3),VORTp_min(3),VORTp_max(3), &
	   &           up_max_abs,VORTp_max_abs
  REAL(kind=kr) :: min_local,max_local,nmin_local(3),nmax_local(3)
  INTEGER(kind=ki) :: ierr



  ! Particle fluctuation
#ifdef PARTICLE_FIELD
  min_local = MINVAL(Part%phys)
  max_local = MAXVAL(Part%phys) 
  CALL MPI_REDUCE(min_local,Part_min,1,PM_MPI_FLOAT_TYPE,MPI_MIN,0,MPI_COMM_WORLD,ierr)
  CALL MPI_REDUCE(max_local,Part_max,1,PM_MPI_FLOAT_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ierr)
#else 
  Part_min = 0._kr
  Part_max = 0._kr
#endif

!-----------------------------------------------
!VELOCITY,VORTICITY,MAX/MIN,ABS (particle below)
#ifdef PARTICLE_FIELD
  ! particle velocity components
  nmin_local(1) = MINVAL(up%phys(:,:,:,vec_x))
  nmax_local(1) = MAXVAL(up%phys(:,:,:,vec_x))
#ifndef TWO_DIMENSIONAL
  nmin_local(2) = MINVAL(up%phys(:,:,:,vec_y))
  nmax_local(2) = MAXVAL(up%phys(:,:,:,vec_y))
#else
  nmin_local(2) = 0._kr 
  nmax_local(2) = 0._kr
#endif
  nmin_local(3) = MINVAL(up%phys(:,:,:,vec_z))
  nmax_local(3) = MAXVAL(up%phys(:,:,:,vec_z))

  CALL MPI_REDUCE(nmin_local,up_min,3,PM_MPI_FLOAT_TYPE,MPI_MIN,0,MPI_COMM_WORLD,ierr)
  CALL MPI_REDUCE(nmax_local,up_max,3,PM_MPI_FLOAT_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ierr)

  ! vorticity components
  nmin_local(2) = MINVAL(up%curl(:,:,:,curl_y))
  nmax_local(2) = MAXVAL(up%curl(:,:,:,curl_y))
#ifndef TWO_DIMENSIONAL
  nmin_local(1) = MINVAL(up%curl(:,:,:,curl_x))
  nmax_local(1) = MAXVAL(up%curl(:,:,:,curl_x))
  nmin_local(3) = MINVAL(up%curl(:,:,:,curl_z))
  nmax_local(3) = MAXVAL(up%curl(:,:,:,curl_z))
#else 
  nmin_local(1) = 0._kr
  nmax_local(1) = 0._kr
  nmin_local(3) = 0._kr
  nmax_local(3) = 0._kr
#endif
  CALL MPI_REDUCE(nmin_local,VORTp_min,3,PM_MPI_FLOAT_TYPE,MPI_MIN,0,MPI_COMM_WORLD,ierr)
  CALL MPI_REDUCE(nmax_local,VORTp_max,3,PM_MPI_FLOAT_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ierr)

  ! maximum absolute value of velocity
#ifdef TWO_DIMENSIONAL
  max_local = MAXVAL( SQRT( up%phys(:,:,:,vec_x)**2 + up%phys(:,:,:,vec_z)**2 ) )
#else
  max_local = MAXVAL( SQRT( up%phys(:,:,:,vec_x)**2 &
              &           + up%phys(:,:,:,vec_y)**2 &
              &           + up%phys(:,:,:,vec_z)**2 ) )
#endif
  CALL MPI_REDUCE(max_local,up_max_abs,1,PM_MPI_FLOAT_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ierr)

  ! maximum absolute value of vorticity
#ifdef TWO_DIMENSIONAL
  max_local = MAXVAL( SQRT( up%curl(:,:,:,curl_y)**2 ) )
#else
  max_local = MAXVAL( SQRT( up%curl(:,:,:,curl_x)**2 &
              &           + up%curl(:,:,:,curl_y)**2 &
              &           + up%curl(:,:,:,curl_z)**2 ) )
#endif
  CALL MPI_REDUCE(max_local,VORTp_max_abs,1,PM_MPI_FLOAT_TYPE,MPI_MAX,0,MPI_COMM_WORLD,ierr)
  
#endif
  
  


END SUBROUTINE compute_upP_minmax
