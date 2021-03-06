#+------------------------------------------------------------------------------+
# | Makefile for PADDI ([Pa]rallel [D]ouble [DI]ffusion) Code                    |
# +------------------------------------------------------------------------------+
#
# To compile, the following flags can be set:
# -------------------------------------------
#
#   -DMPI_MODULE         Set if the mpi implementation provides a f90 MPI Module.
#   -AB_BDF3             Use AB_BDF3 time stepping scheme.
#                        (a RK4_CN time stepping scheme used to exist in pervious versions.)
#                        (Currently, only AB_BDF3 is supported, so always use this flag!)
#   -DTEMPERATURE_FIELD  Include Temperature field
#   -DCHEMICAL_FIELD    Include Chemical field
#   -DPARTICLE_FIELD    Include Particle field
#   -DTWO_DIMENSIONAL   Build 2d code version
#                       (if not specified, the 3D version is generated)
#   -DSINGLE_PRECISION  Build code in single precision
#   -DDOUBLE_PRECISION  Build code in double precision
#
DEFS          = -DDOUBLE_PRECISION -DPARTICLE_FIELD -DMPI_MODULE -DAB_BDF3   

FC            = mpiifort 

F90           = $(FC) 

FFLAGS        = $(DEFS) 

F90FLAGS      = -I/cm/shared/apps/fftw/fftw-3.3.8/include -I/cm/shared/apps/netcdf-f/4.5.2/include -cpp -O0 $(DEFS) 

LD	      = $(FC)

LDFLAGS	      = -L/cm/shared/apps/fftw/fftw-3.3.8/lib -L/cm/shared/apps/netcdf-f/4.5.2/lib -L../../stuff_needed/lib 

LIBS	      = 

ADDLIBS       = -lfftw3 -lfftw3f -ljc -ljpeg -lpnetcdf -lirc

SRCS	      = defprecision_module.f90 defs_2D_3D_module.f90 diagnostics_module.f90 \
		diff_op_module.f90  integral_module.f90 IO_module.f90 \
		main_module.f90 message_passing_module.f90 minit_module.f90 \
		mpi_transf_module.f90 parameter_module.f90 \
		pnetCDF_IO_module.f90 prog.f90 state_module.f90 \
		testing_module.f90 transpose_pencil_module.f90 

PROGRAM	      = ../paper2_lux_conserve

OBJS          = defprecision_module.o defs_2D_3D_module.o diagnostics_module.o \
		diff_op_module.o  integral_module.o IO_module.o main_module.o \
		message_passing_module.o minit_module.o mpi_transf_module.o \
		parameter_module.o pnetCDF_IO_module.o prog.o state_module.o \
		testing_module.o transpose_pencil_module.o 

all:	 $(PROGRAM)

$(PROGRAM):     $(OBJS) $(LIBS)
		@echo "Linking $(PROGRAM) ..."
		@$(LD) $(LDFLAGS) $(OBJS) $(LIBS) $(ADDLIBS) -o $(PROGRAM)
		@echo "done"

state_module.o:	state_module.f90
		@$(F90) -I../../stuff_needed/include -cpp -O0 $(DEFS) -c state_module.f90


main_module.o:	main_module.f90
		@$(F90) -I../../stuff_needed/include -cpp -O0 $(DEFS) -c main_module.f90

clean:
	rm -f $(PROGRAM) $(OBJS) *.mod core

totalclean:
	rm -f $(PROG) $(OBJS) *.mod *.~ *.txt *.gp *.out *.eps

.SUFFIXES: $(SUFFIXES) .f90

.f90.o:
	$(F90) $(F90FLAGS) -c $<

.mod.o:
	$(F90) $(F90FLAGS) -c $*.f90

defs_2D_3D_module.o: defprecision_module.o
diagnostics_module.o: defprecision_module.o defs_2D_3D_module.o \
	integral_module.o mpi_transf_module.o parameter_module.o \
	state_module.o
diff_op_module.o: defprecision_module.o defs_2D_3D_module.o \
	mpi_transf_module.o parameter_module.o
integral_module.o: defprecision_module.o message_passing_module.o \
	mpi_transf_module.o parameter_module.o
IO_module.o: defprecision_module.o defs_2D_3D_module.o \
	message_passing_module.o parameter_module.o pnetCDF_IO_module.o \
	state_module.o testing_module.o
main_module.o: IO_module.o defprecision_module.o defs_2D_3D_module.o \
	diff_op_module.o message_passing_module.o minit_module.o mpi_transf_module.o \
	parameter_module.o pnetCDF_IO_module.o state_module.o 
message_passing_module.o: defprecision_module.o
minit_module.o: defprecision_module.o
mpi_transf_module.o: defprecision_module.o transpose_pencil_module.o
parameter_module.o: defprecision_module.o
pnetCDF_IO_module.o: defprecision_module.o defs_2D_3D_module.o
prog.o: IO_module.o defprecision_module.o main_module.o \
	message_passing_module.o parameter_module.o state_module.o 
state_module.o: defprecision_module.o defs_2D_3D_module.o diff_op_module.o \
	mpi_transf_module.o parameter_module.o
testing_module.o: defprecision_module.o diagnostics_module.o \
	mpi_transf_module.o parameter_module.o state_module.o
transpose_pencil_module.o: defprecision_module.o
diagnostics_module.o: defs_MPI.h diagnostics/rms_scalar.f90 \
	diagnostics/rms_spec_scalar.f90 diagnostics/rms_vector.f90 \
	diagnostics/rms_spec_vector.f90 diagnostics/compute_uTC_rms.f90 diagnostics/compute_upP_rms.f90 \
	diagnostics/compute_uTC_minmax.f90 diagnostics/compute_upP_minmax.f90 \
	diagnostics/compute_average_flux.f90 diagnostics/dissipation_buo.f90 
diff_op_module.o: diff_op/div.f90 diff_op/curl.f90 diff_op/d_by_dx.f90 \
	diff_op/d_by_dy.f90 diff_op/d_by_dz.f90
integral_module.o: defs_MPI.h integral/volume_integral.f90
IO_module.o: defs_MPI.h ../../stuff_needed/include/jcmagic.h \
	IO/open_files.f90 IO/write_output_files.f90 \
	IO/write_compressed_file.f90 IO/write_diagnostics_file.f90 \
	IO/write_horizontal_spectra.f90 IO/write_vertical_spectra.f90 \
	IO/write_z_profile.f90 IO/close_files.f90
main_module.o: defs_MPI.h main/RK2_CN/timestep_RK2_CN.f90 \
	main/RK2_CN/step1_RK2_CN.f90 main/RK2_CN/step2_RK2_CN.f90 \
	main/AB_BDF3/timestep_AB_BDF3.f90 \
	main/AB_BDF3/tmstp_buoyancy_AB_BDF3.f90 \
	main/AB_BDF3/tmstp_particle_AB_BDF3.f90 \
	main/AB_BDF3/tmstp_velocity_AB_BDF3.f90 \
	main/AB_BDF3/tmstp_particle_vel_AB_BDF3.f90 \
	main/AB_BDF3/comp_coeff_AB_BDF3.f90 main/AB_BDF3/CFL_AB_BDF3.f90 \
	main/AB_BDF3/adapt_dt_AB_BDF3.f90 main/read_parameter.f90 \
	main/init.f90 main/free_allocated_memory.f90 main/crhs_buoyancy.f90 main/crhs_particle.f90\
	main/crhs_velocity.f90 main/crhs_part_velocity.f90 main/compute_phys_space_vars.f90
message_passing_module.o: message_passing/start_mpi.f90 \
	message_passing/stop_mpi.f90
minit_module.o: defs_MPI.h minit/compute_particle_mass.f90 
mpi_transf_module.o: defs_MPI.h defs_FFTW.h mpi_transf/init_transforms.f90 \
	mpi_transf/init_transf_pencil.f90 mpi_transf/init_transf_2d.f90 \
	mpi_transf/init_x_transforms_pencil.f90 \
	mpi_transf/init_y_transforms_pencil.f90 \
	mpi_transf/init_z_transforms_pencil.f90 mpi_transf/FFT_r2c.f90 \
	mpi_transf/FFT_c2r.f90 mpi_transf/free_transforms.f90 \
	mpi_transf/f77_wisdom.f mpi_transf/get_FFTW_wisdom.f90 \
	mpi_transf/save_FFTW_wisdom.f90 mpi_transf/collect_phys_double.f90 \
	mpi_transf/collect_phys_single.f90 mpi_transf/collect_spec.f90 \
	mpi_transf/distribute_phys.f90 mpi_transf/distribute_spec.f90 \
	mpi_transf/collect_horizontal_sum_phys.f90
parameter_module.o: parameter/allocate_fft_storage_scheme.f90 \
	parameter/init_fft_storage_scheme.f90 \
	parameter/deallocate_fft_storage_scheme.f90
pnetCDF_IO_module.o: defs_pnetCDF.h ./pnetCDF_IO/pn_open_dump.f90 \
	./pnetCDF_IO/pn_write_dump.f90 \
	./pnetCDF_IO/pn_read_size_and_pa_from_dump.f90 \
	./pnetCDF_IO/pn_read_state_from_dump.f90 \
	./pnetCDF_IO/pn_read_size_and_pa_from_simdat.f90 \
	./pnetCDF_IO/pn_read_state_from_simdat.f90 \
	./pnetCDF_IO/pn_close_dump.f90 ./pnetCDF_IO/pn_open_simdat_file.f90 \
	./pnetCDF_IO/pn_write_step_simdat_file.f90 \
	./pnetCDF_IO/pn_close_simdat_file.f90 ./pnetCDF_IO/pn_check.f90
state_module.o: state/allocate_uTC.f90 state/deallocate_uTC.f90 \
	state/set_initial_condition.f90 state/init_u_phys.f90 state/init_up_phys.f90 \
	state/init_Temp_phys.f90 state/init_Chem_phys.f90 state/init_Part_phys.f90 \
	state/shift_time_pointers.f90 state/decomp_independent_random.f90 \
	state/make_solenoidal.f90 
testing_module.o: defs_MPI.h testing/error_Temp_adv_diff.f90 \
	testing/peak_div_u.f90
transpose_pencil_module.o: defs_MPI.h transpose_pencil/decomp_1d.f90 \
	transpose_pencil/init_transpose_1d_decomp.f90 \
	transpose_pencil/init_transpose_pencil_decomp.f90 \
	transpose_pencil/free_transpose_1d_decomp.f90 \
	transpose_pencil/free_transpose_pencil_decomp.f90 \
	transpose_pencil/pencil_decomp_inquiry_functions.f90 \
	transpose_pencil/create_pencil_mpi_communicators.f90 \
	transpose_pencil/transp1d_123_2dec_to_231_1dec.f90 \
	transpose_pencil/transp1d_231_1dec_to_123_2dec.f90 \
	transpose_pencil/transp1d_123_2dec_to_132_1dec.f90 \
	transpose_pencil/transp1d_132_1dec_to_123_2dec.f90 \
	transpose_pencil/cyclic_pencil_transposes.f90 \
	transpose_pencil/anti_cyclic_pencil_transposes.f90
