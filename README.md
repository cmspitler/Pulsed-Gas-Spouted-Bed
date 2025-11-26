The following files must be included to run the simulations:

particle_input.dat

SBP.mfx

usr1_des.f

particle_input.dat specifies the positions of all particles at the start of the simulation

SBP.mfx contains the settings for the simulations

usr1_des.f specifies a user-defined function for oscillating the gas flow and outputting the particle drag forces.


The following files are modified MFIX source files which are required for outputting the particle drag forces:

allocate_arrays.f

calc_drag_des.f

cfnewvalues.f

des_allocate_mod.f

des_functions.f

des_init_arrays.f

discretelement_mod.f

mpi_funs_des_mod.f

mpi_init_des_mod.f

mpi_pack_des_mod.f

mpi_unpack_des_mod.f
