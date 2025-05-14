This is a plasma periodic along the x-direction bounded between two metal electrodes in the y-direction. 
One electrode (boundary object #3) is grounded.
The electrostatic potential of the other electrode (#1) is a prescribed function of time. 
Actually, it is a product of two functions.
One function gives high frequency oscillations with the shape defined in file init_bo_01_waveform.dat.
The other function controls variation of the amplitude of the oscillating potential with time. 
The shape of this function is defined in file init_bo_01_amplitude_profile.dat. The amplitude profile function is also periodic.
Its period is the time of the last point of the shape function defined in  init_bo_01_amplitude_profile.dat.

If the code does not find file init_bo_NN_amplitude_profile.dat, 
then the amplitude of the oscillating potential of boundary/inner object #NN (if requested) is constant in time.

The simulation runs on 16 MPI processes

To use the updated ion-neutral collision model use the init_neutral_Argon0_with_updated_ion_collision_model.dat file
Do
mv init_neutral_Argon0_with_updated_ion_collision_model.dat init_neutral_Argon0.dat
and pick the correct three parameters. The three parameters are:

1) beta_inf: Defines the impact parameter cutoff. Larger values lead to more accuracy, especially for large values of E/N (e.g. in the sheath or during ignition) while smaller values reduce the frequency of small angle collisions and, thus, increase performance. The value of 9 should work under most circumstances. The parameter is independent of species in might never need to be changed.
2) A: Parameter for the resonant charge-exchange cross-section. Can be obtained from the comparison between simulation and swarm experiments, found in the literature for many species, or might possibly be calculated from charge-exchange cross sections. Values from the original paper: Ar: 2.6, Ne: 2.6, He: 2.6, Kr: 2.5.
3) Polarizability: Can be found tabulated for basically any atom/molecule. Mind the units.

Model is from Nanbu and Kitatani: J. Phys. D: Appl. Phys. 28 (1995) 324-330
