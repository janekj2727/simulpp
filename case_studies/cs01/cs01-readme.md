## Case study 01 – Harmonic oscillator {#cs01}

The simplest system with known exact solution. Useful for the tests of correctness of the integration. Different integrators in NVE are known to give slightly different velocities. Energy conservation was studied in detail in my diploma thesis.

Mass of *atoms* X is 1.0, harmonic bond is used to mimic the spring of the oscillator. In the initial state, atoms are 3 [AA] apart with no velocity. The equilibrium bond length is 2 [AA], giving the total energy 0.5k [K]. The size of the atoms can be modified by setting the sigma parameter of lj potential in the `.field` file.   

### Expected results

For timestep h=0.1[ps] energy is conserved quite well (except for `k3m3` method, where the conservation is generally poor). The value of E_tot should be 0.5 (0.5 k_B [K]). `verlet` (with VERLET=1, check making options) should give identical results with `k3m2e` with energy oscillating below the 0.5 level. Meanwhile, with `k4m2e` (another verlet equivalent method) the total energy is slightly more then 0.5 and the oscillations are smaller. Higher order methods give much smaller fluctuations of E_tot, but do not exactly conserve energy. General trends are observable: `k4m4` – energy decreases a lot, `k5m5` – energy increases, `k5m4e` – energy decreases a little.

For h=0.2[ps] trends are the same, but the higher order methods are no longer advantageous because of the energy drift...

The list of integrators, as well as other options, can be easily modified within the `cs02-run.sh` script.