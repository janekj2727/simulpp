## Case study 06 – Nitrogen + Argon in PBC (simple mixture) {#cs06}

A simple mixture for tests of `simul++`.

Initial configuration was taken from `pokusy/pokus14/out2.configfinal` (result of `mergeconfig` tests).

Argon parameters from [this page](http://stp.clarku.edu/simulations/lj/index.html)

Nitrogen parameters are taken from [this publication](https://www.researchgate.net/publication/257635210_Short-Time_Oxidation_Behavior_of_Low-Carbon_Low-Silicon_Steel_in_Air_at_850-1180_C_II_Linear_to_Parabolic_Transition_Determined_Using_Existing_Gas-Phase_Transport_and_Solid-Phase_Diffusion_Theories).

Mixing rule (LJ): Lorentz–Berthelot

$\varepsilon_{ij} = \sqrt{\varepsilon_i \varepsilon_j}$ and $\sigma_{ij} = \frac {\sigma_i + \sigma_j}{2}$

Note: LJ parameters are probably not consistent and should be changed to some more rational values. Used values are a kind of *first found numbers*.

%Simulation results plotting methods are *echoed* at the end of the `cs05-run.sh` script. Measured quantities can be added directly to `.control` file.