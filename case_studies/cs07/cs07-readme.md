## Case study 07 – Nitrogen in PBC for thermostat tests {#cs07}

A simple system to test the NVT Nosé–Hoover ensemble in the new universal scheme.

Nitrogen parameters are taken from [this publication](https://www.researchgate.net/publication/257635210_Short-Time_Oxidation_Behavior_of_Low-Carbon_Low-Silicon_Steel_in_Air_at_850-1180_C_II_Linear_to_Parabolic_Transition_Determined_Using_Existing_Gas-Phase_Transport_and_Solid-Phase_Diffusion_Theories).

Expected results are obtained. Verlet-family methods give larger total energy fluctuations but no systematic drift. 
They also yields better constraints conservartion.
On the other hand, higher-order Gear methods give very small fluctuations of total energy but also a non-negligible drift.
Concerning energy conserving order, results support theoretical values.

(Edit 2022/06/13: Currently working as expected.)