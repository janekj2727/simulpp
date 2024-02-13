## Case study 02 – dumb-bell {#cs02}

The simplest system to test the integration of constrained bond (SHAKE). No real forces are acting on the system, only fictitious forces are computed that keep the constrained bond length constant.

In contrast with a quick view on preliminary results obtained in `pokus07`, it seems that the energy trend here is the same as in the case of a NVE simulation. Indeed, that makes sense. The fictitious force is eventually a *force* and in spite of absence of real forces, fictitious forces are not integrated reversibly...

The script produces two graphs, each with log axes. The first shows the dependence of energy slope on time step, whereas the second shows the error of constrained bond in the last integration step (being the same as the whole simulation except for the very beginning) as a function of time step. Overall, the estimate of method error order can be done, which is roughly equal to the expected value...

Sometimes, there is a *segmentation fault* during the execution, but the bug is not reproducible and `valgrind` shows no error... weird...