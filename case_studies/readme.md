# Case studies {#case_studies}

This folder contains *benchmark* case studies of `simul++` software. Systems were chosen such that their complexity increases from the simplest harmonic oscillator to the mixtures of two types of molecules. It should cover all systems which can be studied by the current version of `simul++`. 

Each *case study* should dwell in its folder with all input files needed for its launch. A bash script should be included which enables automatic recalculation of results in case of bug fix in the main program. The folder should also contain a brief description. If multiple methods are to be compared, scripts useful for figure generation should be conserved.

If possible, each *case study* should have its DL_POLY duplicate for comparison.

The list of case studies follows:

| CS No. |                     name                     |     atoms | brief description                                        | page          |
|-------:|:--------------------------------------------:|----------:|:---------------------------------------------------------|:--------------|
|      1 |  [harmonic oscillator](cs01/cs01-readme.md)  |         2 | the simplest case study comparable with previous results | @subpage cs01 |
|      2 |       [dump-bell](cs02/cs02-readme.md)       |         2 | the simplest possible test of SHAKE                      | @subpage cs02 |
|      3 |         [argon](cs03/cs03-readme.md)         |        64 | monoatomic LJ fluid                                      | @subpage cs03 |
|      4 |       [nitrogen](cs04/cs04-readme.md)        |       128 | small system with rigid bonds                            | @subpage cs04 |
|      5 |   [nitrogen – bigger](cs05/cs05-readme.md)   |      1024 | medium-size system with rigid bonds                      | @subpage cs05 |
|      6 |    [mixture Ar+N_2](cs06/cs06-readme.md)     | 256 + 256 | simple mixture...                                        | @subpage cs06 |
|      7 | [nitrogen – thermostat](cs07/cs07-readme.md) |       128 | small nitrogen system for thermostat tests               | @subpage cs07 |
|      8 |  [nitrogen – barostat](cs08/cs08-readme.md)  |       128 | small nitrogen system for barostat tests                 | @subpage cs08 |
