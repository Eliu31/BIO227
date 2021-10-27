# BIO227

This is a draft script that attempts to incorporate empty patches into the model shown in Velland Box 6. Extra parameters include Np (initial abundance of species 1) and initial.patch (number of patches initially occupied). All occupied patches are assumed to be completely filled, and non-occupied patches are assumed to be completely empty. Accordingly, the following changes have also been made to the main simulation code:

1. Pr.1 is now calculated by dividing the abundance of species 1 by the abundance of LIVING organisms. This prevents empty spaces from messing up the propagation ratio.
2. Added a variable to differentiate between instances of local reproduction and instances of dispersal.
3. Using the variable above, local reproduction was made population-dependent in incompletely occupied patches. Local reproduction probability increases (from 0 to its base chance of 1) as the number of individuals approaches carrying capacity. This prevents individuals from spontanously appearing in empty patches, or patches with 1 individual from reproducing at the same rate as patches with 100 individuals.

Unfortunately, I haven't spent any time really checking the model or cleaning up the code. There are probably bugs/instances of odd behavior that still need to be fixed.


