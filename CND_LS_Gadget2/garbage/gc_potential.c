#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include "../allvars.h"
#include "../proto.h"

///FILE CREATED MAY 23, 2017
/// LAST MODIFICATIONS: March 9, 2016:
///                     Added modifications to the makefile, and flag to the acceleration computation
///                     Truncate the cluster acceleration: point source at > 8pc
//                      Overall potential = BAR+CLUSTER+DISK

