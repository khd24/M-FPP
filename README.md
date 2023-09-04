# M-FPP
Optimal Placement of Multiple Finite-size Rectangular Facilities in an Existing Layout

@article{date2023P1,
  title={Optimal Placement of M Finite-size Rectangular Facilities in an
Existing Layout: Part I: Theory},
  author={Date, K. and Nagi, R.},
  journal={{International {J}ournal on {P}roduction {R}esearch},
  volume={submitted},
  year={2023}
}

@article{date2023P2,
  title={Optimal Placement of M Finite-size Rectangular Facilities in an
Existing Layout: Part II: Solution Methods},
  author={Date, K. and Nagi, R.},
  journal={{International {J}ournal on {P}roduction {R}esearch},
  volume={submitted},
  year={2023}
}

Clone the repository on your local and run make command to build. Use the following macros while building.
EXTERNAL (to be able to provide command line arguments. If this is not defined, then the code will use problemset_4_3 by default).

Other macros (define while building the executable)
MAX_TIME (time in seconds for Explicit/Implicit enumeration)
MAX_GAP (e.g. 0.01 to specify the optimality gap for implicit enumeration)

**Runtime arguments:**

arg[1]: input file, arg[2]: log file, arg[3]: layout file, arg[4]: mode, arg[5]: heu bucketsize, arg[6]: repair mode, arg[7]: repair bucketsize

**Example:
**
FacilityPlacement.exe problemset_4_3 output_4_3 layout_4_3 3 1 1 1




