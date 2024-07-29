# amrex-kitchen-cpp

This is the c++ version of the amrex-kitchen codes to perform analyses on AMReX plotfiles. 

This was mainly an experiment to see the differences in performance between `python` multiprocessing and c++ with MPI. 

Due to the very well optimized `numpy` package, the python version performs better with the same number of CPUs. 

Due to MPI overhead, the c++ version is only significantly faster with twice as much CPUs as the python version (~80), and AMReX tools are certainly better suited.
