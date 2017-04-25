include_directories(..)

if (NOT RAJA_ENABLE_NESTED)
  add_definitions(-DRAJA_ENABLE_NESTED)
endif ()

set(KRIPKE_SRC 
    Directions.cpp  
    Grid.cpp  
    Input_Variables.cpp
    Kernel.cpp
    Layout.cpp
    Subdomain.cpp  
    Sweep_Solver.cpp
    ParallelComm.cpp  
    Timing.cpp
      
    Kernel/Kernel_3d_GDZ.cpp
    Kernel/Kernel_3d_DGZ.cpp
    Kernel/Kernel_3d_ZDG.cpp
    Kernel/Kernel_3d_DZG.cpp
    Kernel/Kernel_3d_ZGD.cpp
    Kernel/Kernel_3d_GZD.cpp
    
    ParallelComm/BlockJacobiComm.cpp
    ParallelComm/SweepComm.cpp 
    
    Test/TestKernels.cpp)

raja_add_library(
  NAME lib_kripke_raja 
  SOURCES ${KRIPKE_SRC})  
add_dependencies(lib_kripke_raja RAJA)
                                                                                                          
list(APPEND KRIPKE_LIBS lib_kripke_raja) 

# Propagate to parent directory
set(KRIPKE_LIBS ${KRIPKE_LIBS} PARENT_SCOPE) 
 
