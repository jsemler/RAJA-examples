if (RAJA_ENABLE_CALIPER)
  find_package(CALIPER)
  if(CALIPER_FOUND)
    message(STATUS "CALIPER")
    include_directories(${caliper_INCLUDE_DIR})
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DRAJA_USE_CALIPER")
  endif()
endif()

#Used for timing
find_library(RT_LIBRARIES rt)
if (RT_LIBRARIES STREQUAL "RT_LIBRARIES-NOTFOUND")
  message(WARNING "librt not found, some test applications might not link")
  set(RT_LIBRARIES "" CACHE STRING "timing libraries" FORCE)
endif ()
if (CALIPER_FOUND)
    set(RT_LIBRARIES "${RT_LIBRARIES} ${caliper_LIB_DIR}/libcaliper.so" CACHE STRING "testing libraries" FORCE)
endif()

include_directories(include)

add_subdirectory(LULESH-v1.0)
add_subdirectory(LULESH-v2.0)
add_subdirectory(Kripke-v1.1)
