CMake_Minimum_Required(VERSION 2.8)

Project(Prototypes)

If (COMMAND cmake_policy)
    cmake_policy(SET CMP0003 NEW)
Endif (COMMAND cmake_policy)

# ITK :
Find_Package(ITK)
If (ITK_FOUND)
    Include (${ITK_USE_FILE})
Else (ITK_FOUND)
    Message (FATAL_ERROR "ITK not found. Please set ITK_DIR.")
Endif (ITK_FOUND)

# VTK :
Find_Package(VTK)
If (VTK_FOUND)
    If ( "${VTK_MAJOR_VERSION}" LESS 6 )
        Message ("VTK ${VTK_MAJOR_VERSION}.${VTK_MINOR_VERSION} is used.")
	      Include (${VTK_USE_FILE})
    Else ( "${VTK_MAJOR_VERSION}" LESS 6 )
        Find_Package(VTK REQUIRED)
        If (VTK_FOUND)
            Message ("VTK ${VTK_MAJOR_VERSION}.${VTK_MINOR_VERSION} is used.")
            Include (${VTK_USE_FILE})
        Else (VTK_FOUND)
            Message (FATAL_ERROR "Some components of VTKhave not been found.")
        Endif (VTK_FOUND)
    Endif ( "${VTK_MAJOR_VERSION}" LESS 6 )
Else (VTK_FOUND)
    Message (FATAL_ERROR "VTK not found. Please set VTK_DIR.")
Endif (VTK_FOUND)

# OpenMP
If(UNIX AND NOT APPLE)
  set(LINUX TRUE)
endif (UNIX AND NOT APPLE)

If (LINUX)
  message ("OS is Linux")
  Set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
Elseif (APPLE)
  message ("OS is MAC OS X")
  Set(CMAKE_CXX_COMPILER "/usr/local/bin/clang-omp++")
  Set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")
Else ()
  Message (FATAL_ERROR "Code works only on Linux or MAC OS X")
Endif ()

# Release version by default :
If (NOT DEFINED CMAKE_BUILD_TYPE:STRING)
    Set (CMAKE_BUILD_TYPE:STRING Release)
Endif (NOT DEFINED CMAKE_BUILD_TYPE:STRING)

# Include Directories
Include_Directories(${CMAKE_CURRENT_SOURCE_DIR})

# EIGEN
Include_Directories(../../Eigen)

#
# Executable :
#

Add_Executable(Gramiam
compute_gramiam.cpp
)

Add_Executable(MedoidsFinale
Medoids_Finale.cpp
)

Add_Executable(WriteTube
Write_Tube.cpp
)

Add_Executable(WriteTube2
Write_Tube_2.cpp
)

Add_Executable(Density
Density.cpp
)

#
# Target link :
#
Set ( VTK_DEPENDANCIES "" )
If ( "${VTK_MAJOR_VERSION}" LESS 6 )
    List (APPEND VTK_DEPENDANCIES
        vtkCommon
        vtkGraphics
        vtkIO
        vtkRendering
        vtkFiltering
    )
Else ( "${VTK_MAJOR_VERSION}" LESS 6 )
    List (APPEND VTK_DEPENDANCIES
        ${VTK_LIBRARIES}
    )
Endif ( "${VTK_MAJOR_VERSION}" LESS 6 )

Target_Link_Libraries(Gramiam
  ${ITK_LIBRARIES}
  ${VTK_DEPENDANCIES}
)

Target_Link_Libraries(MedoidsFinale
  ${ITK_LIBRARIES}
  ${VTK_DEPENDANCIES}
)

Target_Link_Libraries(WriteTube
  ${ITK_LIBRARIES}
  ${VTK_DEPENDANCIES}
)

Target_Link_Libraries(WriteTube2
  ${ITK_LIBRARIES}
  ${VTK_DEPENDANCIES}
)

Target_Link_Libraries(Density
  ${ITK_LIBRARIES}
  ${VTK_DEPENDANCIES}
)
