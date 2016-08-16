CMake_Minimum_Required(VERSION 2.8)

Project(Prototypes)

If (COMMAND cmake_policy)
    cmake_policy(SET CMP0003 NEW)
Endif (COMMAND cmake_policy)

#
# Find packages (ITK & VTK) :
#
# ITK :
Find_Package(ITK)
If (ITK_FOUND)
    # Include (${USE_ITK_FILE})
    Include (${ITK_USE_FILE})
Else (ITK_FOUND)
    Message (FATAL_ERROR "ITK not found. Please set ITK_DIR.")
Endif (ITK_FOUND)

# VTK :
Find_Package(VTK)
If (VTK_FOUND)
    If ( "${VTK_MAJOR_VERSION}" LESS 6 )
        Message ("VTK ${VTK_MAJOR_VERSION}.${VTK_MINOR_VERSION} is used.")
        #Include (${USE_VTK_FILE})
	Include (${VTK_USE_FILE})
    Else ( "${VTK_MAJOR_VERSION}" LESS 6 )
        #Find_Package(VTK COMPONENTS vtkCommonCore vtkCommonDataModel vtkWrappingTools vtkIOLegacy vtkFiltersCore vtkFiltersGeneral vtkHybrid vtkWidgets NO_MODULE)
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

Include_Directories(
    ${CMAKE_CURRENT_SOURCE_DIR}
)

# OpenMP :
Set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopenmp")

# Release version by default :
If (NOT DEFINED CMAKE_BUILD_TYPE:STRING)
    Set (CMAKE_BUILD_TYPE:STRING Release)
Endif (NOT DEFINED CMAKE_BUILD_TYPE:STRING)


################## EIGEN ##############################
Include_Directories(
    ../../Eigen
)
################################################################################

#
# Executable :
#

Add_Executable(Gramiam
compute_gramiam.cpp
)

Add_Executable(Medoids
Medoids.cpp
)

Add_Executable(MedoidsNormalisation
Medoids_normalisation.cpp
)

Add_Executable(MedoidsFinale
Medoids_Finale.cpp
)

Add_Executable(WriteTube
Write_Tube.cxx
)

Add_Executable(WriteTubeColor
Write_Tube_color.cxx
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

Target_Link_Libraries(Medoids
  ${ITK_LIBRARIES}
  ${VTK_DEPENDANCIES}
)

Target_Link_Libraries(MedoidsNormalisation
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

Target_Link_Libraries(WriteTubeColor
  ${ITK_LIBRARIES}
  ${VTK_DEPENDANCIES}
)

Target_Link_Libraries(Density
  ${ITK_LIBRARIES}
  ${VTK_DEPENDANCIES}
)
