append_source_files(SOURCES)

add_library(PhaseFieldAcid ${SOURCES})
if(BUILD_SHARED_LIBS)
    install(TARGETS PhaseFieldAcid LIBRARY DESTINATION ${CMAKE_INSTALL_LIBDIR})
endif()
target_link_libraries(PhaseFieldAcid PUBLIC ProcessLib PRIVATE ParameterLib)

include(Tests.cmake)
