#--------------------------------------------------------------
set(FEA_TESTS
    test_shellANCF
)


set(FEA_MKL_TESTS
    test_Antonio
    test_ANCFShell_Iso_Joint_Radu_Bullet
    test_ANCFShell_Iso_Joint_RaduLinkDir
)

message(STATUS "Additional test programs for FEA module...")

#--------------------------------------------------------------

foreach(PROGRAM ${FEA_TESTS})
    message(STATUS "...add ${PROGRAM}")

    add_executable(${PROGRAM}  "${PROGRAM}.cpp")
    source_group(""  FILES  "${PROGRAM}.cpp")
    set_target_properties(${PROGRAM} PROPERTIES 
        FOLDER demos
        COMPILE_FLAGS "${CH_CXX_FLAGS}"
        LINK_FLAGS "${CH_LINKERFLAG_EXE}")
    target_link_libraries(${PROGRAM}
        ChronoEngine
        ChronoEngine_fea)
    add_dependencies(${PROGRAM}
        ChronoEngine
        ChronoEngine_fea)
endforeach()

foreach(PROGRAM ${FEA_MKL_TESTS})
    message(STATUS "...add ${PROGRAM}")

    add_executable(${PROGRAM}  "${PROGRAM}.cpp")
    source_group(""  FILES  "${PROGRAM}.cpp")
    set_target_properties(${PROGRAM} PROPERTIES 
        FOLDER demos
        COMPILE_FLAGS "${CH_CXX_FLAGS}"
        LINK_FLAGS "${CH_LINKERFLAG_EXE}")
    target_link_libraries(${PROGRAM}
        ChronoEngine
        ChronoEngine_fea
        ChronoEngine_mkl)
    target_include_directories(${PROGRAM} PUBLIC
        ${CH_MKL_INCLUDES}
	)

endforeach()

