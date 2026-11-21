set(CMAKE_RUNTIME_OUTPUT_DIRECTORY $<0:>${CMAKE_BINARY_DIR}/applications)

file(COPY "${CMAKE_CURRENT_SOURCE_DIR}/mesh"
     DESTINATION "${CMAKE_CURRENT_BINARY_DIR}")

add_executable(hypersonic_intake hypersonic_intake.cpp)

add_executable(volumetric_reaction volumetric_reaction.cpp)

add_executable(heat_conduction heat_conduction.cpp)
target_link_libraries(heat_conduction PRIVATE ust_cfd)

add_executable(flow_energy flow_energy.cpp)
target_link_libraries(flow_energy PRIVATE ust_cfd)

add_executable(incompress_flow incompress_flow.cpp)
target_link_libraries(incompress_flow PRIVATE ust_cfd)
