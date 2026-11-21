add_library(ust_cfd STATIC
    base/Tool.cpp
    field/Zone.cpp
    field/Variable.cpp
    material/Material.cpp
    solver/Solve.cpp
    solver/Solution.cpp
    physics/Equation.cpp
)

target_include_directories(ust_cfd PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
