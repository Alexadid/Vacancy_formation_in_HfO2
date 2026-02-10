cmake_minimum_required(VERSION 3.16)
project(HfO2VacancyMC)

find_package(Geant4 REQUIRED ui_all vis_all)

include(${Geant4_USE_FILE})
include_directories(${PROJECT_SOURCE_DIR}/include)

file(GLOB sources ${PROJECT_SOURCE_DIR}/src/*.cc)
file(GLOB headers ${PROJECT_SOURCE_DIR}/include/*.hh)

add_executable(HfO2VacancyMC main.cc ${sources} ${headers})
target_link_libraries(HfO2VacancyMC ${Geant4_LIBRARIES})

# Copy macros (optional)
file(MAKE_DIRECTORY ${PROJECT_BINARY_DIR}/macros)
file(COPY ${PROJECT_SOURCE_DIR}/macros DESTINATION ${PROJECT_BINARY_DIR})