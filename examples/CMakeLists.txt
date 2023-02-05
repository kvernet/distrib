cmake_minimum_required(VERSION 3.25.1)
project(Examples)

function(add_example_executable name)
	add_executable(${name} ${name}.cpp)
	target_link_libraries(${name} Distrib)
	install(TARGETS ${name} DESTINATION ${CMAKE_INSTALL_PREFIX}/bin)
endfunction()

# add normalDist exe
add_example_executable(normalDist)
add_example_executable(uniformDist)
add_example_executable(solidAngleDist)
add_example_executable(logarithmDist)
add_example_executable(poissonDist)
add_example_executable(uniformBinCenterDist)
