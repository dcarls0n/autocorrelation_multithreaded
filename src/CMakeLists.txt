#Require Ver 3.6 for Intel compiler support
cmake_minimum_required(VERSION 3.6)
project(autocorrelate_threads LANGUAGES CXX)

add_executable(autocorrelate_threads autocorrelate_threads.cpp Barrier.hh)
find_package(Threads REQUIRED)

target_compile_features(autocorrelate_threads INTERFACE cxx_std_14)

set_target_properties(autocorrelate_threads PROPERTIES
	CXX_STANDARD 14
	CXX_EXTENSIONS FALSE
	CXX_STANDARD_REQUIRED TRUE
	LINK_FLAGS "${CMAKE_THREAD_LIBS_INIT}"
)

set(CMAKE_CXX_STANDARD 14)