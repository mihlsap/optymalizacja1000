cmake_minimum_required(VERSION 3.28)
project(optymalizacja9000_main)

set(CMAKE_CXX_STANDARD 20)

include_directories(.)

add_executable(optymalizacja9000_main
        main.cpp
        matrix.cpp
        matrix.h
        ode_solver.cpp
        ode_solver.h
        opt_alg.cpp
        opt_alg.h
        solution.cpp
        solution.h
        user_funs.cpp
        user_funs.h)
