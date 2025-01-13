#pragma once

#include"ode_solver.h"
#include "solution.h"

matrix ff0T(matrix, matrix = NAN, matrix = NAN);

matrix ff0R(matrix, matrix = NAN, matrix = NAN);

matrix df0(double, matrix, matrix = NAN, matrix = NAN);

matrix ff1T(matrix, matrix = NAN, matrix = NAN);

matrix ff1R(matrix x, matrix ud1, matrix ud2);

matrix df1R(double t, matrix Y, matrix ud1, matrix ud2);

matrix solve_simulation(matrix Da, matrix ud1, matrix ud2);

matrix dV_derivative(matrix V, matrix D, matrix P);

matrix dT_derivative(matrix T, matrix V_Vin, matrix Tin);

matrix derivatives(double t, matrix Y, matrix UD1, matrix UD2);

matrix ff2T(matrix, matrix = NAN, matrix = NAN);

matrix ff2R(matrix x, matrix ud1, matrix ud2);

matrix df2R(double t, matrix Y, matrix ud1, matrix ud2);

matrix ff3T(matrix x, matrix ud1, matrix ud2);

matrix ff3R(matrix x, matrix ud1, matrix ud2);

matrix df3R(double t, matrix Y, matrix ud1, matrix ud2);

matrix gradient(matrix x, matrix ud1, matrix ud2);

matrix hessian(matrix x, matrix ud1, matrix ud2);

matrix ff4T(matrix x, matrix ud1, matrix ud2);

matrix ff4R(matrix x, matrix ud1, matrix ud2);

matrix gf(matrix x, matrix ud1, matrix ud2);

void changeSign(std::string &path, char sign1, char sign2);

matrix fT5(matrix x, matrix ud1, matrix ud2);
matrix fR5(matrix x, matrix ud1, matrix ud2);

matrix f5(matrix x, matrix ud1, matrix ud2);
matrix f5_1(double a, matrix x, matrix ud1, matrix ud2);
matrix f5_2(double a, matrix x, matrix ud1, matrix ud2);
