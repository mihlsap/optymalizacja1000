#include"user_funs.h"
#include <functional>

matrix ff0T(matrix x, matrix ud1, matrix ud2) {
    matrix y(1, 1);
    y(0) = pow(x(0), 2) + pow(x(1), 2) - cos(2.5 * M_PI * x(0)) - cos(2.5 * M_PI * x(1)) + 2;
    return y;
}

matrix ff0R(matrix x, matrix ud1, matrix ud2) {
    matrix y;
    matrix Y0 = matrix(2, 1), MT = matrix(2, new double[2]{m2d(x), 0.5});
    matrix *Y = solve_ode(df0, 0, 0.1, 10, Y0, ud1, MT);
    int n = get_len(Y[0]);
    double teta_max = Y[1](0, 0);
    for (int i = 1; i < n; ++i)
        if (teta_max < Y[1](i, 0))
            teta_max = Y[1](i, 0);
    y = abs(teta_max - m2d(ud1));
    Y[0].~matrix();
    Y[1].~matrix();
    return y;
}

matrix df0(double t, matrix Y, matrix ud1, matrix ud2) {
    matrix dY(2, 1);
    double m = 1, l = 0.5, b = 0.5, g = 9.81;
    double I = m * pow(l, 2);
    dY(0) = Y(1);
    dY(1) = ((t <= ud2(1)) * ud2(0) - m * g * l * sin(Y(0)) - b * Y(1)) / I;
    return dY;
}

matrix ff1T(matrix x, matrix ud1, matrix ud2) {
    return {-cos(0.1 * m2d(x)) * exp(-pow((0.1 * m2d(x) - 2 * M_PI), 2)) + 0.002 * pow(0.1 * m2d(x), 2)};
}

matrix ff1R(matrix x, matrix ud1, matrix ud2) {
    matrix Y0(3, 1);
    Y0(0) = 5.0;
    Y0(1) = 1.0;
    Y0(2) = 20.0;
    double t0 = 0.0;
    double dt = 1.0;
    int timesteps = 2000;

    matrix *Y = solve_ode(df1R, t0, dt, timesteps, Y0, x(0), x);

    double Tmax = Y[1](0, 2);
    for (int i = 1; i <= timesteps; ++i) {
        if (Y[1](i, 2) > Tmax) {
            Tmax = Y[1](i, 2);
        }
    }
    delete[] Y;

    return {fabs(Tmax - 50.0)};
}

matrix df1R(double t, matrix Y, matrix ud1, matrix ud2) {
    matrix derivatives(3, 1);
    double a = 0.98, b = 0.63, g = 9.81, PA = 0.5, PB = 1, DB = 0.00365665;
    double F_in = 0.01, T_in = 20, TA = 90.0;
    double DA = m2d(ud1(0));

    double FA_out, FB_out;

    if (Y(0) > 0)
        FA_out = a * b * DA * sqrt(2 * g * Y(0) / PA);
    else
        FA_out = 0;

    if (Y(1) > 0)
        FB_out = a * b * DB * sqrt(2 * g * Y(1) / PB);
    else
        FB_out = 0;

    derivatives(0) = -FA_out;
    derivatives(1) = FA_out + F_in - FB_out;
    derivatives(2) = FA_out / Y(1) * (TA - Y(2)) + F_in / Y(1) * (T_in - Y(2));

    return derivatives;
}

matrix ff2T(matrix x, matrix ud1, matrix ud2) {
    matrix y(1, 1);
    y(0) = pow(x(0), 2) + pow(x(1), 2) - cos(2.5 * M_PI * x(0)) - cos(2.5 * M_PI * x(1)) + 2;
    return y;
}

matrix df2R(double t, matrix Y, matrix ud1, matrix ud2) {
    double alpha = Y(0);
    double omega = Y(1);

    double k1 = ud2(0);
    double k2 = ud2(1);

    double m_arm = 1.0;
    double m_weight = 5.0;
    double l = 1.0;
    double g = 9.81;
    double b = 0.5;

    double I = (m_arm + m_weight) * pow(l, 2);

    double M = k1 * (M_PI - alpha) + k2 * (0.0 - omega);

    matrix dY(2, 1);
    dY(0) = omega;
    dY(1) = (M - m_arm * g * l * sin(alpha) - b * omega) / I;

    return dY;
}


matrix ff2R(matrix x, matrix ud1, matrix ud2) {
    double k1 = x(0);
    double k2 = x(1);

    double m = 1.0, l = 1.0, g = 9.81, b = 0.5;
    double alpha_target = M_PI;
    double omega_target = 0.0;

    matrix Y0(2, 1);
    Y0(0) = 0.0;
    Y0(1) = 0.0;

    matrix control_params(2, 1);
    control_params(0) = k1;
    control_params(1) = k2;

    matrix *result = solve_ode(df2R, 0.0, 0.1, 100.0, Y0, matrix(), control_params);

    double Q = 0.0;
    for (int i = 0; i < get_len(result[0]); ++i) {
        double alpha = result[1](i, 0);
        double omega = result[1](i, 1);
        double M = k1 * (alpha_target - alpha) + k2 * (omega_target - omega);
        Q += 10 * pow(alpha_target - alpha, 2) + pow(omega_target - omega, 2) + pow(M, 2);
    }

    delete[] result;

    return {Q};
}

matrix df3R(double t, matrix Y, matrix ud1, matrix ud2) {
    double C = 0.47;
    double r = 0.12;
    double m = 0.6;
    double ro = 1.2;
    double g = 9.81;

    double S = M_PI * r * r;

    double omega = m2d(ud1);

    double x = Y(0);
    double y = Y(1);
    double vx = Y(2);
    double vy = Y(3);

    // Si�y
    double Dx = 0.5 * C * ro * S * vx * abs(vx);
    double Dy = 0.5 * C * ro * S * vy * abs(vy);
    double Fx = ro * vy * omega * M_PI * r * r * r;
    double Fy = ro * vx * omega * M_PI * r * r * r;

    matrix dY(4, 1);
    dY(0) = vx;
    dY(1) = vy;
    dY(2) = -(Dx + Fx) / m;
    dY(3) = -(Dy + Fy + m * g) / m;

    return dY;
}

matrix ff3R(matrix x, matrix ud1, matrix ud2) {
    matrix y;
    matrix Y_temp(4, 1);
    Y_temp(0) = 0;
    Y_temp(1) = 100;
    Y_temp(2) = x(0);
    Y_temp(3) = 0;
    matrix *Y = solve_ode(df3R, 0, 0.01, 7, Y_temp, x(1));
    int n = get_len(Y[0]);
    int poczatek = 0, koniec = 0;

    for (int i = 0; i < n; i++) {
        if (abs(Y[1](i, 1) - 50) < abs(Y[1](koniec, 1) - 50))
            koniec = i;
        if (abs(Y[1](i, 1)) < abs(Y[1](poczatek, 1)))
            poczatek = i;
    }

    cout << Y[1](koniec, 0) << ";\t" << Y[1](koniec, 1) << ";\n";
    cout << Y[1](poczatek, 0) << ";\t" << Y[1](poczatek, 1) << ";\n";

    y = -Y[1](poczatek, 0);
    if (abs(x(0)) - 10 > 0)
        y = y + ud2 * pow(abs(x(0)) - 10, 2);
    if (abs(x(1)) - 15 > 0)
        y = y + ud2 * pow(abs(x(1)) - 15, 2);
    if (abs(Y[1](koniec, 0) - 5) > 0)
        y = y + ud2 * pow(abs(Y[1](koniec, 0) - 5), 2);

    return y;
}

matrix ff3T(matrix x, matrix ud1, matrix ud2) {
    double temp = M_PI * sqrt(pow(x(0) / M_PI, 2) + pow(x(1) / M_PI, 2));
    matrix y = sin(temp) / temp;

    //Zewn�trzna
    if (ud2(1) > 1) {
        if (-x(0) + 1 > 0)
            y = y + (ud2)(0) * pow(-x(0) + 1, 2);
        if (-x(1) + 1 > 0)
            y = y + (ud2)(0) * pow(-x(1) + 1, 2);
        if (norm(x) - (ud1)(0) > 0)
            y = y + (ud2)(0) * pow(norm(x) - (ud1)(0), 2);
    }
        //Wewn�trzna
    else {
        if (-x(0) + 1 > 0)
            y = 1e12;
        else
            y = y - (ud2)(0) / (-x(0) + 1);

        if (-x(1) + 1 > 0)
            y = 1e12;
        else
            y = y - (ud2)(0) / (-x(1) + 1);

        if (norm(x) - (ud1)(0) > 0)
            y = 1e12;
        else
            y = y - (ud2)(0) / (sqrt(pow(x(0), 2) + pow(x(1), 2)) - (ud1)(0));
    }
    return y;
}

matrix gradient(matrix x, matrix ud1, matrix ud2) {
    matrix gradient(2, 1);
    gradient(0) = 10 * x(0) + 8 * x(1) - 34;
    gradient(1) = 8 * x(0) + 10 * x(1) - 38;
    return gradient;
}

matrix hessian(matrix x, matrix ud1, matrix ud2) {
    matrix hessian(2, 2);
    hessian(0, 0) = 10;
    hessian(0, 1) = 8;
    hessian(1, 0) = 8;
    hessian(1, 1) = 10;
    return hessian;
}

matrix ff4T(matrix x, matrix ud1, matrix ud2) {
    matrix Y;
    if (isnan(ud2(0, 0)))
        Y = pow(x(0) + 2 * x(1) - 7, 2) + pow(2 * x(0) + x(1) - 5, 2);
    else
        Y = ff4T(ud2[0] + x * ud2[1], 0, ud1);
    return Y;
}

matrix ff4R(matrix x, matrix ud1, matrix ud2) {
    int numberOfRows = get_len(x);
    int numberOfColumns = 100;
    static matrix XData(numberOfRows, numberOfColumns), YData(1, numberOfColumns);
    if (solution::f_calls == 1) {
        fstream file;
        file.open(R"(C:\Users\Dell\Downloads\optymalizacja9000-main\XData.txt)");
        file >> XData;
        file.close();
        file.open(R"(C:\Users\Dell\Downloads\optymalizacja9000-main\YData.txt)");
        file >> YData;
        file.close();
    }
    double hypothesisClassificator;
    matrix Y = 0;
    for (int i = 0; i < numberOfColumns; i++) {
        hypothesisClassificator = m2d(trans(x) * XData[i]);
        hypothesisClassificator = 1.0 / (1.0 + exp(-hypothesisClassificator));
        Y = Y - YData(0, i) * log(hypothesisClassificator) - (1 - YData(0, i)) * log(1 - hypothesisClassificator);
    }
    Y = Y / numberOfColumns;
    return Y;
}

matrix gf(matrix x, matrix ud1, matrix ud2) {
    int numberOfRows = get_len(x);
    int numberOfColumns = 100;
    static matrix XData(numberOfRows, numberOfColumns), YData(1, numberOfColumns);
    if (solution::g_calls == 1) {
        fstream file;
        file.open(R"(C:\Users\Dell\Downloads\optymalizacja9000-main\XData.txt)");
        file >> XData;
        file.close();
        file.open(R"(C:\Users\Dell\Downloads\optymalizacja9000-main\YData.txt)");
        file >> YData;
        file.close();
    }
    double hypothesisClassificator;
    matrix result(numberOfRows, 1);
    for (int j = 0; j < numberOfRows; ++j) {
        for (int i = 0; i < numberOfColumns; ++i) {
            hypothesisClassificator = m2d(trans(x) * XData[i]);
            hypothesisClassificator = 1 / (1 + exp(-hypothesisClassificator));
            result(j) = result(j) + XData(j, i) * (hypothesisClassificator - YData(0, i));
        }
        result(j) = result(j) / numberOfColumns;
    }
    return result;
}