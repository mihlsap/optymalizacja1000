#include "opt_alg.h"
#include<iomanip>
#include <format>

void lab0();

void lab1();

void lab2();

void lab3();

void lab4();

void lab5();

void lab6();

int main() {
    try {
        lab6();
    } catch (string EX_INFO) {
        cerr << "ERROR:\n";
        cerr << EX_INFO << endl << endl;
    }
    system("pause");
    return 0;
}

void lab0() {
    //Funkcja testowa
    double epsilon = 1e-2;
    int Nmax = 10000;
    matrix lb(2, 1, -5), ub(2, 1, 5), a(2, 1);
    solution opt;
    a(0) = -1;
    a(1) = 2;
    opt = MC(ff0T, 2, lb, ub, epsilon, Nmax, a);
    cout << opt << endl << endl;
    solution::clear_calls();

    //Wahadlo
    Nmax = 1000;
    epsilon = 1e-2;
    lb = 0;
    ub = 5;
    double teta_opt = 1;
    opt = MC(ff0R, 1, lb, ub, epsilon, Nmax, teta_opt);
    cout << opt << endl << endl;
    solution::clear_calls();

    //Zapis symulacji do pliku csv
    matrix Y0 = matrix(2, 1), MT = matrix(2, new double[2]{m2d(opt.x), 0.5});
    matrix *Y = solve_ode(df0, 0, 0.1, 10, Y0, NAN, MT);
    ofstream Sout("symulacja_lab0.csv");
    Sout << hcat(Y[0], Y[1]);
    Sout.close();
    Y[0].~matrix();
    Y[1].~matrix();
}

void lab1() {
    try {
        // OPTYMALIZACJE

        string FILE_PATH = R"(C:\Users\Ala\Desktop\optyma\Optimalization.csv)";

        fstream csvFile;
        const int numOfElements = 100;
        double minVal = -15, maxVal = 15;
        set<double> uniqueNums;

        default_random_engine randomEngine(random_device{}());
        uniform_real_distribution<double> realDistribution(minVal, maxVal);

        while (uniqueNums.size() < numOfElements) {
            uniqueNums.insert(realDistribution(randomEngine));
        }

        for (double x0: uniqueNums) {
            // Expansion
            double d = 4.25, alpha = 1.00043;
            int Nmax = 1000;

            unique_ptr<double[]> p(expansion(ff1T, x0, d, alpha, Nmax));

            if (!p) {
                cout << "Expansion returned null for x0: " << x0 << endl;
                continue;
            }

            try {
                csvFile.open(FILE_PATH, ios::app);
                if (!csvFile.is_open()) {
                    throw runtime_error("Could not open file");
                }
                csvFile << x0 << "," << p[0] << "," << p[1] << "," << solution::f_calls << ",";
                csvFile.close();
            } catch (const exception &e) {
                cout << "Error writing expansion results: " << e.what() << endl;
                csvFile.close();
                continue;
            }
            solution::clear_calls();

            // Fibonacci
            try {
                solution Fibonacci = fib(ff1T, p[0], p[1], 0.0001);

                csvFile.open(FILE_PATH, ios::app);
                if (!csvFile.is_open()) {
                    throw runtime_error("Could not open file");
                }
                csvFile << m2d(Fibonacci.x) << "," << m2d(Fibonacci.y) << "," << solution::f_calls << ",";
                csvFile.close();
            } catch (const exception &e) {
                cout << "Error in Fibonacci: " << e.what() << endl;
                csvFile.close();
                continue;
            }
            solution::clear_calls();

            // Lagrange
            try {
                solution Lagrange = lag(ff1T, p[0], p[1], 0.0001, 1e-09, Nmax);

                csvFile.open(FILE_PATH, ios::app);
                if (!csvFile.is_open()) {
                    throw runtime_error("Could not open file");
                }
                csvFile << m2d(Lagrange.x) << "," << m2d(Lagrange.y) << "," << solution::f_calls << "\n";
                csvFile.close();
            } catch (const exception &e) {
                cout << "Error in Lagrange: " << e.what() << endl;
                csvFile.close();
                continue;
            }
            solution::clear_calls();
        }

        // SYMULACJA
        double DA0 = 0.001, d = 0.001, alpha = 2.257, epsilon = 1e-10, gamma = 1e-10;
        int Nmax = 10000;

        unique_ptr<double[]> range(expansion(ff1R, DA0, d, alpha, Nmax));

        if (!range) {
            throw runtime_error("Expansion for simulation returned null");
        }

        solution::clear_calls();

        solution FibonacciSolution = fib(ff1R, range[0], range[1], epsilon);
        cout << "Optimal DA value for Fibonacci's method: " << FibonacciSolution.x << "\n";
        cout << "Value of result function: " << FibonacciSolution.y << "\n";
        cout << "Number of calls of result function: " << solution::f_calls << "\n\n";
        solution::clear_calls();

        solution LagrangeSolution = lag(ff1R, range[0], range[1], epsilon, gamma, Nmax);
        cout << "Optimal DA value for Lagrange's method: " << LagrangeSolution.x << "\n";
        cout << "Value of result function: " << LagrangeSolution.y << "\n";
        cout << "Number of calls of result function: " << solution::f_calls << "\n";
        solution::clear_calls();

        double Pa = 0.5, Va0 = 5, Vb0 = 1, Tb0 = 20, t_0 = 0, t_step = 1, t_end = 2000;
        matrix Y0 = matrix(3, 1, Va0);
        Y0(1) = Vb0;
        Y0(2) = Tb0;

        double Da = m2d(FibonacciSolution.x);
        unique_ptr<matrix[]> FibonacciSimulation(solve_ode(df1R, t_0, t_step, t_end, Y0, Da, Pa));

        FILE_PATH =
                R"(C:\Users\Dell\Desktop\stuff\studia\semestr 5\Optymalizacja\project\data\data1\FibonacciSimulation.csv)";
        csvFile.open(FILE_PATH, ios::app);
        if (!csvFile.is_open()) {
            throw runtime_error("Could not open file for Fibonacci simulation");
        }
        csvFile << FibonacciSimulation[1] << "\n";
        csvFile.close();

        solution::clear_calls();

        Da = m2d(LagrangeSolution.x);
        unique_ptr<matrix[]> LagrangeSimulation(solve_ode(df1R, t_0, t_step, t_end, Y0, Da, Pa));

        FILE_PATH =
                R"(C:\Users\Dell\Desktop\stuff\studia\semestr 5\Optymalizacja\project\data\data1\LagrangeSimulation.csv)";
        csvFile.open(FILE_PATH, ios::app);
        if (!csvFile.is_open()) {
            throw runtime_error("Could not open file for Lagrange simulation");
        }
        csvFile << LagrangeSimulation[1] << "\n";
        csvFile.close();

        solution::clear_calls();
    } catch (const exception &e) {
        cout << "Fatal error in lab1: " << e.what() << endl;
    } catch (...) {
        cout << "Unknown error in lab1" << endl;
    }
}


void lab2() {
    try {
        std::string FILE_PATH = R"(C:\Users\Ala\Desktop\optyma\Optimalization.csv)";
        std::ofstream csvFile(FILE_PATH, std::ios::out);
        if (!csvFile.is_open()) throw std::runtime_error("Could not open file");

        double minVal = -1.0, maxVal = 1.0;
        std::default_random_engine randomEngine(std::random_device{}());
        std::uniform_real_distribution<double> realDistribution(minVal, maxVal);

        double epsilon = 1e-6;
        int Nmax = 1000;
        matrix ud1, ud2;
        std::vector<double> step_sizes = {0.44, 0.22, 0.11};


        matrix x0(2, 1);
        x0(0, 0) = realDistribution(randomEngine);
        x0(1, 0) = realDistribution(randomEngine);


        solution Xopt_HJ = HJ(ff2T, x0, 0.22, 0.45, epsilon, Nmax, ud1, ud2);
        solution::clear_calls();

        matrix s0(2, 1, 0.22);
        solution Xopt_Rosen = Rosen(ff2T, x0, s0, 1.22, 0.55, epsilon, Nmax, ud1, ud2);
        solution::clear_calls();


        csvFile.close();
    } catch (const std::exception &e) {
        std::cerr << "Fatal error in lab2: " << e.what() << std::endl;
    } catch (...) {
        std::cerr << "Unknown error in lab2" << std::endl;
    }
}

void lab3() {
    //FUNKCJA TESTOWA
    double penalty_start = 1.23;
    double penalty_scale_ext = 2.0;
    double penalty_scale_int = 1.5;
    double epsilon = 1e-3;
    int Nmax = 10000;
    int num_iterations = 100;
    std::vector<double> a_values = {4.0, 4.4934, 5.0};

    std::ofstream results("C:/Users/Ala/Desktop/optyma/results.csv");
    results << std::fixed << std::setprecision(7);
    results << "a;x_1;x_2;x1_ext;x2_ext;r_ext;y_ext;calls_ext;x1_int;x2_int;r_int;y_int;calls_int\n";

    for (double a_val: a_values) {
        matrix a_value = matrix(1, 1, a_val);

        for (int i = 0; i < num_iterations; ++i) {
            matrix x0 = rand_mat(2, 1);
            x0(0, 0) = x0(0, 0) * 2.0 - 1.0;
            x0(1, 0) = x0(1, 0) * 2.0 - 1.0;

            //Optymalizacja zewnêtrzna
            solution::clear_calls();
            solution external_opt = pen(ff3T, x0, penalty_start, penalty_scale_ext, epsilon, Nmax, a_value,
                                        matrix(2, new double[2]{1.0, 2.0}));
            double x1_ext = external_opt.x(0, 0);
            double x2_ext = external_opt.x(1, 0);
            double r_ext = norm(external_opt.x);
            double y_ext = external_opt.y(0, 0);
            int calls_ext = solution::f_calls;

            //Optymalizacja wewnêtrzna
            solution::clear_calls();
            solution internal_opt = pen(ff3T, x0, penalty_start, penalty_scale_int, epsilon, Nmax, a_value,
                                        matrix(2, new double[2]{1.0, 0.5}));
            double x1_int = internal_opt.x(0, 0);
            double x2_int = internal_opt.x(1, 0);
            double r_int = norm(internal_opt.x);
            double y_int = internal_opt.y(0, 0);
            int calls_int = solution::f_calls;

            // Zapis wyników w jednej linii
            results << a_val << ";"
                    << x0(0, 0) << ";" << x0(1, 0) << ";"
                    << x1_ext << ";" << x2_ext << ";" << r_ext << ";" << y_ext << ";" << calls_ext << ";"
                    << x1_int << ";" << x2_int << ";" << r_int << ";" << y_int << ";" << calls_int << "\n";
        }
    }

    results.close();
    std::cout << "Optymalizacja zakoñczona, wyniki zapisano w pliku results.csv" << std::endl;

    //SYMULACJA PROBLEMU
    matrix x_1 = matrix(2, 1);
    double c0 = 1.23;
    x_1(0) = 0.;     //vx_0  [-10, 10] (m/s)
    x_1(1) = 0.;     //omega [-15, 15] (rad/s)

    x_1 = pen(ff3R, x_1, c0, 2, epsilon, Nmax).x;
    cout << endl;

    ff3R(x_1, c0, 2);
}

void lab4() {
    const double epsilon = 1e-6;
    const int Nmax = 10000;
    matrix points = matrix(2, 1);

    double h[3] = {0.05, 0.12, -1};
    for (auto &value: h) {
        fstream file;
        std::string path = std::format(R"(C:\Users\Dell\Downloads\optymalizacja9000-main\data_{}.csv)", value);
        file.open(path, ios::out);

        for (int i = 0; i < 100; i++) {
            points = 20 * rand_mat(2, 1) - 10;
            file << points(0) << ";" << points(1) << ";";

            solution SD_sol = SD(ff4T, gradient, points, value, epsilon, Nmax);
            file << SD_sol.x(0) << ";" << SD_sol.x(1) << ";" << SD_sol.y[0] << SD_sol.f_calls << ";" << SD_sol.g_calls
                 << ";";
            solution::clear_calls();

            solution CG_sol = CG(ff4T, gradient, points, value, epsilon, Nmax);
            file << CG_sol.x(0) << ";" << CG_sol.x(1) << ";" << CG_sol.y[0] << CG_sol.f_calls << ";" << CG_sol.g_calls
                 << ";";
            solution::clear_calls();

            solution Nw_sol = Newton(ff4T, gradient, hessian, points, value, epsilon, Nmax);
            file << Nw_sol.x(0) << ";" << Nw_sol.x(1) << ";" << Nw_sol.y[0] << Nw_sol.f_calls << ";" << Nw_sol.g_calls
                 << ";" << Nw_sol.H_calls << ";";
            solution::clear_calls();

            file << "\n";
        }
        file.close();

        changeSign(path, ',', '.');
        changeSign(path, ';', ',');
    }
//    points(0) = 6.2231;
//    points(1) = -0.238002;
//
//    Newton(ff4T, gradient, hessian, points, -1, epsilon, Nmax);
//
//    matrix x1(3, 1, 0.0);
//
//    solution Real_sol = CG(ff4R, gf, x1, 0.001, 0.000001, Nmax);
//
//    int m = 100;
//    static matrix X(3, m), Y(1, m);
//    ifstream data(R"(C:\Users\Dell\Downloads\optymalizacja9000-main\XData.txt)");
//    data >> X;
//    data.close();
//    data.open(R"(C:\Users\Dell\Downloads\optymalizacja9000-main\YData.txt)");
//    data >> Y;
//    data.close();
//
//    double P = 0.0;
//
//    for (int i = 0; i < 100; i++) {
//        double h0 = 1.0 / (1 + exp(-(trans(Real_sol.x) * X[i])()));
//        if (lround(h0) == Y(0, i))
//            h0 = 1;
//        else
//            h0 = 0;
//
//        P += h0;
//    }
//    P /= m;
//
//    cout << Real_sol.x(0, 0) << " " << Real_sol.x(1, 0) << " " << Real_sol.x(2, 0) << " " << Real_sol.y(0, 0) << " " << P << " " << Real_sol.g_calls << endl;
}

void lab5() {
    solution results;
    const int iterations = 101;
    vector<double> a = {1.0, 10.0, 100.0};
    double w = 0.0;
    matrix ud1(2, new double[2]{w, a[0]});

    const double epsilon = 1e-6;
    const int Nmax = 10000;

    matrix x0[iterations];
    for (auto &x: x0)
        x = 20 * rand_mat(2, 1) - 10;

    for (auto &value: a) {

        fstream file;
        std::string path = std::format(R"(C:\Users\Dell\Downloads\optymalizacja9000-main\data5\data_{}.csv)", value);
        file.open(path, ios::out);

        for (auto &x: x0) {

            ud1(1) = value;
            results = Powell(ff5T, x, epsilon, Nmax, ud1, NAN);

            if (value == a[0])
                file << x(0) << ";" << x(1) << ";";

            file << results.x(0) << ";" << results.x(1) << ";" << results.y(0) << ";"
                 << results.y(1) << ";" << solution::f_calls << "\n";

            solution::clear_calls();
            ud1(0) += 0.01;
        }
        file.close();
        changeSign(path, ',', '.');
        changeSign(path, ';', ',');
    }

    // Problem rzeczywisty
    fstream file;
    std::string path = R"(C:\Users\Dell\Downloads\optymalizacja9000-main\data5\real_solution_data.csv)";
    file.open(path, ios::out);

    solution result;
    matrix ud1_1(0.0);

    for (int i = 0; i < iterations; i++) {
        double l = (900.0 * ((double) rand() / (double) RAND_MAX)) + 100.0;
        double d = (40.0 * ((double) rand() / (double) RAND_MAX)) - 10.0;
        matrix x(2, new double[2]{l, d});
        result = Powell(ff5R, x, epsilon, Nmax, ud1_1, 0);
        file << l << " " << d << " " << result.x(0) << " " << result.x(1) << " " << result.y(0) << " "
             << result.y(1) << " " << " " << solution::f_calls << endl;
        solution::clear_calls();

        ud1_1(0) += 0.01;
    }

    file.close();
}

void lab6() {
    int N = 2, Nmax = 10000;
    int mi = 20, lambda = 40;
    double epsilon = 1e-5;

    vector<double> sigma = {0.01, 0.1, 1.0, 10.0, 100.0};

    matrix lb(1, N), ub(1, N); //lower i upper boundary
    lb(0, 0) = -5;
    lb(0, 1) = -5;
    ub(0, 0) = 5;
    ub(0, 1) = 5;

    for (auto &x : sigma) {
        fstream file;
        std::string path = std::format(R"(C:\Users\Dell\Downloads\optymalizacja9000-main\data6\data_{}.csv)", x);
        file.open(path, ios::out);

        matrix sigma0(2, 1, x);
        for (int i = 0; i < 100; i++) {
            solution results = EA(ff6T, N, lb, ub, mi, lambda, sigma0, epsilon, Nmax, NAN, NAN);
            file << results.x(0) << ";" << results.x(1) << ";" << results.y << solution::f_calls << "\n";
            changeSign(path, ',', '.');
            changeSign(path, ';', ',');
            solution::clear_calls();
        }
        file.close();
    }

//    matrix sigma0(2, 1, 1);
//
//    solution solEvo = EA(ff6T, N, lb, ub, mi, lambda, sigma0, epsilon, Nmax, NULL, NULL);
//
//    sigma0(0) = sigma0(1) = 100; // sigma - 0.01 / 0.1 / 1 / 10 / 100
//    for (int i = 0; i < 100; ++i) {
//        solution::clear_calls();
//        solEvo = EA(ff6T, N, lb, ub, mi, lambda, sigma0, epsilon, Nmax, NULL, NULL);
//
//        cout << solEvo.x(0) << ";" << solEvo.x(1) << ";" << solEvo.y << solEvo.f_calls << ";" << endl;
//    }


    // Problem rzeczywisty
//    matrix limits2(2, 2), sigma1(2, 1);
//    limits2(0, 0) = limits2(1, 0) = 0.1;
//    limits2(0, 1) = limits2(1, 1) = 3;
//
//    sigma1(0) = sigma1(1) = 1;
//
//    matrix simulation = matrix(0);
//    solution test;
//
//    // rzeczywisty przeprowadzono dla epsilon = 1e-4
//    //test = EA(ff6R, N, limits, mi, lambda, sigma1, epsilon, Nmax, simulation);
//
//    //cout << test << endl;
//    matrix xopt(2, 1);
//    xopt(0) = 2.10424;
//    xopt(1) = 0.0142826;
//
//    ff6R(xopt, NAN, NAN);


}