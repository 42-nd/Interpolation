#include "gen_grid.h"
#include "Cubic_interpolation_spline_1D.h"
#include <iostream> 
#include <cmath>
# define M_PI 3.14159265358979323846

std::vector<Com_Methods::Point> generate_intermediate_points(const std::vector<Com_Methods::Point>& grid, int num_points) {
    std::vector<Com_Methods::Point> intermediate_points;
    for (size_t i = 0; i < grid.size() - 1; ++i) {
        double x_start = grid[i].x();
        double x_end = grid[i + 1].x();
        double step = (x_end - x_start) / (num_points + 1);

        for (int j = 1; j <= num_points; ++j) {
            intermediate_points.push_back(Com_Methods::Point(x_start + step * j, 0.0, 0.0));
        }
    }
    return intermediate_points;
}

int main() {
    double a = -M_PI / 2;
    double b = 0;
    int parts = 10; // количество частей
    double sparse = 1; 

    std::vector<Com_Methods::Point> grid_h = grid_generator(a, b, parts, sparse);

    std::cout << "Grid: h = PI/20" << std::endl;
    std::cout << "X|f(x)" << std::endl;
    std::vector<double> values;
    for (const auto& point : grid_h) {
        std::cout << point.x() << " " << cos(point.x()) + sin(point.x()) << std::endl;
        values.push_back(cos(point.x()) + sin(point.x()));
    }
    std::cout << "                                        " << std::endl;

    Com_Methods::Cubic_Interpolation_Spline_1D spline;
    spline.Update_Spline(grid_h, values);
    std::vector<Com_Methods::Point> intermediate_points_h = generate_intermediate_points(grid_h, 1);

    std::vector<double> errors_f;
    std::vector<double> errors_f_prime;
    std::vector<double> errors_f_double_prime;
    double max_error_f = -1 , max_error_f_prime = -1, max_error_f_double_prime = -1;
    double result[3];
    for (const auto& point : intermediate_points_h) {
        spline.Get_Value(point, result);
        double x = point.x();
        double f_x = cos(x) + sin(x);
        double f_prime_x = -sin(x) + cos(x);
        double f_double_prime_x = -cos(x) - sin(x);
        errors_f.push_back(fabs(f_x - result[0]));
        errors_f_prime.push_back(fabs(f_prime_x - result[1]));
        errors_f_double_prime.push_back(fabs(f_double_prime_x - result[2]));
        max_error_f = std::max(max_error_f, fabs(f_x - result[0]));
        max_error_f_prime = std::max(max_error_f_prime, fabs(f_prime_x - result[1]));
        max_error_f_double_prime = std::max(max_error_f_double_prime, fabs(f_double_prime_x - result[2]));
        std::cout << "x = " << x << ", f(x) = " << f_x << ", f'(x) = " << f_prime_x << ", f''(x) = " << f_double_prime_x << ", s(x) = " << result[0] << ", s'(x) = " << result[1] << ", s''(x) = " << result[2] << std::endl;
    }
    std::cout << "                                        " << std::endl;
    std::cout << "Errors:" << std::endl;
    std::cout << "|f - s|       |f' - s'|    |f'' - s''|" << std::endl;
    for (int i = 0; i < errors_f.size(); i++) {
        std::cout << errors_f[i] << " | " << errors_f_prime[i] << " | " << errors_f_double_prime[i] << std::endl;
    }
    std::cout << "Max:" << std::endl;
    std::cout << max_error_f << " | " << max_error_f_prime << " | " << max_error_f_double_prime << std::endl;

    std::cout << "-----------------------------------------------" << std::endl;


    values.clear();
    errors_f.clear();
    errors_f_prime.clear();
    errors_f_double_prime.clear();
    max_error_f = -1;
    max_error_f_prime = -1;
    max_error_f_double_prime = -1;

    std::vector<Com_Methods::Point> grid_h2 = grid_generator(a, b, parts * 2, sparse);

    std::cout << "Grid: h = PI/40" << std::endl;
    std::cout << "X|f(x)" << std::endl;
    for (const auto& point : grid_h2) {
        std::cout << point.x() << " " << cos(point.x()) + sin(point.x()) << std::endl;
        values.push_back(cos(point.x()) + sin(point.x()));
    }
    std::cout << "                                        " << std::endl;

    Com_Methods::Cubic_Interpolation_Spline_1D spline_h2;
    spline_h2.Update_Spline(grid_h2, values);
    std::vector<Com_Methods::Point> intermediate_points_h2 = generate_intermediate_points(grid_h2, 1);

    for (const auto& point : intermediate_points_h2) {
        spline.Get_Value(point, result);
        double x = point.x();
        double f_x = cos(x) + sin(x);
        double f_prime_x = -sin(x) + cos(x);
        double f_double_prime_x = -cos(x) - sin(x);
        errors_f.push_back(fabs(f_x - result[0]));
        errors_f_prime.push_back(fabs(f_prime_x - result[1]));
        errors_f_double_prime.push_back(fabs(f_double_prime_x - result[2]));
        max_error_f = std::max(max_error_f, fabs(f_x - result[0]));
        max_error_f_prime = std::max(max_error_f_prime, fabs(f_prime_x - result[1]));
        max_error_f_double_prime = std::max(max_error_f_double_prime, fabs(f_double_prime_x - result[2]));
        std::cout << "x = " << x << ", f(x) = " << f_x << ", f'(x) = " << f_prime_x << ", f''(x) = " << f_double_prime_x << ", s(x) = " << result[0] << ", s'(x) = " << result[1] << ", s''(x) = " << result[2] << std::endl;
    }
    std::cout << "                                        " << std::endl;
    std::cout << "Errors:" << std::endl;
    std::cout << "|f - s|       |f' - s'|    |f'' - s''|" << std::endl;
    for (int i = 0; i < errors_f.size(); i++) {
        std::cout << errors_f[i] << " | " << errors_f_prime[i] << " | " << errors_f_double_prime[i] << std::endl;
    }
    std::cout << "Max:" << std::endl;
    std::cout << max_error_f << " | " << max_error_f_prime << " | " << max_error_f_double_prime << std::endl;

    std::cout << "-----------------------------------------------" << std::endl;


    values.clear();
    errors_f.clear();
    errors_f_prime.clear();
    errors_f_double_prime.clear();
    max_error_f = -1;
    max_error_f_prime = -1;
    max_error_f_double_prime = -1;

    std::vector<Com_Methods::Point> grid_h4 = grid_generator(a, b, parts * 4, sparse);

    std::cout << "Grid: h = PI/80" << std::endl;
    std::cout << "X|f(x)" << std::endl;
    for (const auto& point : grid_h4) {
        std::cout << point.x() << " " << cos(point.x()) + sin(point.x()) << std::endl;
        values.push_back(cos(point.x()) + sin(point.x()));
    }
    std::cout << "                                        " << std::endl;

    Com_Methods::Cubic_Interpolation_Spline_1D spline_h4;
    spline_h4.Update_Spline(grid_h4, values);
    std::vector<Com_Methods::Point> intermediate_points_h4 = generate_intermediate_points(grid_h4, 1);

    for (const auto& point : intermediate_points_h4) {
        spline.Get_Value(point, result);
        double x = point.x();
        double f_x = cos(x) + sin(x);
        double f_prime_x = -sin(x) + cos(x);
        double f_double_prime_x = -cos(x) - sin(x);
        errors_f.push_back(fabs(f_x - result[0]));
        errors_f_prime.push_back(fabs(f_prime_x - result[1]));
        errors_f_double_prime.push_back(fabs(f_double_prime_x - result[2]));
        max_error_f = std::max(max_error_f, fabs(f_x - result[0]));
        max_error_f_prime = std::max(max_error_f_prime, fabs(f_prime_x - result[1]));
        max_error_f_double_prime = std::max(max_error_f_double_prime, fabs(f_double_prime_x - result[2]));
        std::cout << "x = " << x << ", f(x) = " << f_x << ", f'(x) = " << f_prime_x << ", f''(x) = " << f_double_prime_x << ", s(x) = " << result[0] << ", s'(x) = " << result[1] << ", s''(x) = " << result[2] << std::endl;
    }
    std::cout << "                                        " << std::endl;
    std::cout << "Errors:" << std::endl;
    std::cout << "|f - s|       |f' - s'|    |f'' - s''|" << std::endl;
    for (int i = 0; i < errors_f.size(); i++) {
        std::cout << errors_f[i] << " | " << errors_f_prime[i] << " | " << errors_f_double_prime[i] << std::endl;
    }
    std::cout << "Max:" << std::endl;
    std::cout << max_error_f << " | " << max_error_f_prime << " | " << max_error_f_double_prime << std::endl;
    std::cout << "-----------------------------------------------" << std::endl;
    values.clear();
    grid_h = grid_generator(a, M_PI/2, parts, sparse);
    for (const auto& point : grid_h) {
        values.push_back(cos(point.x()) + sin(point.x()));
    }
    spline.Update_Spline(grid_h, values);

    double temp[3];
     // Четвертый пункт: табличные значения для x = -π/2, 0, π/2
    std::vector<double> special_points = { -M_PI / 2, 0, M_PI / 2 };
    std::vector<double> special_f_1, special_f_2; 

    std::cout << "\nPoint values -pi/2, 0, pi/2:\n";
    for (size_t i = 0; i < special_points.size(); ++i) {
        double x_k = special_points[i];
        Com_Methods::Point point(x_k, 0.0, 0.0); 
        spline.Get_Value(point, temp); 
        std::cout << "f(" << x_k << ") = " << cos(x_k) + sin(x_k)
            << ", s(" << x_k << ") = " << temp[0] << ", s'(" << x_k << ") = " << temp[1] << std::endl;

        if (i < special_points.size() - 1) {
            // Формула f'(x_k) ≈ (f(x_{k+1}) - f(x_k)) / (x_{k+1} - x_k)
            double f_xkp1 = cos(special_points[i + 1]) + sin(special_points[i + 1]);
            double f_xk = cos(x_k) + sin(x_k);
            double x_kp1 = special_points[i + 1];
            special_f_1.push_back((f_xkp1 - f_xk) / (x_kp1 - x_k));
        }
        if (i > 0){
            // Формула f'(x_k) ≈ (f(x_k) - f(x_{k-1})) / (x_k - x_{k-1})
            double f_xkm1 = cos(special_points[i - 1]) + sin(special_points[i - 1]);
            double f_xk = cos(x_k) + sin(x_k);
            double x_km1 = special_points[i - 1];
            special_f_2.push_back((f_xk - f_xkm1) / (x_k - x_km1));
        }
    }
    std::cout << "\nDerivatives by formulas:\n";
    for (size_t i = 0; i < special_f_1.size() ; ++i) {
        std::cout << "f'(" << special_points[i] << ") = "<< special_f_1[i] << " (first formula)"<< std::endl;
    }
    for (size_t i = 0; i < special_f_2.size(); ++i) {
        std::cout << "f'(" << special_points[i+1] << ") = " << special_f_2[i] <<  " (second formula)" << std::endl;
    }

    return 0;
}