#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <string>
#include <cmath>
#include <limits>

// Структура для хранения координат точки и метки класса
struct Point {
    double x;
    double y;
    int label; // Метка класса: 1 или -1
};

// Структура для хранения области (прямоугольника)
struct Region {
    double x_min;
    double y_min;
    double x_max;
    double y_max;
};

// Функция для вычисления Евклидова расстояния между двумя точками // Добавлено
double euclideanDistance(const Point& p1, const Point& p2) {
    double dx = p1.x - p2.x;
    double dy = p1.y - p2.y;
    return std::sqrt(dx * dx + dy * dy);
}

// Функция для вычисления Манхэттенского расстояния между двумя точками // Добавлено
double manhattanDistance(const Point& p1, const Point& p2) {
    return std::fabs(p1.x - p2.x) + std::fabs(p1.y - p2.y);
}

// Функция для реализации алгоритма Хо-Кашьяпа
bool hoKashyap(const std::vector<Point>& points, std::vector<double>& w, double& b, int max_iterations = 10000, double tolerance = 1e-3) {
    int n = points.size();
    // Формируем матрицу A и вектор b_mat
    // A * [w; b] > 1
    // Здесь A будет матрицей размера n x 3
    std::vector<std::vector<double>> A_matrix(n, std::vector<double>(3, 0.0));
    std::vector<double> b_vec(n, 1.0); // Правая часть неравенств

    for (int i = 0; i < n; ++i) {
        A_matrix[i][0] = points[i].x * points[i].label;
        A_matrix[i][1] = points[i].y * points[i].label;
        A_matrix[i][2] = points[i].label;
    }

    // Инициализация параметров
    std::vector<double> x(3, 0.0); // [w1, w2, b]
    std::vector<double> e(n, 1.0); // Начальные ошибки

    double eta = 0.01; // Шаг обучения

    for (int iter = 0; iter < max_iterations; ++iter) {
        bool all_constraints_satisfied = true;
        for (int i = 0; i < n; ++i) {
            double y = A_matrix[i][0] * x[0] + A_matrix[i][1] * x[1] + A_matrix[i][2] * x[2];
            double error = 1.0 - y;
            if (error > tolerance) {
                all_constraints_satisfied = false;
                for (int j = 0; j < 3; ++j) {
                    x[j] += eta * A_matrix[i][j] * error;
                }
            }
        }
        if (all_constraints_satisfied) {
            w = { x[0], x[1] };
            b = x[2];
            return true;
        }
    }
    return false;
}

int main() {
    setlocale(LC_ALL, "Russian");
    const std::string input_file = "input.txt";
    const std::string output_file = "output.txt";
    const std::string boundary_file = "boundary.txt";

    std::ifstream fin(input_file);
    if (!fin.is_open()) {
        std::cerr << "Не удалось открыть входной файл: " << input_file << std::endl;
        return 1;
    }

    int num_regions;
    fin >> num_regions;
    if (num_regions <= 0) {
        std::cerr << "Количество областей должно быть положительным." << std::endl;
        return 1;
    }

    if (num_regions != 2) {
        std::cerr << "Алгоритм Хо-Кашьяпа реализован только для двух областей (двоичная классификация)." << std::endl;
        return 1;
    }

    std::vector<Region> regions;
    for (int i = 0; i < num_regions; ++i) {
        double x1, y1, x2, y2;
        fin >> x1 >> y1 >> x2 >> y2;
        double xmin = std::min(x1, x2);
        double ymin = std::min(y1, y2);
        double xmax = std::max(x1, x2);
        double ymax = std::max(y1, y2);
        regions.push_back(Region{ xmin, ymin, xmax, ymax });
    }
    fin.close();

    std::random_device rd;
    std::mt19937 gen(rd());
    std::vector<Point> points;
    for (size_t i = 0; i < regions.size(); ++i) {
        const auto& region = regions[i];
        std::uniform_int_distribution<> num_points_dist(20, 50);
        int num_points = num_points_dist(gen);
        std::uniform_real_distribution<> x_dist(region.x_min, region.x_max);
        std::uniform_real_distribution<> y_dist(region.y_min, region.y_max);
        for (int j = 0; j < num_points; ++j) {
            double x = x_dist(gen);
            double y = y_dist(gen);
            int label = (i == 0) ? 1 : -1;
            points.push_back(Point{ x, y, label });
        }
    }

    std::ofstream fout(output_file);
    if (!fout.is_open()) {
        std::cerr << "Не удалось открыть выходной файл: " << output_file << std::endl;
        return 1;
    }

    for (const auto& point : points) {
        fout << point.x << " " << point.y << " " << point.label << "\n";
    }
    fout.close();

    std::cout << "Генерация точек завершена. Результат записан в " << output_file << std::endl;

    if (points.size() < 2) {
        std::cerr << "Недостаточно точек для выполнения алгоритма Хо-Кашьяпа." << std::endl;
        return 1;
    }

    std::vector<double> w;
    double b;

    bool success = hoKashyap(points, w, b);

    if (success) {
        std::cout << "Алгоритм Хо-Кашьяпа успешно завершен." << std::endl;
        std::cout << "Разделяющая прямая: " << w[0] << "x + " << w[1] << "y + " << b << " = 0" << std::endl;

        std::ofstream fb(boundary_file);
        if (!fb.is_open()) {
            std::cerr << "Не удалось открыть файл для записи границы: " << boundary_file << std::endl;
            return 1;
        }

        double x_min = std::numeric_limits<double>::max();
        double x_max = std::numeric_limits<double>::lowest();
        for (const auto& point : points) {
            if (point.x < x_min) x_min = point.x;
            if (point.x > x_max) x_max = point.x;
        }

        double range_padding = (x_max - x_min) * 0.1;
        x_min -= range_padding;
        x_max += range_padding;

        if (fabs(w[1]) > 1e-6) {
            double y1 = (-w[0] * x_min - b) / w[1];
            double y2 = (-w[0] * x_max - b) / w[1];
            fb << x_min << " " << y1 << "\n";
            fb << x_max << " " << y2 << "\n";
        }
        else {
            double x = -b / w[0];
            fb << x << " " << -1e5 << "\n";
            fb << x << " " << 1e5 << "\n";
        }
        fb.close();
        std::cout << "Граница записана в " << boundary_file << std::endl;
    }
    else {
        std::cerr << "Алгоритм Хо-Кашьяпа не смог найти разделяющую прямую за отведенное количество итераций." << std::endl;
        return 1;
    }

    // Демонстрация использования новых функций расстояний // Добавлено
    if (points.size() >= 2) {
        const Point& p1 = points[0];
        const Point& p2 = points[1];
        double dist_euclid = euclideanDistance(p1, p2);
        double dist_manh = manhattanDistance(p1, p2);
        std::cout << "Расстояние между первой и второй точками (Евклидово): " << dist_euclid << std::endl;
        std::cout << "Расстояние между первой и второй точками (Манхэттенское): " << dist_manh << std::endl;
    }

    return 0;
}