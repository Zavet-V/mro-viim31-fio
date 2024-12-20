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

    // Параметры алгоритма
    double eta = 0.01; // Шаг обучения

    for (int iter = 0; iter < max_iterations; ++iter) {
        bool all_constraints_satisfied = true;
        for (int i = 0; i < n; ++i) {
            // Вычисляем y = A_i * x
            double y = A_matrix[i][0] * x[0] + A_matrix[i][1] * x[1] + A_matrix[i][2] * x[2];
            double error = 1.0 - y;
            if (error > tolerance) {
                all_constraints_satisfied = false;
                // Обновляем x
                for (int j = 0; j < 3; ++j) {
                    x[j] += eta * A_matrix[i][j] * error;
                }
            }
        }
        if (all_constraints_satisfied) {
            // Найдено решение
            w = { x[0], x[1] };
            b = x[2];
            return true;
        }
    }
    // Если не удалось найти решение за max_iterations
    return false;
}

int main() {
    setlocale(LC_ALL, "Russian");
    // Имена файлов
    const std::string input_file = "input.txt";
    const std::string output_file = "output.txt";
    const std::string boundary_file = "boundary.txt";

    // Открытие входного файла
    std::ifstream fin(input_file);
    if (!fin.is_open()) {
        std::cerr << "Не удалось открыть входной файл: " << input_file << std::endl;
        return 1;
    }

    // Чтение количества областей
    int num_regions;
    fin >> num_regions;

    // Проверка на корректность считанных значений
    if (num_regions <= 0) {
        std::cerr << "Количество областей должно быть положительным." << std::endl;
        return 1;
    }

    // Проверка, что количество областей равно 2 для двоичной классификации
    if (num_regions != 2) {
        std::cerr << "Алгоритм Хо-Кашьяпа реализован только для двух областей (двоичная классификация)." << std::endl;
        return 1;
    }

    // Чтение областей
    std::vector<Region> regions;
    for (int i = 0; i < num_regions; ++i) {
        double x1, y1, x2, y2;
        fin >> x1 >> y1 >> x2 >> y2;

        // Определение минимальных и максимальных координат
        double xmin = std::min(x1, x2);
        double ymin = std::min(y1, y2);
        double xmax = std::max(x1, x2);
        double ymax = std::max(y1, y2);

        regions.push_back(Region{ xmin, ymin, xmax, ymax });
    }

    fin.close(); // Закрытие входного файла

    // Настройка генератора случайных чисел
    std::random_device rd;
    std::mt19937 gen(rd());
    // Вектор для хранения сгенерированных точек
    std::vector<Point> points;
    // Генерация точек
    for (size_t i = 0; i < regions.size(); ++i) {
        const auto& region = regions[i];
        // Генерация случайного количества точек от 20 до 50 для каждой области
        std::uniform_int_distribution<> num_points_dist(20, 50);
        int num_points = num_points_dist(gen);
        // Распределение для координат внутри области
        std::uniform_real_distribution<> x_dist(region.x_min, region.x_max);
        std::uniform_real_distribution<> y_dist(region.y_min, region.y_max);

        for (int j = 0; j < num_points; ++j) {
            // Генерация точки
            double x = x_dist(gen);
            double y = y_dist(gen);

            // Присвоение метки класса: 1 для первой области, -1 для второй
            int label = (i == 0) ? 1 : -1;

            points.push_back(Point{ x, y, label });
        }
    }

    // Открытие выходного файла для точек
    std::ofstream fout(output_file);
    if (!fout.is_open()) {
        std::cerr << "Не удалось открыть выходной файл: " << output_file << std::endl;
        return 1;
    }

    // Запись точек в файл
    for (const auto& point : points) {
        fout << point.x << " " << point.y << " " << point.label << "\n";
    }

    fout.close(); // Закрытие выходного файла

    std::cout << "Генерация точек завершена. Результат записан в " << output_file << std::endl;

    // Проверка наличия достаточного количества точек для алгоритма
    if (points.size() < 2) {
        std::cerr << "Недостаточно точек для выполнения алгоритма Хо-Кашьяпа." << std::endl;
        return 1;
    }

    // Реализация алгоритма Хо-Кашьяпа
    std::vector<double> w; // Весовые коэффициенты
    double b; // Смещение

    bool success = hoKashyap(points, w, b);

    if (success) {
        std::cout << "Алгоритм Хо-Кашьяпа успешно завершен." << std::endl;
        std::cout << "Разделяющая прямая: " << w[0] << "x + " << w[1] << "y + " << b << " = 0" << std::endl;

        // Открытие файла для записи границы
        std::ofstream fb(boundary_file);
        if (!fb.is_open()) {
            std::cerr << "Не удалось открыть файл для записи границы: " << boundary_file << std::endl;
            return 1;
        }

        // Определение диапазона x для построения границы
        double x_min = std::numeric_limits<double>::max();
        double x_max = std::numeric_limits<double>::lowest();

        for (const auto& point : points) {
            if (point.x < x_min) x_min = point.x;
            if (point.x > x_max) x_max = point.x;
        }

        // Добавляем небольшой отступ для визуализации
        double range_padding = (x_max - x_min) * 0.1;
        x_min -= range_padding;
        x_max += range_padding;

        // Запись границы как (x1, y1) и (x2, y2)
        if (fabs(w[1]) > 1e-6) {
            double y1 = (-w[0] * x_min - b) / w[1];
            double y2 = (-w[0] * x_max - b) / w[1];
            fb << x_min << " " << y1 << "\n";  // (x_min, y1)
            fb << x_max << " " << y2 << "\n";  // (x_max, y2)
        }
        else {
            // Вертикальная прямая
            double x = -b / w[0];
            fb << x << " " << -1e5 << "\n"; // Очень маленькое y
            fb << x << " " << 1e5 << "\n";  // Очень большое y
        }

        fb.close();
        std::cout << "Граница записана в " << boundary_file << std::endl;
    }
    else {
        std::cerr << "Алгоритм Хо-Кашьяпа не смог найти разделяющую прямую за отведенное количество итераций." << std::endl;
        return 1;
    }

    return 0;
}