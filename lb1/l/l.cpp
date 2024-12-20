#include <iostream>
#include <fstream>
#include <vector>
#include <random>
#include <string>

// Структура для хранения координат точки
struct Point {
    double x;
    double y;
};

// Структура для хранения области (прямоугольника)
struct Region {
    double x_min;
    double y_min;
    double x_max;
    double y_max;
};

int main() {
    setlocale(LC_ALL, "Russian");
    // Имена файлов
    const std::string input_file = "input.txt";
    const std::string output_file = "output.txt";

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
    for (const auto& region : regions) {
        // Генерация случайного количества точек от 2 до 20 для каждой области
        std::uniform_int_distribution<> num_points_dist(2, 20);
        int num_points = num_points_dist(gen);

        // Распределение для координат внутри области
        std::uniform_real_distribution<> x_dist(region.x_min, region.x_max);
        std::uniform_real_distribution<> y_dist(region.y_min, region.y_max);

        for (int i = 0; i < num_points; ++i) {
            // Генерация точки
            double x = x_dist(gen);
            double y = y_dist(gen);

            points.push_back(Point{ x, y });
        }
    }

    // Открытие выходного файла
    std::ofstream fout(output_file);
    if (!fout.is_open()) {
        std::cerr << "Не удалось открыть выходной файл: " << output_file << std::endl;
        return 1;
    }

    // Запись точек в файл
    for (const auto& point : points) {
        fout << point.x << " " << point.y << "\n";
    }

    fout.close(); // Закрытие выходного файла

    std::cout << "Генерация точек завершена. Результат записан в " << output_file << std::endl;

    return 0;
}