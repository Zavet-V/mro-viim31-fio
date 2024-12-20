#include <iostream> 
#include <fstream>
#include <vector>
#include <random>
#include <string>
#include <cmath>
#include <limits>
#include <algorithm>
#include <cstdlib>

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

// Функция для вычисления Евклидова расстояния между двумя точками
double euclideanDistance(const Point& p1, const Point& p2) {
    double dx = p1.x - p2.x;
    double dy = p1.y - p2.y;
    return std::sqrt(dx * dx + dy * dy);
}

// Функция для вычисления Манхэттенского расстояния между двумя точками
double manhattanDistance(const Point& p1, const Point& p2) {
    return std::fabs(p1.x - p2.x) + std::fabs(p1.y - p2.y);
}

// Функция для реализации алгоритма Хо-Кашьяпа (двоичная классификация)
bool hoKashyap(const std::vector<Point>& points, std::vector<double>& w, double& b, int max_iterations = 100000, double tolerance = 1e-3) {
    int n = (int)points.size();
    std::vector<std::vector<double>> A_matrix(n, std::vector<double>(3, 0.0));

    for (int i = 0; i < n; ++i) {
        A_matrix[i][0] = points[i].x * points[i].label;
        A_matrix[i][1] = points[i].y * points[i].label;
        A_matrix[i][2] = points[i].label;
    }

    std::vector<double> x(3, 0.0); // [w1, w2, b]
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

// --------------------------------------------------------------------------------------
// Реализация алгоритма FOREL
// --------------------------------------------------------------------------------------
std::vector<std::vector<Point>> forelClustering(std::vector<Point> points, double R, double epsilon = 1e-5) {
    std::vector<std::vector<Point>> clusters;

    while (!points.empty()) {
        // Выберем случайную точку в качестве начального центра
        int idx = rand() % points.size();
        Point center = points[idx];

        while (true) {
            // Найдём все точки, попадающие в радиус R
            std::vector<Point> inSphere;
            for (auto& p : points) {
                if (euclideanDistance(center, p) <= R) {
                    inSphere.push_back(p);
                }
            }

            if (inSphere.empty()) {
                // Если вдруг пусто (что маловероятно), выйдем
                break;
            }

            // Вычислим новый центр масс
            double sumX = 0.0, sumY = 0.0;
            for (auto& p : inSphere) {
                sumX += p.x;
                sumY += p.y;
            }
            Point newCenter;
            newCenter.x = sumX / inSphere.size();
            newCenter.y = sumY / inSphere.size();
            newCenter.label = 0; // Метка не важна для центра

            // Проверим сдвиг
            double shift = euclideanDistance(center, newCenter);
            center = newCenter;

            if (shift < epsilon) {
                // Конвергенция центра кластера
                clusters.push_back(inSphere);

                // Удаляем эти точки из consideration
                std::vector<Point> remaining;
                for (auto& p : points) {
                    bool isInCluster = false;
                    for (auto& cpoint : inSphere) {
                        if (std::fabs(p.x - cpoint.x) < 1e-14 && std::fabs(p.y - cpoint.y) < 1e-14) {
                            isInCluster = true;
                            break;
                        }
                    }
                    if (!isInCluster) {
                        remaining.push_back(p);
                    }
                }
                points = remaining;
                break;
            }
        }
    }
    return clusters;
}

// --------------------------------------------------------------------------------------
// Реализация алгоритма ISODATA (упрощённый вариант)
// --------------------------------------------------------------------------------------
struct Cluster {
    Point centroid;
    std::vector<Point> points;
};

std::vector<Cluster> isodataClustering(std::vector<Point> points, int k, int max_iterations = 100) {
    // Инициализируем k центров случайным образом
    std::vector<Cluster> clusters(k);
    for (int i = 0; i < k; ++i) {
        int idx = rand() % points.size();
        clusters[i].centroid = points[idx];
    }

    for (int iter = 0; iter < max_iterations; ++iter) {
        // Очистим кластеры
        for (auto& c : clusters) {
            c.points.clear();
        }

        // Присвоение каждой точки ближайшему центру
        for (auto& p : points) {
            double minDist = std::numeric_limits<double>::max();
            int closest = -1;
            for (int i = 0; i < k; ++i) {
                double d = euclideanDistance(p, clusters[i].centroid);
                if (d < minDist) {
                    minDist = d;
                    closest = i;
                }
            }
            clusters[closest].points.push_back(p);
        }

        // Пересчёт центров
        for (auto& c : clusters) {
            if (!c.points.empty()) {
                double sumX = 0.0, sumY = 0.0;
                for (auto& p : c.points) {
                    sumX += p.x;
                    sumY += p.y;
                }
                c.centroid.x = sumX / c.points.size();
                c.centroid.y = sumY / c.points.size();
            }
        }
    }
    return clusters;
}

// --------------------------------------------------------------------------------------
// Реализация персептрона
// --------------------------------------------------------------------------------------
bool perceptron(const std::vector<Point>& points, std::vector<double>& w, double& b, double learning_rate = 0.01, int max_iterations = 100000) {
    // Инициализируем веса и смещение нулями
    w = { 0.0, 0.0 };
    b = 0.0;

    for (int iter = 0; iter < max_iterations; ++iter) {
        bool no_error = true;

        // Проходим по всем точкам
        for (auto& p : points) {
            double prediction = w[0] * p.x + w[1] * p.y + b;
            int predicted_label = (prediction >= 0) ? 1 : -1;
            int error = p.label - predicted_label;

            if (error != 0) {
                // Обновляем веса и смещение
                w[0] += learning_rate * p.label * p.x;
                w[1] += learning_rate * p.label * p.y;
                b += learning_rate * p.label;
                no_error = false;
            }
        }

        if (no_error) {
            // Все классифицированы верно, персептрон сошёлся
            return true;
        }
    }
    // Если дошли до сюда, значит за max_iterations не сошлось
    return false;
}

// --------------------------------------------------------------------------------------
// Основная программа
// --------------------------------------------------------------------------------------
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
        std::cerr << "Для демонстрации классификации необходимо ровно две области." << std::endl;
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
        std::cerr << "Недостаточно точек для выполнения классификации." << std::endl;
        return 1;
    }

    // Выполним алгоритм Хо-Кашьяпа
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
        std::cerr << "Возможно, данные не являются линейно разделимыми." << std::endl;
    }

    // Реализация персептрона
    std::vector<double> w_percep;
    double b_percep;
    bool perceptron_success = perceptron(points, w_percep, b_percep);
    if (perceptron_success) {
        std::cout << "Персептрон успешно обучен." << std::endl;
        std::cout << "Разделяющая прямая персептрона: " << w_percep[0] << "x + " << w_percep[1] << "y + " << b_percep << " = 0" << std::endl;
    }
    else {
        std::cerr << "Персептрон не сошёлся за отведённое количество итераций." << std::endl;
        std::cerr << "Скорее всего данные не являются линейно разделимыми." << std::endl;
    }

    // Демонстрация использования функций расстояний
    if (points.size() >= 2) {
        const Point& p1 = points[0];
        const Point& p2 = points[1];
        double dist_euclid = euclideanDistance(p1, p2);
        double dist_manh = manhattanDistance(p1, p2);
        std::cout << "Расстояние между первой и второй точками (Евклидово): " << dist_euclid << std::endl;
        std::cout << "Расстояние между первой и второй точками (Манхэттенское): " << dist_manh << std::endl;
    }

    // Демонстрация работы алгоритма FOREL
    double R = 1.0; // Радиус для кластера
    auto forel_clusters = forelClustering(points, R);
    std::cout << "Количество кластеров, найденных FOREL: " << forel_clusters.size() << std::endl;

    // Демонстрация работы алгоритма ISODATA
    int k = 3;
    auto iso_clusters = isodataClustering(points, k);
    std::cout << "ISODATA: сформировано " << iso_clusters.size() << " кластеров." << std::endl;
    for (int i = 0; i < (int)iso_clusters.size(); ++i) {
        std::cout << "Кластер " << i + 1 << ": " << iso_clusters[i].points.size() << " точек\n";
    }

    return 0;
}
