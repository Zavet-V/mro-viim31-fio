encoding='utf-8'

import matplotlib.pyplot as plt
import numpy as np

# Функция для чтения точек из файла
def read_points(file_name):
    points = []
    with open(file_name, 'r', encoding='utf-8') as file:
        for line in file: 
            x, y, label = map(float, line.split())
            points.append((x, y, label))
    return points

# Функция для чтения границы из файла
def read_boundary(file_name):
    boundary = []
    with open(file_name, 'r', encoding='utf-8') as file:
        for line in file:
            x, y = map(float, line.split())
            boundary.append((x, y))
    return boundary

# Загрузка точек из файла output.txt
points = read_points("output.txt")

# Загрузка границы из файла boundary.txt
boundary = read_boundary("boundary.txt")

# Разделение точек на два класса для визуализации
class_1_points = [(x, y) for x, y, label in points if label == 1]
class_2_points = [(x, y) for x, y, label in points if label == -1]

# Создание графика
plt.figure(figsize=(8, 6))

# Визуализация точек классов
class_1_points = np.array(class_1_points)
class_2_points = np.array(class_2_points)

plt.scatter(class_1_points[:, 0], class_1_points[:, 1], color='red', label='Class 1 (+1)', marker='o')
plt.scatter(class_2_points[:, 0], class_2_points[:, 1], color='blue', label='Class 2 (-1)', marker='x')

# Визуализация границы
boundary = np.array(boundary)
plt.plot(boundary[:, 0], boundary[:, 1], color='black', label='Decision Boundary')

# Настройки графика
plt.title('2D Visualization of Points and Decision Boundary')
plt.xlabel('X-axis')
plt.ylabel('Y-axis')
plt.legend(loc='best')
plt.grid(True)

# Показать график
plt.show()