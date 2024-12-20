import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np

# Чтение точек из файла output.txt
points = []
with open('output.txt', 'r') as f:
    for line in f:
        x, y = map(float, line.strip().split())
        points.append((x, y))

# Преобразование точек в массивы numpy
px, py = zip(*points)
px = np.array(px)
py = np.array(py)

# Вычисление координаты z (например, z = 0 для всех точек)
pz = np.zeros_like(px)

# Если хотите использовать функцию для z, например, z = sin(√(x² + y²)):
# pz = np.sin(np.sqrt(px**2 + py**2))

# Создание 3D-графика
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.scatter(px, py, pz, color='red')

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')

plt.title('3D визуализация точек из output.txt')
plt.show()
