
import numpy as np
import matplotlib.pyplot as plt 
import math
import interpolation
import NMclasses as nm

pi = math.pi
a = -3
b = 3


X = np.linspace(a, b, 500)
f = lambda x: np.cos(x)

x, y, splines = nm.interpolation(f, a, b, 2, 'Normal').splineInterpolation()
# x, y, splines = interpolation.splineInterpolation(f, a, b, 5, 'Normal')

# splines_dicts = list()
# dic = {'a': None, 'b': None, 'c': None, 'd': None}

# for spline in splines:
#     print(spline.create())
x1, y1, splines1 = nm.interpolation(f, a, b, 30, 'Cheb').splineInterpolation()
# x1, y1 = interpolation.splineInterpolation(f, a, b, 30, 'Cheb')
E, N  = interpolation.SPdepN(f, a, b, 4, 34, 'Cheb')
E1, N1 = interpolation.SPdepN(f, a, b, 4, 34, 'Normal')
Err, Nl = interpolation.LNdepN(x, f, a, b, 4, 34, 'Normal')
Err1, Nl1 = interpolation.LNdepN(x, f, a, b, 4, 34, 'Cheb')

nm.Plots(1, [[x,y]], 'График', 'x', 'f(x)').build()
nm.Plots(1, [[x, y], [x1,y1], [X, f(X)]], 'Интерполированный косинус', 'x', 'f(x)', ['Равномерная сетка', 'Сетка Чебышева', "cos(x)"]).build()
nm.Plots(2, [[N1,E1], [N, E]], 'Зависимость макс. ошибки от числа узлов', 'N', 'Err', ['Равномерная сетка сетка Сплайн', 'Сетка Чебышева Сплайн']).build()
plt.show()