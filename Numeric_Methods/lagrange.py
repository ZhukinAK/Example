from tkinter import *
import numpy as np
import matplotlib.pyplot as plt 
import math
import interpolation
import NMclasses as nm 


pi = math.pi
f = lambda x: (np.cos(x))

def lab(a, b, n):
    global entries
    a =  float(a)
    b = float(b)
    n = int(n)
    # a = 0
    # b = 2*pi
    x = np.linspace(a, b, 500)

    Xi, Yi, pn = nm.interpolation(f, a, b, n, 'Normal').lagrange()
    pchel = nm.interpolation(f, a, b, n, 'Cheb').lagrange()[2]
    pn = np.polyval(pn, x)
    pchel = np.polyval(pchel, x)
    [Err, N] = interpolation.LNdepN(x, f, a, b, 4, 26, 'Normal')
    Err1, N = interpolation.LNdepN(x, f, a, b, 4, 26, 'Cheb')
    Er, Len = interpolation.LNdepLen(f, -5, 5, 15, 'Normal')
    Er1, Len1 = interpolation.LNdepLen(f, -5, 5, 15, 'Cheb')
    nm.Plots(1, [[x,pn], [x, pchel], [x, f(x)]], 'Интерполированный косинус', 'x', 'f(x)', ['Нормальная сетка', 'Сетка Чебышева', "cos(x)"]).build()
    nm.Plots(2, [[Len, Er], [Len, Er1]], 'Зависимость макс. ошибки от длины интервала',  'Длина', 'Ошибка', ['Нормальная сетка', 'Сетка Чебышева']).build()
    nm.Plots(3, [[N,Err], [N, Err1]], 'Зависимость макс. ошибки от числа узлов', 'N', 'Err', ['Нормальная сетка', 'Сетка Чебышева']).build()
    
    plt.show()


window=Tk()
window.title('Welcome')
window.geometry('300x200')
Label(window, text = 'Plot').grid(column = 0, row = 0)
Label(window, text = 'Enter a').grid(column = 0, row = 1)
Label(window, text = 'Enter b').grid(column = 0, row = 2)
Label(window, text = 'Number of nodes').grid(column = 0, row = 3)
txt = Entry(window, width = 10)
txt.grid(row = 1, column = 1)
txt1 = Entry(window, width = 10)
txt1.grid(column = 1, row = 2)
txt2 = Entry(window, width = 10)
txt2.grid(column = 1, row = 3)
Button(window, text = 'Start!', command = (lambda: lab(txt.get(), txt1.get(), txt2.get()))).grid(row = 1, column = 2)
Button(window, text = 'Quit', command = (lambda: exit())).grid(row=2, column = 2)
window.mainloop()





# ## 2
# x1 = np.linspace(-2*pi, 2*pi, 500)
# pn1 = lagrange.interpolation(f, -2*pi, 2*pi, 15, 'Normal')[2]
# pchel1 = lagrange.interpolation(f, -2*pi, 2*pi, 15, 'Cheb')[2]
# pn1 = np.polyval(pn1, x1)
# pchel1 = np.polyval(pchel1, x1)

# ## 3
# x2 = np.linspace(-3*pi, 0, 500)
# pn2 = lagrange.interpolation(f, -3*pi, 0, 15, 'Normal')[2]
# pchel2 = lagrange.interpolation(f, -3*pi, 0, 15, 'Cheb')[2]
# pn2 = np.polyval(pn2, x2)
# pchel2 = np.polyval(pchel2, x2)

# plot3 = plots2D.Plots(4, [[x1,pn1], [x1, pchel1], [x1, f(x1)]], 'Интерполированный косинусб (Короткий симметричный промежуток)', 'x', 'f(x)', ['Нормальная сетка', 'Сетка Чебышева', "cos(x)"])
# plot3.build()

# plot4 = plots2D.Plots(5, [[x2,pn2], [x2, pchel2], [x2, f(x2)]], 'Интерполированный косинусб (Корокий несимметричный промежуток)', 'x', 'f(x)', ['Нормальная сетка', 'Сетка Чебышева', "cos(x)"])
# plot4.build()


