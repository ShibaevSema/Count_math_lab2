import decimal

import matplotlib.pyplot as plt
from prettytable import PrettyTable
import numpy as np
from decimal import Decimal, localcontext


def setStart(a, b, x3, x2, x1, k):
    if (func2_check(a, x3, x2, x1, k)):
        return a
    else:
        return b


def func2(x3, x2, x1, k, a, b, eps):
    if (f(a, x3, x2, x1, k) * f(b, x3, x2, x1, k) < 0):
        last_X = 0
        pt = PrettyTable(['Номер шага', 'a', 'b', 'x', 'F(a)', 'F(b)', 'F(x)', '| a - b |'])
        n = 0
        ab = eps
        current_X = a - (f(a, x3, x2, x1, k) * (b - a) / (f(b, x3, x2, x1, k) - f(a, x3, x2, x1, k)))
        while (ab >= eps or abs(f(current_X, x3, x2, x1, k) >= eps)):
            if (f(a, x3, x2, x1, k) * f(current_X, x3, x2, x1, k) < 0):
                b = current_X
            else:
                a = current_X
            last_X = current_X
            current_X = a - f(a, x3, x2, x1, k) * ((b - a) / (f(b, x3, x2, x1, k) - f(a, x3, x2, x1, k)))
            ab = abs(current_X - last_X)
            pt.add_row([n, a, b, current_X, f(a, x3, x2, x1, k), f(b, x3, x2, x1, k),
                        f(current_X, x3, x2, x1, k), ab])
            n = n + 1
    else:
        pt = "Сходимость не выполняется !"
    return pt


def func2_check(a, x3, x2, x1, k):
    return f(a, x3, x2, x1, k) * fDoubleDerivative(a, x3, x2, x1, k) > 0


def func4(x3, x2, x1, k, a, b, eps):
    pt = PrettyTable(['Xk-1', 'F(Xk-1)', 'Xk', 'F(Xk)', 'Xk+1', 'F(Xk+1)', '| Xk - Xk+1 |', "номер итерации"])
    last_X = setStart(a, b, x3, x2, x1, k)
    current_X = last_X + 0.00001
    next_X = current_X - (
            f(current_X, x3, x2, x1, k) * (current_X - last_X) / (f(current_X, x3, x2, x1, k) - f(last_X, x3,
                                                                                                  x2, x1, k)))
    ab = abs(current_X - next_X)
    n = 0
    while (ab >= eps or abs(f(next_X, x3, x2, x1, k)) >= eps):
        last_X = current_X
        current_X = next_X
        next_X = current_X - (
                f(current_X, x3, x2, x1, k) * (current_X - last_X) / (f(current_X, x3, x2, x1, k) - f(last_X, x3,
                                                                                                      x2, x1, k)))
        ab = abs(current_X - next_X)
        n = n + 1
        pt.add_row([last_X, f(last_X, x3, x2, x1, k), current_X, f(current_X, x3, x2, x1, k),
                    next_X, f(next_X, x3, x2, x1, k), ab, n])
    return pt


def func5(x3, x2, x1, k, a, b, eps):
    lambd = getLambda(a, b, x3, x2, x1)
    print(lambd)
    af = fiDerivative(a, x3, x2, x1, lambd)
    bf = fiDerivative(b, x3, x2, x1, lambd)
    print("fi'(a)= " + str(af) + '\n' + "fi'(b) = " + str(bf) + '\n')
    pt = PrettyTable(['Ответ', '|X-Xo|', 'F(x)', 'Количество итераций', ])
    count = 0
    if (abs(fiDerivative(a, x3, x2, x1, lambd)) < 1 and abs(fiDerivative(b, x3, x2, x1, lambd)) < 1):
        if (f1(a, x3, x2, x1) < f1(b, x3, x2, x1)):
            x0 = b
        else:
            x0 = a

        x = fi(x0, x3, x2, x1, k, lambd)
        dif = abs(x - x0)
        while ( f(x, x3, x2, x1, k) > eps):
            x0 = x
            x = fi(x0, x3, x2, x1, k, lambd)
            count += 1
            dif = abs(x - x0)
            pt.add_row([x, abs(x - x0), f(x, x3, x2, x1, k), count])
        return pt

    else:
        return 'Метод не сходится !'





def f(x, x3, x2, x1, k):
    return x3 * pow(x, 3) + x2 * pow(x, 2) + x1 * x + k


def maxQ(a, b, x3, x2, x1, lambd):
    return max(abs(fiDerivative(a, x3, x2, x1, lambd)), abs(fiDerivative(b, x3, x2, x1, lambd)))


def fDoubleDerivative(x, x3, x2, x1, k):
    return x3 * 6 * x + x2 * 2


def f1(x, x3, x2, x1):
    return 3 * x3 * pow(x, 2) + 2 * x2 * x + x1


def fi(x, x3, x2, x1, k, lambd):
    return x + lambd * f(x, x3, x2, x1, k)


def fiForGraph(x, x3, x2, x1, k, lambd):
    return lambd * f(x, x3, x2, x1, k)


def fiDerivative(x, x3, x2, x1, lambd):
    return 1 + lambd * (3 * x3 * pow(x, 2) + 2 * x2 + x1)


def f2(x, x3, x2):
    return 6 * x3 * x + 2 * x2


def searchX(min_range, max_range, x, x3, x2, x1):
    a = f1(min_range, x3, x2, x1)
    b = f1(max_range, x3, x2, x1)
    c = f1(x, x3, x2, x1)
    if a >= b and a >= c:
        return min_range
    else:
        if b >= a and b >= c:
            return max_range
        else:
            return x


def getLambda(a, b, x3, x2, x1):
    if (f1(a, x3, x2, x1) < f1(b, x3, x2, x1)):
        l = -1 / f1(b, x3, x2, x1)
    else:
        l = -1 / f1(a, x3, x2, x1)
    return l


def printGraph(a, b):
    fig, ax = plt.subplots()
    x = np.linspace(a, b, 100)
    y = x3 * pow(x, 3) + x2 * pow(x, 2) + x1 * x + k
    ax.plot(x, y)
    plt.show()


def printGraphIterazia(a, b):
    fig, ax = plt.subplots()
    x = searchX(a, b, b, x3, x2, x1)
    lambd = getLambda(a, b, x3, x2, x1)
    x = np.linspace(a, b, 100)
    y = fi(x, x3, x2, x1, k, lambd)
    ax.plot(x, y)
    ax.plot(x, x)
    func = f(x, x3, x2, x1, k)
    ax.plot(x, func)
    plt.show()


if __name__ == '__main__':
    answerGiven = True
    scannerline = True

    while scannerline:
        print('Ввод из файла/из строки (1/0): ')
        mes = input()
        if mes == '1':
            try:
                pathh = open('lol', 'r')
                x3, x2, x1, k = map(float, pathh.readline().split(' '))
                a, b = map(float, pathh.readline().split(' '))
                eps = float(pathh.readline())

                answerGiven = True
                scannerline = False
            finally:
                pathh.close()
        else:
            print('Коэффициент перед x^3: ')
            x3 = float(input())
            print('Коэффициент перед x^2: ')
            x2 = float(input())
            print('Коэффициент перед x^1: ')
            x1 = float(input())
            print('Свободный член: ')
            k = float(input())
            print('Левая граница приближения: ')
            a = float(input())
            print('Правая граница приближения: ')
            b = float(input())
            print('Погрешность: ')
            eps = float(input())
            answerGiven = True
            scannerline = False

    print('Выберите метод: \n' +
          '1. Метод хорд \n' +
          '2. Метод секущих \n' +
          '3. Метод простой итерации \n' +
          'Ваш ответ: ')
    answer = ''

    while answerGiven:
        give = input()
        if give == '1':
            printGraph(a, b)
            print(func2(x3, x2, x1, k, a, b, eps))
            answerGiven = False
        elif give == '2':
            printGraph(a, b)
            print(func4(x3, x2, x1, k, a, b, eps))
            answerGiven = False
        elif give == '3':
            printGraphIterazia(a, b)
            print(func5(x3, x2, x1, k, a, b, eps))
            answerGiven = False
        else:
            print('Ошибка: не тот номер \n' +
                  'попробуйте еще раз')
            continue
