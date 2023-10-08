import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation
from mpl_toolkits import mplot3d
import sympy as sp
import math
import random

a_0 = .529e-10

radius, pheda = sp.symbols("radius pheda")

def nthDerivative(n, func):
    if n <= 1: return sp.diff(func)
    return nthDerivative(n-1, sp.diff(func))

# ritorna n!
def factorial(n):
    if n <= 1: return 1
    return n*factorial(n-1)

# ritorna il polinomio di laguerre per r s e x

def sphericalArmonic(m, l):
    negative = -1 ** ( ( abs(m) + m ) / 2 )
    first_fact = np.sqrt( ( (2*l +1) / 4*np.pi) * factorial (l - abs(m)) / factorial (l + abs(m)) ) 
    res = first_fact * negative * sp.assoc_legendre(abs(m), l, sp.cos(pheda))
    del negative, first_fact
    
    return res


#ritorna il fattore di normalizzazione 
def normalizationFactor(n ,l):
    first_term = (2/(a_0*n+1e-7))**3
    second_term = factorial(n - l -1)/( 2*n * (factorial(n+l))**3 +1e-7)
    return np.sqrt(first_term*second_term)


def RadialFunction(n,l):
    p = 2*radius / (n*a_0+1e-7)
    return -normalizationFactor(n,l) * (p)**l * np.e ** (-p/2) * sp.assoc_laguerre(n-l+1, 2*l+1, 2*p)

function: int

def wavefunction(angle, r):
    return function.subs(pheda, angle).subs(radius, r).evalf()

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#+++++++++++++++++++++++++++++++++++++++ Simulazione ++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# 1.5s

resolution = 150
size = 1.5e-3


n,m,l = 4, 3, 0

function = RadialFunction(n, l) * sphericalArmonic(m, l)
function = sp.simplify(function)

linspc = np.linspace(-size, size, resolution)
graph = np.array([np.array([float(wavefunction(math.atan2(x, y), x**2+y**2)) for x in linspc]) for y in linspc])

        


fig, ax = plt.subplots()

heatmap = ax.imshow(graph, cmap = "inferno")

plt.show()
