import numpy as np
import matplotlib.pyplot as plt
import threading as tr
import sympy as sp
import math
import random
from matplotlib.animation import FuncAnimation
from mpl_toolkits import mplot3d
from tqdm import tqdm

a_0 = .529e-10

radius, pheda = sp.symbols("radius pheda")

class HydrogenAtom:
    def __init__(self, n, m, l, debugger = True):
        super(HydrogenAtom, self).__init__()
        self.n, self.m, self.l = n, m, l
        print(r"[+] Hydrogen atom model created!"*debugger, sep='')
        self.debugger = debugger
        self.function: int 
        self.supportsObjects = False
            
    def nthDerivative(self, n, func):
        if n <= 1: return sp.diff(func)
        return self.nthDerivative(n-1, sp.diff(func))

    # ritorna n!
    def factorial(self, n):
        if n <= 1: return 1
        return n*self.factorial(n-1)

    def sphericalArmonic(self, m, l):
        negative = -1 ** ( ( abs(m) + m ) / 2 )
        first_fact = np.sqrt( ( (2*l +1) / 4*np.pi) * self.factorial (l - abs(m)) / self.factorial (l + abs(m)) ) 
        res = first_fact * negative * sp.assoc_legendre(abs(m), l, sp.cos(pheda))
        del negative, first_fact
        
        return res

    #ritorna il fattore di normalizzazione 
    def normalizationFactor(self, n ,l):
        first_term = (2/(a_0*n+1e-7))**3
        second_term = self.factorial(n - l -1)/( 2*n * (self.factorial(n+l))**3 +1e-7)
        return np.sqrt(first_term*second_term)


    def RadialFunction(self, n,l):
        p = 2*radius / (n*a_0+1e-7)
        return -self.normalizationFactor(n,l) * (p)**l * np.e ** (-p/2) * sp.assoc_laguerre(n-l+1, 2*l+1, 2*p)

    def psi(self, angle, r):
        return abs(self.function.subs(pheda, angle).subs(radius, r).evalf())**2

    def plot(self, resolution, size = 2e-3, simType= '2d') -> None:
        if simType != '2d': return self.plot3d(resolution, size)
        print(r'[+] Solving Psi...'*self.debugger, sep="")

        self.function = self.RadialFunction(self.n, self.l) * self.sphericalArmonic(self.m, self.l)

        print(r'[+] Generating the heatmap...'*self.debugger, sep='')
        linspc = np.linspace(-size, size, resolution)
        
        if self.debugger: 
            print()
            with tqdm(total=resolution**2) as pbar:
                graph = np.array([np.array([float(self.psi(math.atan2(x, y), x**2+y**2))+bool(pbar.update(1)) for x in linspc]) for y in linspc])
            print()
        else:graph = np.array([np.array([float(self.psi(math.atan2(x, y), x**2+y**2)) for x in linspc]) for y in linspc])

        fig, ax = plt.subplots()

        heatmap = ax.imshow(graph, cmap='inferno')
        plt.show()

    def plot3d(self, resolution, size = 2e-3):
        self.function = self.RadialFunction(self.n, self.l) * self.sphericalArmonic(self.m, self.l)

        linspc = np.linspace(-size, size, resolution)
        X, Y = np.meshgrid(linspc, linspc)
        graph = np.array([np.array([float(self.psi(math.atan2(x, y), x**2+y**2)) for x in linspc]) for y in linspc])

        fig = plt.figure()
        ax = fig.add_subplot(111, projection='3d')

        heatmap = ax.plot_surface(X, Y, graph, cmap='inferno')
        plt.show() 
