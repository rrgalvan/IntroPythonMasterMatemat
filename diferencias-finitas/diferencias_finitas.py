# -*- coding: utf-8 -*-
"""
Módulo de diferencias finitas en Python 3.x
"""

# Escribire los objetos que me interese importar
from numpy import diag, ones, linspace, array

def laplace1d(f, ua, ub, n, a=0, b=1):
    """
    Esta función resuelve el problema de Laplace 1d
        u''  = f     en [a,b]
        u(a) = ua
        u(b) = ub
    mediante el método de las diferencias finitas
    """
    h = (b-a)/(n+1) # Tamaño de la partición
    # 1. Matriz
    A_h = (1./h**2) * (
            2*diag( ones(n) ) 
            - diag( ones(n-1), +1 ) 
            - diag( ones(n-1), -1)
            )
    # 2. Segundo miembro
    f_h = []
    x = linspace(0, 1, num=n+2) # x_0, ..., x_{n-1}
    x_interior = x[1:n+1]
    
    f_h = f(x_interior) # f_h es el array resultante de aplicar f a cada elmento del array x
    f_h[0]  += ua/h**2
    f_h[-1] += ub/h**2
    
    # 3. Resolver sistema
    from numpy.linalg import solve
    u_h = solve(A_h, f_h)

    # Concatenamos la solución con los datos en los extremos del intervalo
    u = array( [ua] + list(u_h) + [ub] )
    
    # Devolvemos la partición x y la solución obtenida
    return x, u

def calor_1d_EulerImp(f, ua, ub, u0, n, k, a=0, b=1):
    """
    Esta función resuelve el problema de Laplace 1d
        u''  = f     en [a,b]
        u(a) = ua
        u(b) = ub
    mediante el método de las diferencias finitas
    """
    h = (b-a)/(n+1) # Tamaño de la partición
    # 1. Matriz
    A_h = (k/h**2) * (
            2*diag( ones(n) ) 
            - diag( ones(n-1), +1 ) 
            - diag( ones(n-1), -1)
            )
    # 2. Segundo miembro
    f_h = []
    x = linspace(0, 1, num=n+2) # x_0, ..., x_{n-1}
    x_interior = x[1:n+1]
    
    f_h = f(x_interior) # f_h es el array resultante de aplicar f a cada elmento del array x
    f_h[0]  += ua/h**2
    f_h[-1] += ub/h**2
    
    # 3. Resolver sistema
    from numpy.linalg import solve
    u_h = solve(A_h, f_h)

    # Concatenamos la solución con los datos en los extremos del intervalo
    u = array( [ua] + list(u_h) + [ub] )
    
    # Devolvemos la partición x y la solución obtenida
    return x, u