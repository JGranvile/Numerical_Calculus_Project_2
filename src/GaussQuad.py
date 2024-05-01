# %%
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import scipy
import math

# %%
class BasisPolynomials():
    """
    Uma classe para calcular polinômios básicos e suas derivadas,
    especificamente para polinômios de Legendre e Chebyshev.
    """

    def iterative_legendre(self, n, x):
        """
        Calcule o polinômio de Legendre de grau 'n' no ponto 'x' iterativamente.

        Parâmetros:
        n (int): Grau do polinômio de Legendre.
        x (float): O ponto em que o polinômio é avaliado.

        Retorna:
        float: O valor do polinômio de Legendre de grau 'n' em 'x'.
        """
        P0, P1 = 1, x
        for k in range(2, n + 1):
            Pk = ((2*k - 1)*x*P1 - (k - 1)*P0) / k
            P0, P1 = P1, Pk
        return P0 if n==0 else P1 if n==1 else Pk

    def iterative_legendre_derivative(self, n, x):
        """
        Calcule a derivada do polinômio de Legendre de grau 'n' no ponto 'x' iterativamente.

        Parâmetros:
        n (int): Grau do polinômio de Legendre.
        x (float): O ponto em que a derivada é avaliada.

        Retorna:
        float: A derivada do polinômio de Legendre de grau 'n' em 'x'.
        """
        if n == 0:
            return 0
        else:
            return n * (x * self.iterative_legendre(n, x) - self.iterative_legendre(n - 1, x)) / (x**2 - 1)

    def newton_method(self, n, initial_guess):
        """
        Aplique o método de Newton para encontrar raízes do polinômio de Legendre de grau 'n'.

        Parâmetros:
        n (int): Grau do polinômio de Legendre.
        inicial_guess (float): estimativa inicial para a raiz.

        Retorna:
        float: Uma raiz aproximada do polinômio de Legendre de grau 'n'.
        """
        x = initial_guess
        for _ in range(100):
            Pn = self.iterative_legendre(n, x)
            Pn_prime = self.iterative_legendre_derivative(n, x)
            dx = -Pn / Pn_prime
            x += dx
            if abs(dx) < 1e-10:
                break
        return x

    def legendre(self, n, a, b):
        """
        Calcule os nós e pesos para a quadratura de Gauss-Legendre no intervalo [a, b].

        Parâmetros:
        n (int): Número de nós.
        a (float): Limite inferior do intervalo.
        b (float): Limite superior do intervalo.

        Retorna:
        tupla: Uma tupla contendo duas listas, os nós e seus pesos correspondentes.
        """
        nodes = []
        weights = []

        for i in range(1, n + 1):
            
            initial_guess = math.cos(math.pi * (i - 0.25) / (n + 0.5))

            root = self.newton_method(n, initial_guess)
            
            
            transformed_root = 0.5 * ((b - a) * root + a + b)
            nodes.append(transformed_root)
            
            
            weight = 2 / ((1 - root**2) * (self.iterative_legendre_derivative(n, root)**2))
            adjusted_weight = weight * 0.5 * (b - a)
            weights.append(adjusted_weight)

        return nodes, weights

    def chebyshev(self, n, a, b):
        """
        Calcule os nós e pesos para a quadratura de Gauss-Chebyshev no intervalo [a, b].

        Parâmetros:
        n (int): Número de nós.
        a (float): Limite inferior do intervalo.
        b (float): Limite superior do intervalo.

        Retorna:
        tupla: Uma tupla contendo duas listas, os nós e seus pesos correspondentes.
        """
        nodes = [(0.5 * ((b - a) * math.cos(math.pi * (k + 0.5) / n) + a + b)) for k in range(n)]
        weights = [(math.pi / n) * 0.5 * (b - a)] * n  # Adjusted weights
        return nodes, weights



# %%
class GaussQuadrature(BasisPolynomials):
    """
    Uma classe para realizar Quadratura de Gauss usando polinômios de base
    para integração numérica.

    Atributos:
    func (chamável): A função a ser integrada.
    n (int): Número de nós.
    método (str): O método a ser usado ('legendre' ou 'chebyshev').
    intervalo (tupla): O intervalo durante o qual integrar.
    """
    
    def __init__(self, func, n, method="legendre", interval=(-1, 1)):
        """
        Inicialize o objeto GaussQuadrature.

        Parâmetros:
        func (chamável): A função a ser integrada.
        n (int): Número de nós.
        método (str): O método a ser usado ('legendre' ou 'chebyshev').
        intervalo (tupla): O intervalo durante o qual integrar.
        """
        super().__init__()
        self.func = func
        self.n = n
        self.nodes = None
        self.weights = None
        self.method = method
        self.interval = interval

    def generate_and_save(self, filename = 'nodes_and_weights.csv'):
        """
        Gere nós e pesos e salve-os em um arquivo CSV.

        Parâmetros:
        filename (str): O nome do arquivo para salvar os nós e pesos.

        Retorna:
        GaussQuadrature: Retorna a própria instância para encadear métodos.
        """
        
        nodes, weights = self.nodes_and_weights()
    
        df = pd.DataFrame({'Nodes': nodes, 'Weights': weights})
        df.to_csv(filename, index=False)
        return self

    def load(self, filename = 'nodes_and_weights.csv'):
        """
        Carregue nós e pesos de um arquivo CSV.

        Parâmetros:
        filename (str): O nome do arquivo a partir do qual carregar os nós e pesos.

        Retorna:
        GaussQuadrature: Retorna a própria instância para encadear métodos.
        """

        df = pd.read_csv(filename)
        self.nodes, self.weights = df['Nodes'].values, df['Weights'].values
        return self
    
    def nodes_and_weights(self):
        """
        Gere nós e pesos com base no método e intervalo escolhido.

        Retorna:
        tupla: Uma tupla contendo duas listas, os nós e seus pesos correspondentes.
        """

        if self.method == "legendre":
            return self.legendre(self.n, self.interval[0], self.interval[1])
        elif self.method == "chebyshev":
            return self.chebyshev(self.n, self.interval[0], self.interval[1])
        else:
            raise ValueError("Unknown quadrature type. Please choose 'legendre' or 'chebyshev'.")
    
    def gauss_quadrature(self):
        """
        Execute a quadratura de Gauss para aproximar a integral da função.

        Retorna:
        float: A integral aproximada da função.
        """

        return sum(self.weights * self.func(self.nodes))

def Gaussian_Quad(n, interval, func, method='legendre', filename='nodes_and_weights.csv'):
    """
    Execute a quadratura gaussiana para integração numérica.

    Parâmetros:
    n (int): Número de nós.
    intervalo (tupla): O intervalo durante o qual integrar.
    func (chamável): A função a ser integrada.
    método (str): O método a ser usado ('legendre' ou 'chebyshev').
    filename (str): Nome do arquivo para salvar/carregar nós e pesos.

    Retorna:
    tupla: valor integral aproximado, nós e pesos.
    """
    
    quadrature = GaussQuadrature(func, n, method, interval)
    quadrature.generate_and_save(filename).load(filename)
    return quadrature.gauss_quadrature(), quadrature.nodes, quadrature.weights



