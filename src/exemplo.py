
# Importe a classe GaussQuadrature ou a função Gaussian_Quad
from GaussQuad import GaussQuadrature, Gaussian_Quad

# Crie uma função para integrar
def func(x):
    return x*10;  # Função de exemplo para integrar

# Use o Gaussian_Quad para realizar a integração
num_nodes = 5
intervalo = (574, 1314)  # Intervalo de exemplo
método = 'legendre'  # ou 'chebyshev'

# Usando a função Gaussian_Quad
resultado = Gaussian_Quad(num_nodes, intervalo, func, método)
print("Resultado da integração:", resultado)

"""
    Realiza quadratura gaussiana para integração numérica.

    Parâmetros:
    n (int): Número de nós.
    intervalo (tuple): O intervalo sobre o qual integrar.
    func (callable): A função a ser integrada.
    método (str): O método a ser usado ('legendre' ou 'chebyshev').
    filename (str): Nome do arquivo para salvar/carregar nós e pesos.

    Retorna:
    tuple: Valor aproximado da integral, nós e pesos.
    """
