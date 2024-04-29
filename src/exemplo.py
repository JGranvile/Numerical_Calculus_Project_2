# Import GaussQuadrature class or Gaussian_Quad function
from GaussQuad import GaussQuadrature, Gaussian_Quad

# Create a function to integrate
def func(x):
    return x**2  # Example function to integrate

# Use Gaussian_Quad to perform the integration
num_nodes = 5
interval = (0, 1)  # Example interval
method = 'legendre'  # or 'chebyshev'

# Using Gaussian_Quad function
result = Gaussian_Quad(num_nodes, interval, func, method)
print("Integration result:", result)

