def calcular_Q():
    x = float(input("Digite o valor de x: "))
    y = float(input("Digite o valor de y: "))
    Q = 199870000 * (x/1.225)**0.5 * (y/7.9053)**3.15
    return Q

# Chamando a função e exibindo o resultado
resultado = calcular_Q()
print("O valor de Q é:", resultado)