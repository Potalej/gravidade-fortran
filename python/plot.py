from numpy import genfromtxt
import matplotlib.pyplot as plt
from sys import argv
from os import listdir

diretorio = argv[1] # Diretorio passado como parametro
diretorio_completo = f'./plot/{diretorio}'
arquivos_corpos = listdir(diretorio_completo) # Captura os arquivos no diretorio

for arquivo in arquivos_corpos:
  with open(f'{diretorio_completo}/{arquivo}', 'r') as arq:
    linhas = [
      [float(valor) for valor in linha.split()]
      for linha in arq.read().split("\n")
      if linha.strip() != ""
    ]
    x,y = list(zip(*linhas))
    plt.plot(x,y)
plt.show() # Exibe    
    