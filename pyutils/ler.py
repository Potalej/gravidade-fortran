import pandas as pd
import numpy as np
import json

def parse_coordenadas(x, N):
  """
  Dado um vetor x de dimensao 6N, retorna os vetores q e p
  contidos em x
  """
  R = [
    [float(a) for a in x[0:N]],
    [float(a) for a in x[N:2*N]],
    [float(a) for a in x[2*N:3*N]]
  ]
  P = [
    [float(a) for a in x[3*N:4*N]],
    [float(a) for a in x[4*N:5*N]],
    [float(a) for a in x[5*N:6*N]]
  ]
  R = np.array(list(zip(*R)))
  P = np.array(list(zip(*P)))
  return R, P

def ler_simulacao_csv (diretorio:str, arquivo:str='data.csv', chunksize=100, saltos=1):
  """
  Le um arquivo de dados de simulacao (trajetorias) em .csv
  """
  valores_R, valores_P = [], []
  for chunk in pd.read_csv(f"{diretorio}/{arquivo}", chunksize=chunksize, skiprows=3):
    # Aplica a funcao
    def func_res (x):
      [R, P] = [*parse_coordenadas(x, int(len(x)/6))]
      valores_R.append(R)
      valores_P.append(P)
    chunk.apply(func_res, axis=1)
  return valores_R[::saltos], valores_P[::saltos]

def ler_simulacao_bin (diretorio:str, arquivo:str='data.bin'):
  """
  Le um arquivo de dados de simulacao (trajetorias) em .bin
  """
  arq = open(f"{diretorio}/{arquivo}")

  # A primeira linha tem o formato: h, G, N
  tipo_1a_linha = np.dtype([
    ('col1', 'f8'), # h
    ('col2', 'f8'), # G
    ('col3', 'f8'), # eps
    ('col4', 'i'),  # N
    ('col5', 'i'),  # checkpoints
  ])
  primeira_linha = np.fromfile(file=arq, count=1, dtype=tipo_1a_linha)[0]

  # Agora captura as massas (2a linha)
  N = primeira_linha[3]
  tipo_2a_linha = np.dtype([(f'col{i}', 'f8') for i in range(N)])
  massas = np.fromfile(file=arq, count=1, dtype=tipo_2a_linha)[0]

  # Agora le o arquivo inteiro
  tipo_dado = np.dtype([(f'col{i}', 'f8') for i in range(6*N)])
  dados = np.fromfile(file=arq, dtype=tipo_dado)

  dados = [list(x) for x in dados]
  valores_R, valores_P = [], []
  # counter = 0
  for linha in dados:
    [R, P] = [*parse_coordenadas(linha, N)]
    valores_R.append(R)
    valores_P.append(P)
    # counter += 1
    # if counter == 20: break
  return valores_R, valores_P

def ler_valores_iniciais (diretorio:str, arquivo:str='vi.json'):
  with open(diretorio + '/' + arquivo, 'r') as arq:
    return json.load(arq)