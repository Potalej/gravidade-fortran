import socket
import numpy as np
import pyqtgraph as pg
from pyqtgraph.Qt import QtWidgets
import sys

app = QtWidgets.QApplication(sys.argv)
win = pg.GraphicsLayoutWidget()
plot = win.addPlot()
sc = plot.plot([], [], pen=None, symbol='o', symbolSize=5)
plot.setXRange(-10,10)
plot.setYRange(-10,10)
plot.enableAutoRange(x=False, y=False)
plot.setAspectLocked(True)

win.show()

# Criando socket e ligando servidor
s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
s.bind(('127.0.0.1', 50007))
s.listen()

# Aceita conexao
conn, addr = s.accept()

lista_dados = []
with conn:
    print(f"Conectado por {addr}")

    data = conn.recv(4)
    N = int.from_bytes(data, byteorder='little', signed=False)
    print(f"N = {N}")
    tamanho = N * 60

    # Agora começa a receber os dados em tempo real
    while True:

        try: data = conn.recv(tamanho).decode()
        except: break
        
        try: 
            pontos = [tuple([float(x) for x in linha.split()]) for linha in data.strip().splitlines()]
            x,y,z = zip(*pontos)
        except: continue
        lista_dados.append(pontos)


        sc.setData(x,y)
        QtWidgets.QApplication.processEvents()

# Fechando a conexao
conn.close()