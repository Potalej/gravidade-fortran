import socket
import numpy as np
import pyqtgraph as pg
from pyqtgraph.Qt import QtWidgets
import sys

app = QtWidgets.QApplication(sys.argv)
win = pg.GraphicsLayoutWidget()
plot = win.addPlot()

# Particulas
sc = plot.plot([], [], pen=None, symbol='o', symbolSize=5)
plot.setXRange(-10,10)
plot.setYRange(-10,10)
plot.enableAutoRange(x=False, y=False)
plot.setAspectLocked(True)

# Texto
texto = pg.TextItem(text="t=0", color="w", anchor=(0, 1))
plot.addItem(texto)
texto.setPos(-9, 9)  # posição em coordenadas do gráfico

# Exibicao
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
    i = 0
    while True:
        if i < 5:
            i = i + 1
            continue
        i = 0 
        try: data = conn.recv(tamanho).decode()
        except: break
        
        try:
            valores = data.strip().splitlines()
            t = float(valores[0])
            pontos = [tuple([float(x) for x in linha.split()]) for linha in valores[1:]]
            x,y,z = zip(*pontos)
        except: continue
        lista_dados.append(pontos)

        texto.setText(f"t={t}")

        # z = np.array(z)
        # z_norm = (z - z.min()) / (z.max() - z.min())
        # cmap = pg.colormap.get("viridis")
        # cores = cmap.map(z_norm, mode='qcolor')

        # sc.setData(x,y, symbolBrush=cores)
        sc.setData(x,y)
        QtWidgets.QApplication.processEvents()

# Fechando a conexao
conn.close()