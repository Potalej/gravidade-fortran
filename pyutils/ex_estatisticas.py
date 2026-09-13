"""
Exemplo de script para lidar com o estatisticas.nml
Basta informar o nome da pasta (eg 20260913_001) que ele faz o trabalho
e plota.
"""
import f90nml
import json
import matplotlib.pyplot as plt
import numpy as np
import sys
import os

diretorio_padrao = "../out/data/"
cor = "#bb8cd1"

if len(sys.argv) < 2:
  print("Informe um subdiretorio")
  sys.exit(1)

diretorio = sys.argv[1]
diretorio = diretorio if diretorio[-1] == "/" else diretorio + "/"
diretorio = diretorio_padrao + diretorio

os.makedirs(diretorio + "img", exist_ok=True)

with open(f"{diretorio}vi.json", "r") as arq_json:
  infos_json = json.load(arq_json)

t0 = infos_json["integracao"]["t0"]
tf = infos_json["integracao"]["tf"]
checkpoints = infos_json["integracao"]["checkpoints"]
eixo_t = np.linspace(t0, tf, num=checkpoints+1, endpoint=True)

nml = f90nml.read(f"{diretorio}estatisticas.nml")

## energia
energia = nml["energia_nml"]

print(f"E inicial: {energia["E"][0]}")
print(f"E final: {energia["E_final"]}")
print(f"E medio: {energia["E_medio"]}")

fig, axs = plt.subplots(2, 1)

axs[0].set_title(r"Energia total")
axs[0].plot(eixo_t, energia["E"], c=cor)
axs[0].set_ylabel(r"$E(q, p)$")
axs[0].set_xlabel(r"$t$")
axs[0].grid(True)

axs[1].set_title("Erro relativo")
axs[1].plot(eixo_t[1:], energia["E_err_rel"], c=cor)
axs[1].set_ylabel(r"$|E - E_0|/E_0$")
axs[1].set_xlabel(r"$t$")
axs[1].grid(True)

plt.tight_layout()
plt.savefig(f"{diretorio}img/energia.png")
plt.show()

## integrais
integrais = nml["integrais_nml"]

fig, axs = plt.subplots(3, 1)

axs[0].plot(eixo_t, integrais["Jx"], label=r"$J_x$")
axs[0].plot(eixo_t, integrais["Jy"], label=r"$J_y$")
axs[0].plot(eixo_t, integrais["Jz"], label=r"$J_z$")
axs[0].set_ylabel(r"$\vec J$")
axs[0].set_xlabel(r"$t$")
axs[0].grid(True)
axs[0].legend()

axs[1].plot(eixo_t, integrais["Px"], label=r"$P_x$")
axs[1].plot(eixo_t, integrais["Py"], label=r"$P_y$")
axs[1].plot(eixo_t, integrais["Pz"], label=r"$P_z$")
axs[1].set_ylabel(r"$\vec P$")
axs[1].set_xlabel(r"$t$")
axs[1].grid(True)
axs[1].legend()

axs[2].plot(eixo_t, integrais["Qcmx"], label=r"$Q_x$")
axs[2].plot(eixo_t, integrais["Qcmy"], label=r"$Q_y$")
axs[2].plot(eixo_t, integrais["Qcmz"], label=r"$Q_z$")
axs[2].set_ylabel(r"$\vec Q_{cm}$")
axs[2].set_xlabel(r"$t$")
axs[2].grid(True)
axs[2].legend()

plt.tight_layout()
plt.savefig(f"{diretorio}img/integrais.png")
plt.show()

## dinamica
dinamica = nml["dinamica_nml"]

fig, axs = plt.subplots(1, 1)
axs = [axs]

axs[0].set_title("Raio de meia massa")
axs[0].plot(eixo_t, dinamica["rmm"], c=cor)
axs[0].set_ylabel(r"$\vec r_{hm}$")
axs[0].set_xlabel(r"$t$")
axs[0].grid(True)

plt.tight_layout()
plt.savefig(f"{diretorio}img/dinamica.png")
plt.show()