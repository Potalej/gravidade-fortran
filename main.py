import platform
from os import system

def rodar (comando:str)->None:
  system(comando)

def build (buildar:bool=False):
  rodar('cmake -B build -G Ninja && ninja -C build')

def sortear_salvar (buildar:bool=False):
  if buildar: build()
  rodar('gravidade.exe -sv ./presets/condicionar/exemplo.txt')

def sortear_simular (buildar:bool=False):
  if buildar: build()
  rodar('gravidade.exe -s ./presets/condicionar/exemplo.txt')

def lemniscata (buildar:bool=False):
  if buildar: build()
  rodar('gravidade.exe -vi ./presets/valores_iniciais/lemniscata.txt')

def trajetorias (buildar:bool=False):
  if buildar: build()
  arq = input('> ')
  if arq == "":
    rodar('gravidade.exe -e ./presets/data/lemniscata.csv')
  else:
    rodar(f'gravidade.exe -e ./out/data/{arq}/data.csv')

def docker (buildar:bool=False):
  rodar('docker buildar -t local:gravidade .')
  rodar('echo "Para acessar o diretorio do gravidade-fortran, use `cd /src`"')
  rodar('docker run -it -v .:/src local:gravidade')

helpers = {
  1: {'nome': 'build',             'func': build},
  2: {'nome': 'Sortear e salvar',  'func': sortear_salvar},
  3: {'nome': 'Sortear e simular', 'func': sortear_simular},
  4: {'nome': 'Lemniscata (v.i.)', 'func': lemniscata},
  5: {'nome': 'Trajetórias',       'func': trajetorias},
  6: {'nome': 'Iniciar Docker',    'func': docker}
}

def main ():

  # Seleciona a funcao a partir do SO
  os = platform.system()
  if os == "Windows": cmd = lambda arq: system(f"cd helpers/windows && {arq}.bat")
  elif os == "Linux": cmd = system(f"ls {arq}.sh")
  else:
    print("Algo de errado nao esta certo...")
    return

  # Lista as possibilidades
  print("Escolha o que fazer:")
  for i in helpers:
    helper=helpers[i]
    print(f"({i}) {helper['nome']}")
  escolha = input('> ')
  try: 
    buildar = not (len(escolha) == 2 and escolha[0] == '0')
    escolha = int(escolha)
  except:
    print("Escolha inválida!")
    return

  helper_escolhido = helpers[escolha]
  helper_escolhido['func'](buildar)

if __name__ == "__main__":
  main()