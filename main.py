import platform
from os import system

def rodar (comando:str)->None:
  system(comando)

def build (buildar:bool=False, os:str="Linux"):
  rodar('cmake -B build -G Ninja && ninja -C build')

def sortear_salvar (buildar:bool=False, os:str="Linux"):
  if buildar: build()
  if os == "Windows":
    rodar('gravidade.exe -sv ./presets/condicionar/exemplo.json')
  elif os == "Linux":
    rodar('./gravidade -sv ./presets/condicionar/exemplo.txt')

def sortear_simular (buildar:bool=False, os:str="Linux"):
  if buildar: build()
  if os == "Windows":
    rodar('gravidade.exe -s ./presets/condicionar/exemplo.json')
  elif os == "Linux":
    rodar('./gravidade -s ./presets/condicionar/exemplo.txt')

def lemniscata (buildar:bool=False, os:str="Linux"):
  if buildar: build()
  if os == "Windows":
    rodar('gravidade.exe -vi ./presets/valores_iniciais/exemplo_vi.json')
  elif os == "Linux":
    rodar('./gravidade -vi ./presets/valores_iniciais/lemniscata.txt')

def trajetorias (buildar:bool=False, os:str="Linux"):
  if buildar: build()
  arq = input('> ')
  if arq == "":
    if os == "Windows":
      rodar('gravidade.exe -e ./presets/data/lemniscata.csv')
    elif os == "Linux":
      rodar('./gravidade -e ./presets/data/lemniscata.csv')
  else:
    if os == "Windows":
      rodar(f'gravidade.exe -e ./out/data/{arq}/data.csv')
    elif os == "Linux":
      rodar(f'./gravidade -e ./out/data/{arq}/data.csv')

def docker (buildar:bool=False, os:str="Linux"):
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
  elif os == "Linux": cmd = lambda arq: system(f"ls {arq}.sh")
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
  helper_escolhido['func'](buildar, os)

if __name__ == "__main__":
  main()