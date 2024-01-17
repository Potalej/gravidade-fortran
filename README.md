# 🪐 gravidade-fortran 🪐

A mesma gravidade, só que em Fortran. Vamo que vamo :)

## ⚙️ Compilando

Para gerar uma build, basta rodar
```
cmake -B build -G Ninja
```
e para compilar com o `ninja`:

```
ninja -C build
```

Um executável será gerado no diretório `build`.

Outra possibilidade é com o uso dos helpers. Use:

```
sh helpers/build.sh
```

para compilar o programa. Se quiser compilar e rodar um exemplo de preset, use:

```
sh helpers/compilar_gerar.sh
```

## 🧮 Simulando com Presets

É possível utilizar presets (i.e., arquivos com configurações pré-definidas) para gerar valores iniciais aleatórios com determinadas condições, como com uma energia total desejada, por exemplo.

Com base em um preset modelo disponível em "presets/", escreva seu preset e rode com:

```
./gravidade preset=\"SEU_ARQUIVO.txt\"
```

Atenção às barras invertidas, que são necessárias para a leitura do nome do arquivo no Fortran.

## 🧮 Simulando

Para determinar as condições iniciais do problema, é necessário escrever um "main.f90" conforme os exemplos na pasta "exemplos". É possível gerar condições iniciais aleatórias para qualquer valor de N, sendo estes condicionados ou não, ou inserir os dados manualmente. O "main.f90" padrão é um problema de 3 corpos cujas trajetórias formam uma lemniscata, com o método de Verlet.

Os métodos de integração disponíveis são:
- [Runge-Kutta de ordem 4 (RK4)](https://pt.wikipedia.org/wiki/M%C3%A9todo_de_Runge-Kutta#O_m%C3%A9todo_Runge%E2%80%93Kutta_cl%C3%A1ssico_de_quarta_ordem)
- [Runge-Kutta-Fehlberg (RKF45)](https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta%E2%80%93Fehlberg_method)
- [Velocity Verlet (Simplético)](https://en.wikipedia.org/wiki/Verlet_integration#Velocity_Verlet)

Também pode ser utilizada uma correção numérica baseada na [condição de 1ª ordem de Karush-Kuhn-Tucker](https://en.wikipedia.org/wiki/Karush%E2%80%93Kuhn%E2%80%93Tucker_conditions). Seu uso, porém, gera um custo computacional 7x maior se aplicado a cada passo. Use com moderação.

As massas, posições e momentos lineares são armazenados em arquivos .csv no diretório "data". Para evitar sobreescrita de dados, o nome do arquivo captura a data corrente no formato "aaaammdd_vv.csv", onde "v" se refere à versão do dia, iniciando em 10 e indo até 99.

A análise dos dados pode ser feita com Python através de [gravidade-analise](https://github.com/Potalej/gravidade-analise).