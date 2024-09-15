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

Um executável será gerado no diretório raiz.

Outra possibilidade é com o uso dos helpers. Use:

```
sh helpers/build.sh
```

para compilar o programa. Se quiser compilar e rodar um exemplo de preset, use:

```
sh helpers/compilar_gerar.sh
```

## 🧮 Simulando

### Presets para geração
É possível utilizar presets (i.e., arquivos com configurações pré-definidas) para gerar valores iniciais aleatórios com determinadas condições, como com uma energia total desejada, por exemplo.

Com base em um preset modelo disponível em "presets/", escreva seu preset e rode com:

```
./gravidade -s SEU_ARQUIVO.txt
```

Há dois modos para geração aleatória no momento:
- `modo sorteio`: Gera valores aleatórios no intervalo desejado e os condiciona conforme as integrais primeiras;
- `modo sorteio_henon`: Gera valores aleatórios no intervalo desejado, mas aplica as condições iniciais de Henón, ignorando as integrais primeiras informadas. Para mais informações, veja a [wiki](https://github.com/Potalej/gravidade-fortran/wiki/Rodando).

### Valores iniciais
Se já tiver os valores iniciais do problema e quiser utilizá-los, é possível através da opção `-vi`. Um modelo de valores iniciais de um problema de três corpos com trajetória em formato de lemniscata está disponível no diretório PRESETs.

Para rodar este caso, utilize:

```
./gravidade -vi SEU_ARQUIVO.txt
```

### Visualizando
Além disso, é possível também visualizar as simulações feitas através da opção `-e` ou `--exibir`. É uma visualização simples das trajetórias, mas pode ajudar a conferir resultados rapidamente.

```
./gravidade -e SEU_ARQUIVO.csv
```

## Métodos
Os métodos de integração disponíveis são:
- Métodos de Runge-Kutta:
    - [Runge-Kutta de ordem 4 (RK4)](https://pt.wikipedia.org/wiki/M%C3%A9todo_de_Runge-Kutta#O_m%C3%A9todo_Runge%E2%80%93Kutta_cl%C3%A1ssico_de_quarta_ordem);
    - [Runge-Kutta-Fehlberg (RKF45)](https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta%E2%80%93Fehlberg_method);
- Métodos simpléticos:
    - [Euler simplético](https://en.wikipedia.org/wiki/Symplectic_integrator#A_first-order_example);
    - [Velocity-Verlet](https://en.wikipedia.org/wiki/Verlet_integration#Velocity_Verlet);
    - [Ruth de 3ª Ordem (RUTH3)](https://en.wikipedia.org/wiki/Symplectic_integrator#A_third-order_example);
    - [Ruth de 4ª Ordem (RUTH4)](https://en.wikipedia.org/wiki/Symplectic_integrator#A_fourth-order_example).

Também pode ser utilizada uma correção numérica baseada na [condição de 1ª ordem de Karush-Kuhn-Tucker](https://en.wikipedia.org/wiki/Karush%E2%80%93Kuhn%E2%80%93Tucker_conditions). Seu uso, porém, gera um custo computacional até 7x maior se aplicado a cada passo. Use com moderação.

As massas, posições e momentos lineares são armazenados em arquivos .csv no diretório "data". Para evitar sobreescrita de dados, o nome do arquivo captura a data corrente no formato "aaaammdd_vv.csv", onde "v" se refere à versão do dia, iniciando em 10 e indo até 99.

A análise dos dados pode ser feita com Python através de [gravidade-analise](https://github.com/Potalej/gravidade-analise).