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
./gravidade -s SEU_ARQUIVO.json
```

Há três modos para geração aleatória no momento, que geram valores iniciais aleatórios e os condicionam:
- `sorteio_ip_iterativo`: Condicionamento iterativo;
- `sorteio_ip_direto`: Condicionamento direto;
- `sorteio_aarseth`: Condiciona diretamente e depois aplica o proposto por (Aarseth, 2003) para obter o equilíbrio de virial.

Confira [aqui](/presets/condicionar/exemplo.json) um exemplo de arquivo de valores para sorteio.

### Valores iniciais
Se já tiver os valores iniciais do problema e quiser utilizá-los, é possível através da opção `-vi`. Um modelo de valores iniciais de um problema de três corpos com trajetória em formato de lemniscata está disponível no diretório PRESETs.

Para rodar este caso, utilize:

```
./gravidade -vi SEU_ARQUIVO.json
```

Confira [aqui](/presets/valores_iniciais/exemplo_vi.json) um exemplo de arquivo de valores iniciais.

### Visualizando
Além disso, é possível também visualizar as simulações feitas através da opção `-e` ou `--exibir`. É uma visualização simples das trajetórias, mas pode ajudar a conferir resultados rapidamente.

```
./gravidade -e SEU_ARQUIVO.csv
```

## Métodos
Os métodos de integração implementados são:

- Métodos de Runge-Kutta:
    - [Runge-Kutta de ordem 2 (RK2)](https://pt.wikipedia.org/wiki/M%C3%A9todo_de_Runge-Kutta);
    - [Runge-Kutta de ordem 3 (RK3)](https://pt.wikipedia.org/wiki/M%C3%A9todo_de_Runge-Kutta);
    - [Runge-Kutta de ordem 4 (RK4)](https://pt.wikipedia.org/wiki/M%C3%A9todo_de_Runge-Kutta#O_m%C3%A9todo_Runge%E2%80%93Kutta_cl%C3%A1ssico_de_quarta_ordem);
    - [Runge-Kutta-Fehlberg (RKF45)](https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta%E2%80%93Fehlberg_method) (INDISPONÍVEL);
- Métodos simpléticos:
    - [Euler simplético](https://en.wikipedia.org/wiki/Symplectic_integrator#A_first-order_example);
    - [Velocity-Verlet](https://en.wikipedia.org/wiki/Verlet_integration#Velocity_Verlet);
    - [Ruth de 3ª Ordem (RUTH3)](https://en.wikipedia.org/wiki/Symplectic_integrator#A_third-order_example);
    - [Ruth de 4ª Ordem (RUTH4)](https://en.wikipedia.org/wiki/Symplectic_integrator#A_fourth-order_example).
    - Stormer-Verlet Composto de 8ª Ordem e 15 Estágios (svcp8s15);.
    - Stormer-Verlet Composto de 10ª Ordem e 35 Estágios (svcp10s35).

Também há três aplicações para modificar soluções:

1. Um corretor numérico baseado na [condição de 1ª ordem de Karush-Kuhn-Tucker](https://en.wikipedia.org/wiki/Karush%E2%80%93Kuhn%E2%80%93Tucker_conditions). Bastante útil para simulações com muitos corpos;
2. Colisões perfeitamente elásticas entre os corpos.
3. Um amortecedor no potencial, que impede aproximações muito intensas. Para poucos corpos pode gerar instabilidades, mas é útil para grandes quantidades de corpos.

As massas, posições e momentos lineares são armazenados em arquivos .csv no diretório "data". Para evitar sobreescrita de dados, o nome do arquivo captura a data corrente no formato "aaaammdd_vv.csv", onde "v" se refere à versão do dia, iniciando em 001 e indo até 999.

A análise dos dados pode ser feita com Python através de [gravidade-analise](https://github.com/Potalej/gravidade-analise).


## Referências

* AARSETH, Sverre. Gravitational N-Body Simulations: Tools and Algorithms. Cambridge: Cambridge University Press, 2003.
* VOLCHAN, Sérgio. Uma Introdução à Mecânica Celeste. Rio de Janeiro: Instituto Nacional de Matemática Pura e Aplicada - IMPA, 2007.
* BERTSEKAS, Dmitri Panteli. Nonlinear Programming. 3ed. Nashua: Athena Scientific, 2016.
* HAIRER, Ernst; WANNER, Gerhard; LUBICH, Christian. Geometric Numerical Integration: Structure-Preserving Algorithms for Ordinary Differential Equations. Heidelberg: Springer-Verlag, 2006. DOI: 10.1007/3-540-30666-8. Disponível em: https://doi.org/10.1007/3-540-30666-8.
* ROMA, Alexandre et al. Métodos para a solução numérica de equações diferenciais ordinárias a valores iniciais. São Paulo: Notas de aula, 2019.
* OKUNBOR, D. I.; SKEEL, R. D. Canonical Runge—Kutta—Nyström methods of orders five and six. Journal of Computational and Applied Mathematics, v. 51, n. 3, p. 375–382, jun. 1994. 