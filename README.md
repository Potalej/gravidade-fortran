<h1 align="center">
    <img src="img/gf-logo.png" width="500">
</h1>

<p align="center">
    <img src="https://img.shields.io/github/languages/top/potalej/gravidade-fortran?label=Fortran&logo=fortran&labelColor=%236f4c91&color=gray" alt="fortran">
    <img src="https://img.shields.io/github/languages/code-size/potalej/gravidade-fortran?label=tamanho&labelColor=%236f4c91&color=gray" alt="tamanho">
    <img src="https://img.shields.io/github/last-commit/potalej/gravidade-fortran?label=alterado&labelColor=%236f4c91&color=gray" alt="ultimo commit">
    <img src="https://img.shields.io/github/issues/potalej/gravidade-fortran?label=issues&labelColor=%236f4c91&color=gray" alt="issues">
</p>

<p align="center">
    <img src="https://img.shields.io/badge/JSON--Fortran-v9.2.0-dark?logo=fortran&color=blue" alt="json-fortran">
    <img src="https://img.shields.io/badge/ncorpos--utilidades-v1.0.1-brightgreen?logo=fortran&color=blue" alt="ncorpos-utilidades">
    <img src="https://img.shields.io/badge/ncorpos--valores--iniciais-v0.2.5-dark?logo=fortran&color=blue" alt="ncorpos-valores-iniciais">
</p>

<p align="center">A mesma gravidade, só que em Fortran. Vamo que vamo :)</p>

> *"Meus N-corpos, minhas regras"*   
> — ros

## ⚙️ Compilando

Depois de clonar o repositório, rode

```
git submodule update --init --recursive
```

para as bibliotecas [ncorpos-valores-iniciais](https://github.com/potalej/ncorpos-valores-iniciais) e [ncorpos-utilidades](https://github.com/potalej/ncorpos-utilidades) ficarem disponíveis.

Para compilar:

```
cmake -B build
cd build 
make
cd ..
```

ou usando, por exemplo, o Ninja:

```
cmake -B build -G Ninja
ninja -C build
```

Um executável será copiado para o diretório raiz.

### 🚩 Flags

#### `PRECISAO=64`

Pode ser 32, 64 (padrão) ou 128. Todas as variáveis do tipo `REAL(pf)` terão a precisão definida. Exemplo de uso:

```
cmake -B build -DPRECISAO=128
```

> [!WARNING]
> Ao utilizar uma precisão diferente de 64, a opção `FETCHJSONFORTRAN` abaixo será desconsiderada e o JSON-Fortran será compilado localmente com a precisão correspondente.

<h4 id="fetchjsonfortran">

`FETCHJSONFORTRAN=OFF`

</h4>

Se ativada, baixa uma versão recente da API [JSON-Fortran](https://github.com/jacobwilliams/json-fortran) e a compila localmente; isso pode demorar. Por padrão, vem desativada e tenta usar uma versão disponibilizada via gerenciador de pacotes (e.g.: conda). Exemplo de uso:

```
cmake -B build -DFETCHJSONFORTRAN=ON
```

#### `USAR_GPU=OFF` (experimental)

Ativa a paralelização em GPU utilizando o OpenMP offload. Foi implementado para experimentar e aprender, não utilize.

#### `GPROF=OFF`

Ativa o [GNU Profiler](https://ftp.gnu.org/old-gnu/Manuals/gprof-2.9.1/html_mono/gprof.html), utilizado para análise de desempenho. Para ativar, basta usar:

```
cmake -B build -DGPROF=ON ...
```

Após rodar uma simulação com o programa compilado com gprof, será gerado um arquivo "gmon.out" que fornece um relatório usando o comando:

```
gprof gravidade gmon.out > relatorio.txt
```

## 🧮 Simulando

### Configurações

Para todos os casos de simulação que seguem, os seguintes parâmetros podem ser passados para modificar o comportamento do programa:

#### `-ps, --pasta-saida`

Por padrão, a pasta de saída do programa é "out". Caso queira utilizar outra pasta, basta informar o nome:

```bash
./gravidade [...] -ps minha_pasta
```

#### `-es, --extensao-saida`

O programa suporta dois tipos de saída: ".csv" e ".bin" (padrão). Para escolher, basta informar:

```bash
./gravidade [...] -es .csv
```

### Presets para geração

É possível utilizar presets (i.e., arquivos com configurações pré-definidas) para gerar valores iniciais aleatórios com determinadas condições, como com uma energia total desejada, por exemplo.

Com base nos presets disponíveis em "exemplos/", escreva seu preset e rode com:

```
./gravidade -s SEU_ARQUIVO.json
```

Os seguintes modos de condicionamento de valores iniciais estão disponíveis:
- `sorteio_ip_iterativo`: Condicionamento iterativo. Suporta potencial amortecido.
- `sorteio_ip_direto`: Condicionamento direto. Não suporta potencial amortecido.
- `sorteio_aarseth`: Condicionamento direto para $E=-1/4$ e outras integrais nulas, começando em equilíbrio. Não suporta potencial amortecido. Proposto por (Aarseth, 2003).
- `sorteio_aarseth_modificado`: Condicionamento direto ou iterativo (a depender) para $E=-1/4$ e outas integrais nulas, começando em equilíbrio. Suporta potencial amortecido, utilizando o condicionamento iterativo se for o caso.

Confira [aqui](/exemplos/sortear/exemplo.json) um exemplo de arquivo de valores para sorteio.

### Valores iniciais

Se já tiver os valores iniciais do problema e quiser utilizá-los, é possível através do modo `-vi`. Um modelo de valores iniciais de um problema de três corpos com trajetória em formato de lemniscata está disponível nos exemplos.

Para rodar este caso, utilize:

```
./gravidade -vi SEU_ARQUIVO.json
```

Confira [aqui](/exemplos/valores_iniciais/exemplo.json) um exemplo de arquivo de valores iniciais.

### Testes

O programa contém testes que podem ser rodados via:

```
./gravidade -t
```

Os testes rodarão e exibirão uma saída que indicará se o programa está funcionando corretamente ou não.

### Visualizando 🧐
É possível visualizar as simulações de duas formas.

#### Gnuplot ou Matplotlib
Através da opção `-e` ou `--exibir`. É uma visualização simples das trajetórias, mas pode ajudar a conferir resultados rapidamente.

```
./gravidade -e SEU_ARQUIVO.csv
```

#### Visualização em tempo real
Também podemos visualizar as simulações em tempo real via sockets, ativando a opção `"exibir": true` no arquivo de valores iniciais utilizado. Para isso, é preciso ter um servidor local aberto para receber as informações, que são enviadas via TCP.

Um exemplo de script em Python que recebe e exibe os dados em tempo real pode ser encontrado em [aqui](/pyutils/ex_servidor.py). Basta rodar este script e em seguida rodar sua simulação com a opção ativada. É importante observar que o envio de dados prejudica consideravelmente o desempenho, então não é bom usar isso para simulações grandes e sérias.

![Exemplo de visualização em tempo real via socket](/img/exemplo_visualizacao.gif)

## 📈 Métodos
Os métodos de integração implementados são:

- Métodos gerais:
    - Euler explícito e implícito;  
    - Métodos de Runge-Kutta:
        - [Runge-Kutta de ordem 2 (RK2)](https://pt.wikipedia.org/wiki/M%C3%A9todo_de_Runge-Kutta);
        - [Runge-Kutta de ordem 3 (RK3)](https://pt.wikipedia.org/wiki/M%C3%A9todo_de_Runge-Kutta);
        - [Runge-Kutta de ordem 4 (RK4)](https://pt.wikipedia.org/wiki/M%C3%A9todo_de_Runge-Kutta#O_m%C3%A9todo_Runge%E2%80%93Kutta_cl%C3%A1ssico_de_quarta_ordem);
        - [Runge-Kutta-Fehlberg (RKF45)](https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta%E2%80%93Fehlberg_method) (DESATIVADO TEMPORARIAMENTE);
- Métodos simpléticos:
    - [Euler simplético](https://en.wikipedia.org/wiki/Symplectic_integrator#A_first-order_example);
    - [Velocity-Verlet](https://en.wikipedia.org/wiki/Verlet_integration#Velocity_Verlet);
    - [Ruth de 3ª Ordem (RUTH3)](https://en.wikipedia.org/wiki/Symplectic_integrator#A_third-order_example);
    - [Ruth de 4ª Ordem (RUTH4)](https://en.wikipedia.org/wiki/Symplectic_integrator#A_fourth-order_example).
    - Runge-Kutta-Nystrom de 5ª Ordem e 5 passos (rkn551);
    - Runge-Kutta-Nystrom de 6ª Ordem e 7 passos (rkn671);
    - Stormer-Verlet Composto de 8ª Ordem e 15 Estágios (svcp8s15);.
    - Stormer-Verlet Composto de 10ª Ordem e 35 Estágios (svcp10s35).

É possível calcular as forças de duas formas:
- Diretamente, sem parâmetros (custo $O(N^2)$);
- Barnes-Hut, com parâmetro $\theta > 0$ (custo $O(N \log N)$ para $N$ grande).

Também há três aplicações para modificar soluções:

1. Um corretor numérico que aplica correções a partir dos desvios de energia total (e outras integrais primeiras). Bastante útil para simulações com muitos corpos;
2. Colisões perfeitamente elásticas entre os corpos.
3. Um amortecedor no potencial, que impede aproximações muito intensas. Para poucos corpos pode gerar instabilidades, mas é útil para grandes quantidades de corpos.

As massas, posições e momentos lineares são armazenados em arquivos .bin (padrão) ou .csv no diretório "data". Para evitar sobreescrita de dados, o nome do arquivo captura a data corrente no formato "aaaammdd_vvv.*", onde "v" se refere à versão do dia, iniciando em 001 e indo até 999.

## 📚 Dependências

- [OpenBLAS](https://github.com/OpenMathLib/OpenBLAS): Rotinas numéricas de álgebra linear. Instale via algum gerenciador de pacotes, como conda ou pacman. Os diretórios onde OpenBLAS é procurado [estão aqui](https://github.com/Potalej/gravidade-fortran/blob/main/cmake/FindOpenBLAS.cmake); se der erro, provavelmente o erro está aqui.
- [JSON-Fortran](https://github.com/jacobwilliams/json-fortran/): Rotinas para manipulação de arquivos JSON via Fortran. O CMake por padrão procurará a biblioteca no sistema, então instale via pacotes. Se não for possível, ative a flag [`-DFETCHJSONFORTRAN=ON`](#fetchjsonfortran) que uma versão recente do JSON-Fortran será compilada localmente junto do programa.

## Referências

* AARSETH, Sverre. Gravitational N-Body Simulations: Tools and Algorithms. Cambridge: Cambridge University Press, 2003.
* VOLCHAN, Sérgio. Uma Introdução à Mecânica Celeste. Rio de Janeiro: Instituto Nacional de Matemática Pura e Aplicada - IMPA, 2007.
* BERTSEKAS, Dmitri Panteli. Nonlinear Programming. 3ed. Nashua: Athena Scientific, 2016.
* HAIRER, Ernst; WANNER, Gerhard; LUBICH, Christian. Geometric Numerical Integration: Structure-Preserving Algorithms for Ordinary Differential Equations. Heidelberg: Springer-Verlag, 2006. DOI: 10.1007/3-540-30666-8. Disponível em: https://doi.org/10.1007/3-540-30666-8.
* ROMA, Alexandre et al. Métodos para a solução numérica de equações diferenciais ordinárias a valores iniciais. São Paulo: Notas de aula, 2019.
* OKUNBOR, D. I.; SKEEL, R. D. Canonical Runge—Kutta—Nyström methods of orders five and six. Journal of Computational and Applied Mathematics, v. 51, n. 3, p. 375–382, jun. 1994.
* POTALEJ, O. A. Simulação numérica do problema de N-corpos gravitacional. 2024. Trabalho de Conclusão de Curso (Graduação) – Instituto de Matemática e Estatística, Universidade de São Paulo, São Paulo, 2024. Disponível em: https://doi.org/10.11606/003256212
