# gravidade-fortran

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

## 🧮 Simulando

**:warning: Uma interface para inserir as condições iniciais via arquivos está sendo construída :warning:**

Para determinar as condições iniciais do problema, é necessário escrever um "main.f90" no formato do deixado no repositório como exemplo. É possível gerar os valores como no "main_gerando.f90", no qual se determina intervalos para as massas, posições e momentos, além da quantidade de corpos. Outra possibilidade é inserir os dados manualmente, como no "main.f90". 

As massas, posições e momentos lineares são armazenados em arquivos .csv no diretório "data". Para evitar perda de dados, o nome do arquivo captura a data atual no formato "aaaammdd_vv.csv", onde "v" se refere à versão do dia, iniciando em 10 e indo até 99.

Uma forma de trabalhar os dados gerados é através do [gravidade-python](https://github.com/Potalej/gravidade-python), a versão inicial e em Python desse projeto. A ideia é que mais para frente seja possível acoplar esse visualizador em Python aqui, mas por enquanto segue a maracutaia.
