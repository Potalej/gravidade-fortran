# LOG.md - Diário de Desenvolvimento

## [2025-11-10] v0.8.1: Opção de pasta de saída (issue #15)

Agora é possível escolher a pasta de "out". Basta passar o parâmetro `-ps destino` ou `--pasta-saida destino`.

---

## [2025-10-14] v0.8.0: Nova forma de armazenar dados

Após sugestão, troquei o uso de arquivos CSV por arquivos binários na saída de dados. Isso traz uma redução considerável no tamanho dos arquivos de saída, bem como tem um impacto positivo significativo em simulações com muitos salvamentos.

---

## [2025-08-13] v0.7.2: Correção na inicialização da simulação

Na modularização não adicionei as variáveis que guardam se as massas são iguais, o que impacta no output de informações da simulação.

---

## [2025-08-11] v0.7.1: Enviando o tempo na exibição

Para facilitar comparações, o instante $t$ agora é enviado pelo socket quando o "exibir" estiver habilitado.

---

## [2025-08-08] v0.7.0: Maior modularização do programa

- Boa parte do programa foi refatorada e reorganizada para tudo ficar mais modularizado. Isso incluiu mudanças no CMakeLists principal e na criação de um secundário para os integradores.
- Alguns módulos foram transformados em bibliotecas, o que aumenta a centralização e ao mesmo tempo reduz a dependência entre módulos, o que permite testar melhor e futuramente fazer alterações sem tanta dificuldade.
- O integrador-base dos métodos de Runge-Kutta também foi suprimido e os métodos de Runge-Kutta disponíveis (de passo fixo) foram adaptados para a estrutura comum de integração.
- Agora adicionar integradores ficou mais fácil:
  - Adicione o integrador em algum subdiretório de "src/calcs/integracao";
  - Adicione a referência em "src/calcs/integracao/CMakeLists.txt".

---


## [2025-08-03] v0.6.2: Correções e melhorias

- "src/libs/aleatorio.f90": Para a geração de valores iniciais, agora o "raio" não é mais um parâmetro, sendo calculado em função do intervalo passado. Além disso, se o intervalo não for simétrico na origem, os valores gerados também não serão, com a devida translação.
- "src/calcs/mecanica/colisao.f90": Não é necessário calcular a norma do vetor normal, apenas o seu quadrado, então isso foi removido.
- "src/calcs/condicionamento/condicoesIniciais.f90": Código adaptado para o módulo "aleatorio", e também adiciona a opção `sorteio`, que não condiciona os valores iniciais, simplesmente os gera.
- "src/calcs/condicionamento/condicionamento.f90": Adaptação para o módulo "aleatorio", melhorias no processo de condicionamento direto, inclusão de condições de Sundman e de Delta.

---

## [2025-07-26] v0.6.1: Mais sobre o potencial amortecido e o virial

Na v0.5.5 corrigi a geração de valores iniciais no sentido de usar o potencial amortecido quando necessário. Ocorre que o cálculo do virial estava incorreto, porque o potencial deixa de ser homogêneo, então a expressão correta para o virial é:

$$
\dfrac{dG}{dt} = 2T + \sum_{a=1}^N \langle \vec F_a(\vec q), \vec q_a \rangle
$$

Para mais detalhes, [veja a wikipedia](https://en.m.wikipedia.org/wiki/Virial_theorem#Statement_and_derivation).

Isso também me obrigou a mudar o condicionamento do Aarseth, já que o interesse principal é no equilíbrio, e não exatamente nos valores de $T$ e $V$. Adicionei a rotina do Aarseth Modificado para fazer isso, sendo idêntico ao Aarseth no caso em que não há amortecimento.

---

## [2025-07-24] v0.6.0: GPU + N-corpos? Vamos experimentar

Adicionei um primeiro teste de uso de GPU neste programa, usando o OpenMP offloads para placas NVidia. Ainda não entendo direito como isso funciona, mas pelo menos já funciona. Com o tempo vou aprendendo e melhorando. Evidentemente, está em fase experimental. Mais para frente penso em tentar usar CUDA se parecer valer a pena. Antes, preciso aprender o básico.

Por ora, para compilar usando a GPU é só passar a flag

```shell
-DUSAR_GPU=ON
```

que só vai funcionar também para placas da NVidia.

Ah, também coloquei as funções de força numa pasta separada, finalmente.

---

## [2025-07-24] v0.5.6: Melhorias - Usando melhor o json-fortran e massas normalizadas

Quando adicionei o JSON-Fortran algumas versões atrás, o utilizei quase que somente no script de simulação, e estava passando os parâmetros manualmente para o integrador e outras rotinas. Isso não faz sentido. Simplesmente melhorei a escrita do código, o que acaba facilitando adicionar funcionalidades no futuro.

Das massas normalizadas, estava confusa a entrada de dados especialmente para sorteio de valores iniciais. Estava algo assim:

```json
"massas_iguais": true,
"sorteio":{
  "massas": {
    "intervalo": [1.0, 2.0],
    "distribuicao": "uniforme"
  }
  ...
}
```

Não era claro o que isso deveria fazer. Agora melhorei:

```json
"massas_iguais": true,
"sorteio":{
  "massas": {
    "normalizadas": false,
    "intervalo": [1.0, 2.0],
    "distribuicao": "uniforme"
  }
  ...
}
```

1. A opção `massas_iguais` se aplica somente à otimização das funções de força e dos integradores, e não depende das massas serem de fato iguais.

2. Se a opção `normalizadas` estiver ativada, as outras opções de massas serão ignoradas e as massas serão tomadas como 1/N.

3. Se estiver desativada, o intervalo será levado em conta. Se seu mínimo e seu máximo forem iguais, as massas serão iguais a esse valor e a distribuição não será levada em conta.

4. Se o intervalo tiver mínimo e máximo diferentes, a distribuição será utilizada para sortear os valores nesse intervalo.

---

## [2025-07-21] v0.5.5: Correção - Amortecimento do potencial no condicionamento

Na v0.5.4 incorporei corretamente a energia amortecida em quase todo o código, mas faltou incorporá-la de fato no processo de condicionamento de valores iniciais, pois me baseei completamente na hipótese do potencial ser homogêneo de grau -1, o que não ocorre no amortecimento. Corrigi isso no condicionamento indireto e no de Aarseth adicionando a subrotina `condicionar_potencial_amortecido`, que aplica o método de Newton para resolver o problema

$$f(h) = V_1 (hq_0) - \varepsilon (E - T),$$

e a posição condicionada final é \(h \epsilon q_0\).

O condicionamento direto depende fortemente do potencial ser homogêneo, e ao que me consta trata-se de um problema implícito tridimensional (alpha, beta, omega), o que necessita de um Newton complicado demais por agora. Assim, por enquanto, o condicionamento indireto lança um aviso quando usado com \(\varepsilon \neq 0\), indicando que não funciona, e o de Aarseth só o usa quando não há amortecimento.

---

## [2025-07-20] v0.5.4: Correção - Energia amortecida

Quando adicionamos um amortecimento na força também precisamos incluir isso no potencial, mas por algum motivo isso não estava feito na versão sincronizada, então adicionei. Agora o potencial tem um \(\varepsilon\) também!

---

## [2025-07-12] v0.5.3: Correções nos RK's, melhorias nos INTENTs e atualizações

- Por falta de uso não percebi que os Runge-Kuttas estavam com problemas. Corrigi.

- Arrumei o uso de INTENTs, que principalmente no caso dos OUTs estava definido como INOUT. Isso faz alguma diferença, provavelmente.

---

## [2025-07-12] v0.5.2: Plot em tempo real em Linux também!

Pronto, agora a atualização da v0.5.0 também funciona em Linux, hehe

---

## [2025-07-11] v0.5.1: Correções e remoção de variáveis inúteis

Ativei uma flag de depuração por curiosidade e centenas de linhas de avisos apareceram. Fiz questão de resolver a maioria deles, principalmente removendo variáveis que declarei mas não usei ou parei de usar em algum momento e esqueci lá. É isso.

---

## [2025-07-11] v0.5.0: Plotando em tempo real via sockets!

As pessoas sempre ficam curiosas sobre o que exatamente faço, e dificilmente consigo mostrar já que é chatíssimo o processo de ler os dados de uma simulação e gerar um vídeo usando Python.

Para resolver isso, implementei o valor booleano `exibir` nos dados iniciais, que quando ativado se conecta a um socket e envia dados de posição (futuramente de velocidade também) continuamente. Com o python, por exemplo, conseguimos capturar esses dados e plotá-los usando alguma biblioteca otimizada para esse tipo de coisa.

É evidente que o programa fica muito mais lento. Mas é legal para ver mais ou menos o que acontece em algum determinado instante da simulação. E as pessoas gostam de ver centenas ou milhares de bolinhas saltitando pela tela. Confesso que também gosto :)

Ah, e por enquanto só funciona no Windows. Vou resolver isso daqui uns dias.

---

## [2025-07-10] v0.4.7: Otimização significativa no cálculo da força e GPROF

- Por algum motivo, separar o cálculo da força entre dois corpos em uma outra função reduziu o tempo do RKN551 entre 30% e 50%. Inclusive, o desempenho utilizando uma `FUNCTION` é melhor que o utilizando uma `SUBROUTINE`. Vai entender...

- Para analisar o desempenho do programa incorporei o GPROF. Para habilitá-lo, basta compilar o programa utilizando a flag `-DGPROF=ON`.

---

## [2025-07-10] v0.4.6: Vetor global de distâncias

As distâncias entre os corpos são utilizadas em pelo menos quatro locais diferentes: na integração, nos choques elásticos, na correção e no cálculo do potencial, geralmente nessa ordem. Assim, o módulo `integracao` agora tem um vetor de distâncias que é atualizado durante o cálculo das forças (na integração) e é utilizado pelas outras funções. Isso reduz consideravelmente o custo das outras operações, que são todas \(N(N-1)/2\) a princípio.

---

## [2025-06-09] v0.4.5: Agora a flag PRECISAO também se aplica ao JSON-Fortran

O problema mencionado nas versão 0.4.0 e corrigido com certa tramóia na versão 0.4.2 sobre o JSON-Fortran por padrão vir instalado no sistema em REAL64 agora foi corrigido completamente. A flag PRECISAO (que pode ser utilizada na compilação como `-DPRECISAO=X`) agora também se aplica ao JSON-Fortran, o que permite utilizar precisão dupla ou quádrupla (acho que simples também) sem problemas.

---

## [2025-06-09] v0.4.4: Não é muito inteligente não manter uma versão local de uma biblioteca

Baixei a versão 9.0.3 (fevereiro/2025) do JSON-Fortran e a coloquei localmente para ser compilada junto do resto do programa. Não era possível usá-la na Rede IME por ser necessário instalá-la como pacote, e assim também faz mais sentido. Tudo certo.

---

## [2025-06-03] v.0.4.3: Maiores cuidados com os coeficientes de integração e implementação do Apêndice A de (FUKUSHIMA,2000)

- Coeficientes definidos em 128 bits para garantir a máxima precisão dos métodos com ordem maior que 2;

- Implementação da sugestão de Toshio Fukushima em "Reduction of round-off error in symplectic integrators" no método de Verlet. Futuramente cabe fazer tal análise e implementação nos outros métodos simpléticos.

---

## [2025-06-01] v0.4.2: Otimizações para simulações com massas iguais

Quando tomamos massas iguais (geralmente \(m=1/N\), mas nem sempre) muitos cálculos podem ser feitos mais diretamente, e às vezes até vetorialmente, então foram adicionadas funções alternativas que são chamadas quando as massas são todas iguais.

---

## [2025-05-31] v0.4.1: Tipos globais, versionamento e mais

- Foi criado um módulo "tipos" que contém os tipos principais utilizados e também define o tipo global REAL(pf)
  a partir da flag -DPRECISAO=X, onde X=32, 64 ou 128 (o padrão é 64);
- O json-fortran por padrão do MinGW vem compilado para REAL(64), então foi necessário adaptar o programa para
  lidar com qualquer tipo de entrada e saída com o json-fortran. A solução temporária foi interagir com o
  json-fortran somente através de REAL(64) e prévia ou posteriormente transformar os dados de ou para REAL(128).
  Isso pode limitar a precisão em simulações que exigem muito cuidado, mas por hora é suficiente.
- Versionamento agora aparece no "help".

---

## [2025-05-09] v0.4.0: Agora usamos JSON-Fortran!

Implementado o uso da biblioteca JSON-Fortran para a leitura e escrita de arquivos de condicionamento e de valores iniciais. Foi necessário atualizar toda a estrutura do programa, o que justifica a mudança de versão. Ainda falta permitir o uso de diferentes tipos de FLOAT.

---

## [2025-05-01] v0.3.2: Testes de relaxamento e anisotropia

Implementado o cálculo de anisotropia do tensor de inércia a partir de seus autovalores.

---

## [2025-04-08] v0.3.1: Melhorias no cálculo das colisões elásticas

Foi implementada uma melhoria significativa no cálculo das colisões elásticas, substituindo a necessidade do cálculo da base contendo o versor normal por um cálculo direto e mais eficiente. Para mais detalhes, ver a documentação técnica.

---

## [2025-03-16] v0.3.0: Octree e mais

Implementada a possibilidade de uso de árvores octanárias para verificação de colisões.

Embora teoricamente o método via octrees seja mais eficiente que a verificação direta,
são necessárias melhorias no código para que a eficiência possa ser observada. Mas funciona! Futuramente farei as melhorias devidas.

Também foi implementado o método descrito por Aarseth para geração de valores iniciais
sob condição de virial. Esta seção será reorganizada futuramente mas também funciona.

---

## [2025-03-02] v0.2.1: Correção na paralelização das forças

Neste dia comecei a fazer algum tipo de versionamento no programa, então vou contar a partir daqui.

A primeira atualização de 0.1 para 0.2.1 foi uma correção na paralelização ads forças, fazendo os \(N^2\) cálculos em vez de \(N(N-1)/2\) no caso paralelizado para garantir que não haja sobreescrição de dados.

---

## [2023-02-22]

- Commit inicial! Vamos lá.
