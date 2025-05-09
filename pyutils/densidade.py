"""
Para facilitar a conversao entre densidade e raio.
A relacao entre a massa m, a densidade u e o volume V eh

    u = m / V.

Considernado particulas esfericas de raio r, seu volume eh

    V = 4 pi r^3 / 3.

Substituindo:

    u = 3 m / 4 pi r^3.

E vale a volta;

    r = cbrt(3 m / 4 pi u)

"""
from math import pi

print("< Calculo de densidade ou de raio, um em funcao do outro. >\n")
escolha = input("(1) Densidade em funcao do raio \n(2) Raio em funcao da densidade\n\n> ")

if escolha not in ["1", "2"]:
    print("\nEscolha inválida!")
    exit()

escolha = int(escolha)

if escolha == 1:
    densidade = input("\nInsira o valor da densidade:\n> ")
    try: densidade = float(densidade)
    except:
        print("\nValor estranho de densidade!")
        exit()
    
    massa = input("Insira o valor da massa:\n> ")
    try: massa = float(massa)
    except:
        print("\nValor estranho de massa!")
        exit()

    raio = (3 * massa / (4 * pi * densidade))**(1/3)
    print("\nRaio: ", raio)

elif escolha == 2:
    raio = input("\nInsira o valor do raio:\n> ")
    try: raio = float(raio)
    except:
        print("\nValor estranho de raio!")
        exit()
    
    massa = input("\nInsira o valor da massa:\n> ")
    try: massa = float(massa)
    except:
        print("\nValor estranho de massa!")
        exit()

    densidade = 3 * massa / (4 * pi * (raio**3))
    print("\nDensidade: ", densidade)