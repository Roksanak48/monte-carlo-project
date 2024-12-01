import random
from math import gcd

#funkcja do skalowania otrzymanych liczb z genaratorów do przedziału (0,1)
def skaluj_do_0_1(lista_liczb, M):
    return [x / M for x in lista_liczb]


#generatory
#RC4(32) - K to lista
def RC4(m, K, n):
    x = []
    S = []
    L = len(K)
    #KSA(K)
    for i in range(m):
        S.append(i)
    j = 0
    for i in range(m):
        j = (j + S[i] + K[i % L]) % m
        S[i], S[j] = S[j], S[i]
    #PRGA
    i, j = 0, 0
    for k in range(n):
        i = (i + 1) % m
        j = (j + S[i]) % m
        S[i], S[j] = S[j], S[i]
        x.append(S[(S[i] + S[j]) % m])
    return skaluj_do_0_1(x, m)

#LCG(13,1,5)
def LCG(M,a,c, x_0, n):
    x = [x_0]
    for i in range(1, n):
        x.append((a*x[-1] + c) % M)
    return skaluj_do_0_1(x, M)

#GLCG(2^10,{3,7,68})
#sprobuj zrobić gdy lista_a ma więcej niż 3 elementy. może użyć słownik??
def GLCG(M, lista_a, lista_x, n):
    for i in range(3,n):
        lista_x.append((lista_a[0]*lista_x[-1] + lista_a[1]*lista_x[-2] + lista_a[2]*lista_x[-3]) % M)
    return skaluj_do_0_1(lista_x, M)

def GLCG2(M, lista_a, lista_x, n):
    k = len(lista_a)
    if len(lista_x) != k:
        raise ValueError("Długość lista_x musi być równa długości lista_a")

    for i in range(k, n):
        suma = sum(lista_a[j] * lista_x[-(j+1)] for j in range(k))
        lista_x.append(suma % M)

    return skaluj_do_0_1(lista_x, M)

def BBSG(p, q, ziarno, liczba_bitow, n):
    if p % 4 != 3 or q % 4 != 3:
        raise ValueError("p i q muszą być liczbami pierwszymi i ≡ 3 (mod 4)")

    M = p * q

    if ziarno is None:
        raise ValueError("Ziarno nie może być puste")
    if gcd(ziarno, M) != 1:
        raise ValueError("Ziarno musi być względnie pierwsze z M=p*q.")

    wyniki = []

    for i in range(n):
        ziarno = (ziarno ** 2) % M
        wynik = ziarno & ((1 << liczba_bitow) - 1)
        wyniki.append(wynik)

    return skaluj_do_0_1(wyniki, (2 ** liczba_bitow - 1))


#sprawdzenie generatorów
print(LCG(13,1,5,1,10))
print(GLCG(1024, [3,7,68], [1,670,356], 10))
print(GLCG2(1024, [3,7,68], [1,670,356], 10))
print(RC4(32,[i for i in range(12)],10))
print(BBSG(499, 547, 12345, 16, 10))


#testy



