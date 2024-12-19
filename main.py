import random
from math import gcd
import numpy as np

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
#print(LCG(13,1,5,1,10))
#print(GLCG2(1024, [3,7,68], [1,670,356], 10))
#print(RC4(32,[i for i in range(12)],10))
#print(BBSG(499, 547, 12345, 16, 10))

dane_LCG = LCG(13,1,5,1,10000)

#testy

from scipy import stats
from scipy.stats import chisquare

#Kolmogorov-Smirnov test
p_ks_LCG = stats.kstest(dane_LCG, "uniform", alternative = 'two-sided')
print(p_ks_LCG)

#chi2
n_bins = np.linspace(0,1, 11)
observed, bins = np.histogram(dane_LCG, bins=n_bins)

expected = np.full(len(observed), len(dane_LCG) / len(observed))
chi2_stat, p_chi_LCG = chisquare(f_obs=observed, f_exp=expected)
print(p_chi_LCG)

#wczytanie liczb pi, e ,pierwiastek z 2 z pliku
import urllib.request

def read_digits(url):
    data = []
    with urllib.request.urlopen(url) as f:
        for line in f:
            data.append(line.strip())
    datastring = []

    for line in data:
        datastring.append(line.decode("utf-8"))

    datastring = ''.join(datastring)
    datastring = list(map(int, list(datastring)))

    return (np.array(datastring))


digits_pi = read_digits('http://www.math.uni.wroc.pl/~rolski/Zajecia/data.pi')
digits_e = read_digits('http://www.math.uni.wroc.pl/~rolski/Zajecia/data.e')
digits_sqrt2 = read_digits('http://www.math.uni.wroc.pl/~rolski/Zajecia/data.sqrt2')
print("wczytane")

#frequency mobil test
import math
from scipy.stats import norm
from scipy.special import erfc
def frequency_monobit_test(binary_digits):
    n = len(binary_digits)
    ones = np.sum(binary_digits == 1)
    zero = np.sum(binary_digits == 0)
    #s = np.array([1 if bit == '1' else -1 for bit in binary_digits])
    #S_n = np.sum(s)
    S_n = ones-zero
    test_statistic = S_n / math.sqrt(n)
    # Calculate the p-value using the complementary error function
    p_value = 2*(1-norm.cdf(abs(test_statistic)))
    return p_value


file_paths = {
    "pi": digits_pi,
    "e": digits_e,
    "sqrt2": digits_sqrt2
}

results_2 = {}
for name, digits_number in file_paths.items():
    p_value = frequency_monobit_test(digits_number)
    results_2[name] = [p_value]

for number, p_value in results_2.items():
    print(f"{number}: p-value = {p_value}")