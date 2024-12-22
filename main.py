import random
from math import gcd
import numpy as np
import csv

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
def GLCG(M, lista_a, lista_x, n):
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

#testy
from scipy import stats
from scipy.stats import chisquare
from scipy.stats import norm

#Kolmogorov-Smirnov test
def KS(dane):
    statystyki = stats.kstest(dane, "uniform", alternative = "two-sided")
    return statystyki.pvalue

#chi2
def chi_2(dane):
    n_bins = np.linspace(0, 1, 11) ##sprawdzić jaka liczba zamiast 11 bedzie optymalna
    observed, bins = np.histogram(dane, bins=n_bins)

    expected = np.full(len(observed), len(dane) / len(observed))
    chi2_stat, p_chi = chisquare(f_obs=observed, f_exp=expected)
    return float(p_chi)


def count_zeros_ones(array):
    count_zeros = np.sum(array == 0)
    count_ones = np.sum(array == 1)
    return count_zeros, count_ones


from scipy.stats import cramervonmises

def cramer_von_mises_test(dane):
    result = cramervonmises(dane, 'uniform')
    return result.pvalue


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

#frequency monobit test
import math
from scipy.stats import norm
from scipy.special import erfc
def frequency_monobit_test(binary_digits):
    n = len(binary_digits)
    ones = np.sum(binary_digits == 1)
    zero = np.sum(binary_digits == 0)
    S_n = ones-zero
    test_statistic = S_n / math.sqrt(n)
    p_value = 2*(1-norm.cdf(abs(test_statistic)))
    return p_value


file_paths = {
    "pi": digits_pi,
    "e": digits_e,
    "sqrt2": digits_sqrt2
}
#monobit na digits_number
results_2 = {}
for name, digits_number in file_paths.items():
    p_value = frequency_monobit_test(digits_number)
    results_2[name] = [p_value]

for number, p_value in results_2.items():
    print(f"{number}: p-value = {p_value}")

from scipy.stats import chi2

def second_level_testing(p_values, s=10):
    R = len(p_values)
    intervals = [(i / s, (i + 1) / s) for i in range(s)]
    E = R / s
    O = [sum(1 for p in p_values if interval[0] <= p < interval[1]) for interval in intervals]
    chi_squared = sum((Oi - E) ** 2 / E for Oi in O)
    p_final = 1 - chi2.cdf(chi_squared, df=s - 1)

    return p_final

################################
###secon-level, monobit na fragmentach
def monobit_fragmenty(fragment_size = 1000):
    results_3 = {}
    for name, digits_number in file_paths.items():
        p_values = []

        for i in range(0, len(digits_number), fragment_size):
            fragment = digits_number[i:i + fragment_size]
            if len(fragment) < fragment_size:
                break
            p_value = frequency_monobit_test(fragment)
            p_values.append(p_value)
        # Zapisanie p-wartości do pliku CSV`
        output_file = f"{name}_p_values_{fragment_size}.csv"
        with open(output_file, mode='w', newline='') as file:
            writer = csv.writer(file)
            writer.writerow(["p_value"])
            writer.writerows([[p] for p in p_values])

        # Second-level testing na wektorze p-wartości
        result = second_level_testing(p_values)
        results_3[name] = result

    # Wyświetlanie wyników
    for number, stats in results_3.items():
        print(f"{number}: mean p-value = {stats:.4f}")

print('1000')
monobit_fragmenty(1000)
print('5000')
monobit_fragmenty(5000)
print('10000')
monobit_fragmenty(10000)


def generate_and_test_generators(n=2**20, repetitions=1000):
    results = {
        "python": {"KS": [], "CHI2": [], "LC": []},
        "RC4": {"KS": [], "CHI2": [], "LC": []},
        "LCG": {"KS": [], "CHI2": [], "LC": []},
        "GLCG": {"KS": [], "CHI2": [], "LC": []},
        "BBSG": {"KS": [], "CHI2": [], "LC": []},
    }
    # Generate sequences
    python_all = np.random.random(n*repetitions)
    rc4_data_all = RC4(32, [i for i in range(12)], n * repetitions)
    lcg_data_all = LCG(13, 1, 5, 0, n * repetitions)
    glcg_data_all = GLCG(1024, [3, 7, 68], [1, 670, 356], n * repetitions)
    bbsg_data_all = BBSG(499, 547, 12345, 16, n * repetitions)
    for i in range(repetitions):
        python_data = python_all[i * n:(i+1) * n]
        rc4_data = rc4_data_all[i * n:(i + 1) * n]
        lcg_data = lcg_data_all[i * n:(i + 1) * n]
        glcg_data = glcg_data_all[i * n:(i + 1) * n]
        bbsg_data = bbsg_data_all[i * n:(i + 1) * n]

        # Run tests and collect p-values
        results["python"]["KS"].append(KS(python_data))
        results["python"]["CHI2"].append(chi_2(python_data))
        results["python"]["LC"].append(cramer_von_mises_test(python_data))

        results["RC4"]["KS"].append(KS(rc4_data))
        results["RC4"]["CHI2"].append(chi_2(rc4_data))
        results["RC4"]["LC"].append(cramer_von_mises_test(rc4_data))

        results["LCG"]["KS"].append(KS(lcg_data))
        results["LCG"]["CHI2"].append(chi_2(lcg_data))
        results["LCG"]["LC"].append(cramer_von_mises_test(lcg_data))

        results["GLCG"]["KS"].append(KS(glcg_data))
        results["GLCG"]["CHI2"].append(chi_2(glcg_data))
        results["GLCG"]["LC"].append(cramer_von_mises_test(glcg_data))

        results["BBSG"]["KS"].append(KS(bbsg_data))
        results["BBSG"]["CHI2"].append(chi_2(bbsg_data))
        results["BBSG"]["LC"].append(cramer_von_mises_test(bbsg_data))

    # Save results to files
    for generator, tests in results.items():
        for test_name, p_values in tests.items():
            filename = f"{generator}_{test_name}_pvalues.txt"
            with open(filename, "w") as f:
                for p in p_values:
                    f.write(f"{p}\n")

    print("P-values saved to files.")


#generate_and_test_generators(2**15, 1000)

#wyniki tabela 1
import os
import glob

folder_path = "./"
files = glob.glob(os.path.join(folder_path, "*_pvalues.txt"))
data = {}

for file in files:
    generator_test = os.path.basename(file).replace("_pvalues.txt", "")
    with open(file, "r") as f:
        p_values = [float(line.strip()) for line in f]
    data[generator_test] = p_values

for key, values in data.items():
    result = second_level_testing(values)
    print(f"Klucz: {key}, Wynik: {result}")
