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

dane_LCG = LCG(13,1,5,1,10000)

#testy

from scipy import stats
from scipy.stats import chisquare
from scipy.stats import norm

#Kolmogorov-Smirnov test
def KS(dane):
    statystyki = stats.kstest(dane, "uniform", alternative = "two-sided")
    return statystyki.pvalue

#p_ks_LCG = stats.kstest(dane_LCG, "uniform", alternative = 'two-sided')
#k_odp = p_ks_LCG.pvalue
#print(p_ks_LCG)
#print(k_odp)


#chi2
def chi_2(dane):
    n_bins = np.linspace(0, 1, 11) ##sprawdzić jaka liczba zamiast 11 bedzie optymalna
    observed, bins = np.histogram(dane, bins=n_bins)

    expected = np.full(len(observed), len(dane) / len(observed))
    chi2_stat, p_chi = chisquare(f_obs=observed, f_exp=expected)
    return float(p_chi)

print(chi_2(dane_LCG))

#n_bins = np.linspace(0,1, 11)
#observed, bins = np.histogram(dane_LCG, bins=n_bins)

#expected = np.full(len(observed), len(dane_LCG) / len(observed))
#chi2_stat, p_chi_LCG = chisquare(f_obs=observed, f_exp=expected)
#print({'chi': p_chi_LCG})
#############
def count_zeros_ones(array):
    count_zeros = np.sum(array == 0)
    count_ones = np.sum(array == 1)
    return count_zeros, count_ones

from scipy.stats import norm

def run_test_continuous(data):
    # Tworzenie serii: 1 (rosnąca) i 0 (malejąca lub równa)
    runs = [1 if data[i] > data[i - 1] else 0 for i in range(1, len(data))]

    # Liczba serii
    total_runs = 1 + sum(runs[i] != runs[i - 1] for i in range(1, len(runs)))
    # Obliczanie średniej i odchylenia standardowego liczby serii
    n1 = sum(runs)
    n0 = len(runs) - n1
    mean_runs = (2 * n0 * n1) / (n0 + n1) + 1
    std_runs = np.sqrt((2 * n0 * n1 * (2 * n0 * n1 - n0 - n1)) / ((n0 + n1) ** 2 * (n0 + n1 - 1)))

    # Statystyka z
    z = (total_runs - mean_runs) / std_runs
    p_value = 2 * (1 - norm.cdf(abs(z)))  # Dwustronny p-value

    return p_value


#linear complexity test
def convert_to_binary(input_list):
    binary_list = [0 if 0 <= x < 0.5 else 1 for x in input_list]
    return binary_list

def berlekamp_massey_algorithm(sequence):
    n = len(sequence)
    c = np.zeros(n, dtype=int)
    b = np.zeros(n, dtype=int)
    c[0], b[0] = 1, 1
    l, m, d = 0, -1, 0

    for i in range(n):
        d = sequence[i] ^ np.dot(c[1:l + 1], sequence[i - l:i][::-1]) % 2
        if d == 1:
            t = c.copy()
            c[i - m:i - m + l + 1] ^= b[:l + 1]
            if l <= i // 2:
                l = i + 1 - l
                m = i
                b = t

    return l

def linear_complexity_test(input_list, block_size):
    binary_sequence = convert_to_binary(input_list)
    n = len(binary_sequence)
    num_blocks = n // block_size
    lc_values = []

    # Divide the sequence into blocks and calculate linear complexity
    for i in range(num_blocks):
        block = binary_sequence[i * block_size:(i + 1) * block_size]
        lc = berlekamp_massey_algorithm(block)
        lc_values.append(lc)

    # Expected linear complexity
    mean_lc = (block_size / 2) + ((-1) ** block_size) / 36
    variance = block_size * (1 / 2 - 1 / 18 + ((-1) ** block_size) * (1 / 36))

    # Calculate test statistic and p-value
    chi_squared = np.sum([(lc - mean_lc) ** 2 / variance for lc in lc_values])
    p_value = norm.sf(chi_squared)

    return p_value


# Przykładowe dane binarne
#block_size = 500
#p_value = linear_complexity_test(dane_LCG, block_size)
#print(f"P-value: {p_value}")

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

from scipy.stats import chi2

def second_level_testing(p_values, s=10):
    R = len(p_values)

    # Define the intervals Pi
    intervals = [(i / s, (i + 1) / s) for i in range(s)]

    # Calculate Ei (expected number of p-values in each interval)
    E = R / s

    # Calculate Oi (observed number of p-values in each interval)
    O = [sum(1 for p in p_values if interval[0] <= p < interval[1]) for interval in intervals]

    # Calculate chi-squared statistic
    chi_squared = sum((Oi - E) ** 2 / E for Oi in O)

    # Calculate p-value for the second-level test
    p_final = chi2.sf(chi_squared, df=s - 1)

    return p_final

################################
rc4_1 = RC4(31, [i for i in range(12)], 1000)
print(chi_2(rc4_1))
lcg_1 = LCG(13, 1, 5, 1, 1000)
print(chi_2(lcg_1))
print()

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
        results["python"]["LC"].append(linear_complexity_test(python_data, block_size=500))

        results["RC4"]["KS"].append(KS(rc4_data))
        results["RC4"]["CHI2"].append(chi_2(rc4_data))
        results["RC4"]["LC"].append(linear_complexity_test(rc4_data, block_size=500))

        results["LCG"]["KS"].append(KS(lcg_data))
        results["LCG"]["CHI2"].append(chi_2(lcg_data))
        results["LCG"]["LC"].append(linear_complexity_test(lcg_data, block_size=500))

        results["GLCG"]["KS"].append(KS(glcg_data))
        results["GLCG"]["CHI2"].append(chi_2(glcg_data))
        results["GLCG"]["LC"].append(linear_complexity_test(glcg_data, block_size=500))

        results["BBSG"]["KS"].append(KS(bbsg_data))
        results["BBSG"]["CHI2"].append(chi_2(bbsg_data))
        results["BBSG"]["LC"].append(linear_complexity_test(bbsg_data, block_size=500))

    # Save results to files
    for generator, tests in results.items():
        for test_name, p_values in tests.items():
            filename = f"{generator}_{test_name}_pvalues.txt"
            with open(filename, "w") as f:
                for p in p_values:
                    f.write(f"{p}\n")

    print("P-values saved to files.")


generate_and_test_generators(500, 20)

#wczytanie danych
import os
import glob

# Folder, gdzie znajdują się pliki
folder_path = "./"

# Wyszukanie wszystkich plików z końcówką "_pvalues.txt"
files = glob.glob(os.path.join(folder_path, "*_pvalues.txt"))

# Wczytanie danych z każdego pliku do słownika
data = {}

for file in files:
    generator_test = os.path.basename(file).replace("_pvalues.txt", "")  # Wyodrębnij nazwę generatora i testu
    with open(file, "r") as f:
        p_values = [float(line.strip()) for line in f]  # Konwersja każdej linii na float
    data[generator_test] = p_values  # Zapisz dane w słowniku

for key, values in data.items():
    result = second_level_testing(values)
    print(f"Klucz: {key}, Wynik: {result}")
