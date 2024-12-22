library(ggplot2)
setwd('C:/Users/roksa/OneDrive/Pulpit/magisterka 1/semestr3/monte/monte-carlo-project')

files <- list.files(pattern = "_pvalues\\.txt$")

# Funkcja do wczytywania plików
read_pvalues <- function(file) {
  data <- read.table(file, header = FALSE)
  test_info <- unlist(strsplit(file, "_"))
  generator <- test_info[1]
  test_name <- gsub("\\.txt$", "", test_info[2])
  data.frame(Generator = generator, Test = test_name, PValue = data$V1)
}

pvalues_data <- do.call(rbind, lapply(files, read_pvalues))

pvalues_data$Test[pvalues_data$Test == 'LC'] = 'CVM'
pvalues_data$Generator[pvalues_data$Generator == 'python'] = 'MT'

# Tworzenie wykresu
p_wartosi_wykres <- ggplot(pvalues_data, aes(x = PValue, fill = Test)) +
  geom_histogram(binwidth = 0.05, position = "dodge", alpha = 0.7) +  
  facet_wrap(~ Generator, scales = "free", ncol = 2) +  
  labs(title = "Histogramy p-wartości w zależności od testu",
       x = "PValue",
       y = "Częstotliwość") +
  theme_minimal() +
  theme(legend.position = "top") +
  xlim(-0.1,1)

ggsave("p_wartosc_a_test_2.png", plot = p_wartosi_wykres)
  
################pi, e, sqrt2
##dla 5000
pi_p_values5000 <- data.frame('p'=read.csv("pi_p_values_5000.csv"), 'generator' = 'pi')
e_p_values5000 <-  data.frame('p'=read.csv("e_p_values_5000.csv"), 'generator' = 'e')
sqrt2_p_values5000 <- data.frame('p'=read.csv("sqrt2_p_values_5000.csv"), 'generator' = 'sqrt2')

p_values_data5000 <- rbind(pi_p_values5000, e_p_values5000, sqrt2_p_values5000)


liczby_5000 <- ggplot(p_values_data5000, aes(x = p_value, fill = generator)) +
  geom_histogram(binwidth = 0.1, alpha = 0.7, position = "identity") +
  facet_wrap(~ generator, scales = "free_y") +  # Oddzielny histogram dla każdego generatora
  labs(title = "Histogramy p-wartości dla różnych generatorów liczbowych",
       subtitle = "dla próbek długości 5000",
       x = "p-wartość",
       y = "Częstość") +
  theme_minimal() +
  theme(legend.position = "none")
ggsave("p_liczby_5000.png", plot = liczby_5000)