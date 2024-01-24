# "weather_seoul" este un fi??ier CSV cu starile vremii in Seoul, Coreea de Sud,
# in perioada 31.12.2019 - 20.04.2020

# 1. Stabiliti care este spatiul starilor S.

# Acesta este spatiul starilor:
stari_col <- unique(weather_seoul)
stari_col

# Acesta este numarul de stari:
numar_stari <- nrow(stari_col)
numar_stari

################################################################################

# 2. Estimati matricea de trecere P (probabilitatile sunt calculate in EXCEL).

# Introducem matricea de trecere P:
P <- matrix(c(0.564102564, 0.153846154, 0.076923077, 0.179487179, 0,           0.025641026,
              0.303030303, 0.575757576, 0.030303030, 0.090909091, 0,           0,
              0,           0.285714286, 0.428571429, 0.285714286, 0,           0,
              0.137931034, 0.206896552, 0,           0.586206897, 0.034482759, 0.034482759,
              1,           0,           0,           0,           0,           0,
              0.5,         0,           0,           0,           0,           0.5),
            nrow = 6,
            byrow = TRUE)
P

# Verificam ca suma pe linii este 1:
rowSums(P)

# Definim un vector care sa contina starile lantului Markov:
stari <- c("clear", "partly-cloudy", "rain", "cloudy", "snow", "wind")
stari

# Adaugam etichete liniilor si coloanelor matricei P:
rownames(P) <- stari
colnames(P) <- stari
P

################################################################################

# 10. Stabiliti daca lantul Markov este ergodic.

# Un lant Markov este ergodic daca si numai daca matricea de trecere este 
# regulara.
matrixpower <- function(matrice,k) {
  if (k == 0) return (diag(dim(matrice)[1]))
  if (k == 1) return(matrice)
  if (k > 1) return( matrice %*% matrixpower(matrice, k-1))
}

# Verificam daca matricea P este regulara, deci daca exista un numar natural
# s a.i. P^s > 0:
matrixpower(P, 10)

# P ridicata la puterea a zecea este o matrice pozitiva (P^10 > 0), astfel
# putem afirma ca P este o matrice regulara => lantul este ergodic.

################################################################################

# 11. Discutati existenta distributiei stationare, determinati-o.

# Daca matricea de tranzitie P este regulara, atunci lantul are o distributie
# limita, care este unica distributie stationara a lantului.

# Am aratat mai sus ca matricea P este regulara.

# Gasim distributia stationara folosind functia stationary.
stationary <- function(mat) {
  x = eigen(t(mat))$vectors[,1]
  as.double(x/sum(x))
}
distrib_stationara <- stationary(P)
distrib_stationara

################################################################################

# 12. Cercetati existenta distributiei limita a lantuluiMarkov.

# P este regulara => distrib. limita = distrib. stationara
distrib_limita <- distrib_stationara
distrib_limita

# Functie pentru a determina cea mai mica putere m a unei matrice de trecere  
# a unui lant Markov, cu proprietatea ca P^m=P^(m+1):

power_func <- function(matrice){
  for(m in 1:1000){
    mat1 <- round(matrixpower(matrice, m), digits = 9)
    mat2 <- round(matrixpower(matrice, m + 1), digits = 9)
    if(all.equal(mat1, mat2, tolerance = 0) == TRUE){
      return(m)
      break
    }
  }
}

power <- power_func(P)
power

# Matricea limita este:
P_lim <- matrixpower(P, 37)
P_lim

################################################################################

# 14. Propuneti o distributie initiala, iar apoi simulati o traiectorie a 
# lantului pe un anumit numar de pasi (un numar destul de mare de pasi sau 
# tranzitii), determinati distributia de frecvente absolute si relative a 
# intregii traiectorii, comparati distributia de frecvente relative a intregii 
# traiectorii simulate cu distributia stationara si cu distributia limita.

# Functia sample este o functie cu ajutorul careia se realizeaza sau se obtine
# o selectie cu revenire sau fara revenire, de volum precizat,obtinandu-se 
# valori sau realizari ale unei variabile aleatoare cu distributie de 
# probabilitate data sau cunoscuta.

set.seed(123454)
rezultate_simulare <- sample(x = c("clear", "partly-cloudy", "rain", "cloudy", "snow", "wind"),
                             size = 100,
                             replace = TRUE,
                             prob = c(1/6, 1/6, 1/6, 1/6,1/6, 1/6))
rezultate_simulare

# Distributia de frecvente absolute:
distr_frecv_abs <- table(rezultate_simulare)
distr_frecv_abs

# Distributia de frecvente relative:
distr_frecv_rel <- table(rezultate_simulare)/length(rezultate_simulare)
distr_frecv_rel 

################################################################################

# 15. Replicati simularea traiectoriei de 10000 de ori pentru acelasi orizont 
# de timp, considerand aceeasi distributie initiala ca mai sus, comparati 
# distributia de frecvente relative a ultimului moment considerat in simulare 
# cu distributia stationara si cu distributia limita.

# Introducem vectorul cu distributia de probabilitate initiala:
initial <- c(1/6, 1/6, 1/6, 1/6, 1/6, 1/6) 
initial

# Atribuim vectorului de probabilitate initial etichetele starilor:
names(initial) <- stari
initial

markov <- function(init, matrice, n, labels) {
  # daca starile nu sunt precizate, atunci starile sunt numerotate 1,....,m
  # numerele intregi consecutive de la 1 la lungimea lui init, adica numarul de stari
  if (missing(labels)) labels <- 1:length(init)
  # construieste un vector denumit simlist cu n+1 elemente, initializat cu toate elementele 0
  # functia numeric(numar intreg) creeaza un vector de lungime numar intreg cu toate elementele 0
  simlist <- numeric(n+1) 
  # starile sunt etichetele 1, 2, ..., k, unde k este numarul de stari al lantului
  # starile lantului sunt etichetele, adica numerele intregi consecutive
  states <- 1:length(init) 
  # primul element al vectorului simlist este o realizare sau simularea unei singure stari 
  # din distributia lui X0, adica din distributia initiala a lantului
  simlist[1] <- sample(states,
                       1,
                       replace = TRUE,
                       prob = init)
  # celelalte elemente ale vectorului simlist, adica simlist[2], ..., simlist[n+1]
  # reprezinta fiecare cate o realizare sau o simulare
  # din distributia de probabilitate data de linia simlist[i-1] a matricei P
  for (i in 2:(n+1)) 
  { simlist[i] <- sample(states,
                       1,
                       replace = TRUE,
                       prob = matrice[simlist[i-1],]) }
  # ataseaza starile corespunzatoare elementelor vectorului simlist:
  labels[simlist] 
}

# replicate() este o functie prin care se repeta o anumita actiune de un anumit 
# numar precizat de ori.

# replicam functia Markov definita anterior de 10000 ori:
sim_replicata <- replicate(10000,
                       markov(initial, P, 100, stari))
sim_replicata 

# sim_replicata reprezinta simularea a 10000 lanturi Markov, adica traiectoriile
# a 10000 traiectorii, de la momentul initial 0, pana peste 100 de zile.
# traiectoriile se gasesc pe fiecare din cele 100 de coloane.

# Distributia de frecvente absolute:
distr_frecv_abs_2 <- table(sim_replicata)
distr_frecv_abs_2

# Distributia de frecvente relative:
distr_frecv_rel_2 <- table(sim_replicata)/length(sim_replicata)
distr_frecv_rel_2 
