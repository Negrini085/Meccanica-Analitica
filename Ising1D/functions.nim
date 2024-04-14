import std/random
import strutils
import math

type 
    Ising* = object
        spins*: array[100, int]
        energia*: float32
        J*: float32
        h*: float32

# Procedura per l'implementazione delle Periodic Boundary Conditions
proc pbc(indice: int, n_spin: int): int =

    # Unico caso problematico è se mi trovo oltre (al di sopra) della lunghezza del vettore che stiamo prendendo in 
    # considerazione: per le PbC dobbiamo interagire con il primo spin del vettore contenitore
    if indice == n_spin:
        result = 0
    elif indice == -1:
        result = n_spin - 1
    else:
        result = indice



# Procedure per la creazione di un nuovo modello di Ising caratterizzato da temperatura infinita, ossia con spin orientati
# in ordine casuale: per fare questo inizializzo tutti gli spin up e poi faccio un numero elevato di ribaltamenti per ottenere
# una configurazione completamente casuale ---> scelgo casualmente sia l'indice della casella che l'orientazione dello spin
proc newIsing*(regime_J: float32, campo: float32): Ising = 
    # Inizializzo tutti gli spin ad uno
    for i in 0..<len(result.spins):
        result.spins[i] = 1
    # Effettuo le permutazioni necessarie    
    for i in 0..<3*len(result.spins):
        result.spins[rand(len(result.spins)-1)] = -1 + 2 * rand(1)

    # Calcolo l'energia del sistema nella configurazione iniziale
    result.h = campo; result.J = regime_J
    result.energia = 0
    for i in 0..<len(result.spins):
        # Devo fare attenzione all'ultima iterazione, in quanto vogliamo che vengano applicate delle boundaries conditions per
        # far sì che i risultati siano confrontabili con quelli che si ottengono nel limite termodinamico, ossia per numero di 
        # particelle facenti parte del sistema in analisi tendente ad infinito.
        result.energia = result.energia - result.J * float32(result.spins[pbc(i, len(result.spins))] * result.spins[pbc(i+1, len(result.spins))]) - 0.5 * result.h * float32((result.spins[pbc(i, len(result.spins))] + result.spins[pbc(i+1, len(result.spins))]))



# Procedura utile per lavorare con le variazioni di energia legate agli spinflips
proc boltzmann(modello: Ising, ind, sval : int): float32 = 

    # Tengo in considerazione sono i valori che cambiano in base allo spin flip, quindi i primi vicini dato che stiamo
    # considerando una tale tipologia di interazione
    result = - modello.J * float32(sval * (modello.spins[pbc(ind - 1, len(modello.spins))] + modello.spins[pbc(ind + 1, len(modello.spins))])) - modello.h * float32(sval)

# Per evolvere il sistema lavoriamo con l'algoritmo di Metropolis, che consente di campionare il peso statistico di Gibbs:
# scegliamo casualmente uno degli spin facenti parte del modello di Ising 1D in analisi e ne valutiamo il costo energetico
# d'inversione, ossia quanto cambia il valore dell'energia. Il flip viene accettato con una probabilità che è il minimo fra
# 1 oppure il peso di Gibbs: una mossa MonteCarlo di tipo Metropolis consiste nel tentare N_spins flip.
proc metropolisMove(modello: var Ising, beta: float32): void =

    var
        ind, spin: int
        en1, en2, prob: float32

    # Cerco di effettuare N flip
    for i in 0..<len(modello.spins):
        ind = rand(len(modello.spins)-1)
        # Testo che cosa implica effettuare lo spin (che variazione ho dell'energia interna?)
        # Per fare questo provo a fare il flip e vedo come cambia il contributo di energia legato al
        # singolo spin
        spin = -1 * modello.spins[ind]
        en1 = modello.boltzmann(ind, modello.spins[ind])
        en2 = modello.boltzmann(ind, spin)

        # Valuto la probabilità e scelgo se accettare lo spin flip oppure no
        prob = min(1.0, exp(- beta * (en2 - en1)));
        if rand(1.0) < prob:
            modello.spins[ind] = spin
            modello.energia = modello.energia + en2 - en1


var
    modello: Ising = newIsing(1.0, 0.0)
    beta: float32 = 2

for i in 0..<300:
    modello.metropolisMove(beta)

echo modello