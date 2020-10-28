import pandas as pd
import numpy as np
from viterbi import viterbi
from Bio import SeqIO

input_dir = "./input/"
# leggo la sequenza fasta del gene APT
record = SeqIO.read(input_dir + "Q28748.fa", "fasta")
seq = list(record.seq)

# leggo il mapping tra AA idrofobici e idrofilici
amap = pd.read_csv(input_dir + "aa-hydro.txt", sep = ",", names = ["tmp"])
am = amap.T 

# stati nascosti e stati osservabili
S = np.array(["E", "M", "C"])
SY = np.array(["H", "L"])

# matrice di transizione
M = pd.DataFrame([[0.7, 0.3, 0.0], [0.25, 0.5, 0.25], [0.0, 0.3, 0.7]], columns = S, index = S)

# matrice probabilita' di emissione
E = pd.DataFrame([[0.4, 0.6], [0.85, 0.15], [0.5, 0.5]], columns = SY, index = S)

# probabilita' iniziali
pinizio = pd.DataFrame([[1., 0., 0.]], columns = S)

o = am[seq].to_numpy()[0]
m = viterbi(M,E,S,pinizio,o)

# confronto con i siti annotati della proteina Q28748
o0 = np.array(["E"] * 33 + ["M"] * (57 - 34 + 1) + ["C"])
somma = np.sum(o0 == m) / len(m)
tab = pd.crosstab(o0, m, rownames = ['o0'], colnames = ['m'])







'''


# leggo il recettore dell'insulina
record = SeqIO.read(input_dir + "insulin.fa", "fasta")
seq = list(record.seq)

# sequenza osservata
o = am[seq].to_numpy()[0]
m = viterbi(M,E,S,pinizio,o)

# confronto con i siti annotati
insulin = np.array(["E"] * 956 + ["M"] * (979-957+1) + ["C"] * (1382-980+1))
somma = np.sum(insulin == m) / len(m)
tab = pd.crosstab(insulin, m, rownames = ['insulin'], colnames = ['m'])

# all'inizio delle proteine di membrana c'e' quasi sempre un sito di binding
# proviamo ad eliminare i primi 50 aa
m = viterbi(M,E,S,pinizio,o[50:])
somma = np.sum(insulin[50:] == m) / len(m)
tab = pd.crosstab(insulin[50:], m, rownames = ['insulin'], colnames = ['m'])

# proviamo con un altra matrice di transizione
# in cui C stato finale e da M non si va mai in E
# matrice di transizione
M = pd.DataFrame([[0.99, 0.01, 0.], [0., 0.95, 0.05], [0., 0., 1.]], columns = S, index = S)

m = viterbi(M,E,S,pinizio,o[50:])
somma = np.sum(insulin[50:] == m) / len(m)
tab = pd.crosstab(insulin[50:], m, rownames = ['insulin'], colnames = ['m'])








# riconsidero la sequenza iniziale
# leggo la sequenza fasta del gene APT
record = SeqIO.read(input_dir + "Q28748.fa", "fasta")
seq = list(record.seq)

# sequenza osservata
o = am[seq].to_numpy()[0]
m = viterbi(M,E,S,pinizio,o)

# confronto con i siti annotati della proteina Q28748
o0 = np.array(["E"] * 33 + ["M"] * (57 - 34 + 1) + ["C"])
m[-1] = "C"
somma = np.sum(o0 == m) / len(m)
tab = pd.crosstab(o0, m, rownames = ['o0'], colnames = ['m'])
'''