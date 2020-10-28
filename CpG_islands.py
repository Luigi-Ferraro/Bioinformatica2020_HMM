import numpy as np 
import pandas as pd 
from viterbi import viterbi


def find_CpG_islands_example1(gene):
    # Stati nascosti e osservabili
    S = np.array(["I", "N"])
    SY = np.array(["A", "C", "G", "T"])

    # Matrice transizione
    M = pd.DataFrame([[0.8, 0.2], [0.2, 0.8]], columns = S, index = S)

    # Matrice probabilita' di emissione
    E = pd.DataFrame([[0.1, 0.4, 0.4, 0.1], [1/4] * 4], columns = SY, index = S)

    # Probabilita' iniziali
    pinizio = pd.DataFrame([[0.5, 0.5]], columns = S)

    path = viterbi(M, E, S, pinizio, gene)

    for i in range(len(gene)):
        print(gene[i] + "  ", end =" ")
    print("\n")
    for i in range(len(path)):
        print(path[i] + "  ", end =" ")
    print("\n")


def find_CpG_islands_example2(gene):
    # Stati nascosti e osservabili
    S = np.array(["AI", "CI", "GI", "TI", "AN", "CN", "GN", "TN"])
    SY = np.array(["A", "C", "G", "T"])

    # Matrice transizione
    m = [[1.85152516e-01, 2.75974026e-01, 4.00289017e-01, 1.37026750e-01, 3.19045117e-04, 3.19045117e-04, 6.38090233e-04, 2.81510397e-04],
         [1.89303979e-01, 3.58523577e-01, 2.52868527e-01, 1.97836007e-01, 4.28792308e-04, 5.72766368e-04, 3.75584503e-05, 4.28792308e-04],
	     [1.72369088e-01, 3.29501650e-01, 3.55446538e-01, 1.40829292e-01, 3.39848138e-04, 4.94038497e-04, 7.64658311e-04, 2.54886104e-04],
     	 [9.38783432e-02, 3.40823149e-01, 3.75970400e-01, 1.86949063e-01, 2.56686367e-04, 5.57197235e-04, 1.05804868e-03, 5.07112091e-04],
	     [0.00000000e+00, 3.78291020e-05, 0.00000000e+00, 0.00000000e+00, 2.94813496e-01, 1.94641138e-01, 2.86962055e-01, 2.23545482e-01],
	     [0.00000000e+00, 7.57154865e-05, 0.00000000e+00, 0.00000000e+00, 3.26811872e-01, 2.94079570e-01, 6.17258712e-02, 3.17306971e-01],
	     [0.00000000e+00, 5.73810399e-05, 0.00000000e+00, 0.00000000e+00, 2.57133507e-01, 2.33483327e-01, 2.94234944e-01, 2.15090841e-01],
	     [0.00000000e+00, 3.11417347e-05, 0.00000000e+00, 0.00000000e+00, 1.79565378e-01, 2.32469115e-01, 2.94623408e-01, 2.93310958e-01]]
    M = pd.DataFrame(m, columns = S, index = S)

    # Matrice probabilita' di emissione
    d = np.eye(4)
    E = pd.DataFrame(np.concatenate([d,d]), columns = SY, index = S)

    # Probabilita' iniziali
    pinizio = pd.DataFrame([[1/8] * 8], columns = S)

    path = viterbi(M, E, S, pinizio, gene)

    for i in range(len(gene)):
        print(gene[i] + "  ", end =" ")
    print("\n")
    for i in range(len(path)):
        print(path[i][-1] + "  ", end =" ")
    print("\n")



    
if __name__ == "__main__":
    input_dir = "./input/"
    with open(input_dir + "gene1.txt", "r") as f:
        gene = np.array(list(f.read()))

    print("Primo esempio")
    find_CpG_islands_example1(gene)
    print("Secondo esempio")
    find_CpG_islands_example2(gene)
