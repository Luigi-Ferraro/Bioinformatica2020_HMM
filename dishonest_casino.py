import numpy as np 
import pandas as pd 
from viterbi import viterbi

def dishonest_casino():
    # Stati nascosti e osservabili
    S = np.array(["F", "L"])
    SY = np.arange(1,7)

    # Matrice transizione
    M = pd.DataFrame([[0.95, 0.05], [0.1, 0.9]], columns = S, index = S)

    # Matrice probabilita' di emissione
    E = pd.DataFrame([[1/6] * 6, [1/10] * 5 + [1/2]], columns = SY, index = S)

    # Probabilita' iniziali
    pinizio = pd.DataFrame([[1., 0.]], columns = S)

    # Sequenza di lanci osservata
    o = np.array([1,2,1,2,2,3,4,5,3,4,6,6,5,4,2,1,3,4,6,6,1,6,6,6,6,2,6,5,1,3,2,2,1])

    path = viterbi(M, E, S, pinizio, o)
    print(path)



if __name__ == "__main__":
    dishonest_casino()
