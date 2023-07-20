import sys

if __name__ == '__main__':
    with open(sys.argv[1]) as f:
        content = f.read()
        content = content.replace("\n", '')
        cts = content.split("=")
    T = []
    P = []
    S = []
    ind = 2
    while ind < (len(cts) - 6):
        T.append(float(cts[ind]))
        ind = ind + 2
        P.append(float(cts[ind]))
        ind = ind + 2
        S.append(float(cts[ind]))
        ind = ind + 2

    T.append(float(cts[ind]))
    ind = ind + 2
    P.append(float(cts[ind]))
    ind = ind + 2
    S.append(float(cts[ind].split('-')[0]))
    ind = ind + 2

    averageT = sum(T) / len(T)
    averageP = sum(P) / len(P)
    averageS = sum(S) / len(S)
    print("average TOTAL Time was: " + str(averageT) + "\naverage Time in PARALLEL was: " + str(averageP)+"\naverage time in SINGLEs was: "+str(averageS))