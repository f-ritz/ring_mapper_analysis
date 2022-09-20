import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


def main():
    f = open(
            "processed/minittr-6-2HP-ref-16.5-NaC-0-Mg_C0117_DMS/output/BitVector_Files/minittr-6-2HP-ref_bitvectors.txt")
    lines = f.readlines()
    f.close()
    dict = {}
    lines = lines[4:]
    a = np.zeros((173, 173))
    b = np.zeros((173, 173))
    c = np.zeros((173, 173))
    d = np.zeros((173, 173))
    n = 0
    for l in lines:
        n += 1
        if n > 1000:
            break
        data = np.array(list(l.split('\t')[1]))
        for i, e in enumerate(data):
            mut_i = 0
            if e == 'A' or e == 'C' or e == 'G' or e == 'T':
                mut_i = 1
            for j in range(i + 1, len(data)):
                mut_j = 0
                if data[j] == 'A' or data[j] == 'C' or data[j] == 'G' or data[j] == 'T':
                    mut_j = 1
                if mut_i == 0 and mut_j == 0:
                    a[i][j] += 1
                elif mut_i == 1 and mut_j == 0:
                    b[i][j] += 1
                elif mut_i == 0 and mut_j == 1:
                    c[i][j] += 1
                else:
                    d[i][j] += 1

    for i in range(0, 173):
        for j in range(i + 1, 173):
            dim = (a[i][j] + b[i][j]) * (c[i][j] + d[i][j]) * (a[i][j] + c[i][j]) * (b[i][j] + d[i][j])
            n = a[i][j] + b[i][j] + c[i][j] + d[i][j]
            if dim == 0:
                continue
            chisq = (n * (abs(a[i][j] * d[i][j] - b[i][j] * c[i][j]) - 0.5 * n) ** 2) / dim
            if chisq > 20:
                print(i + 1, j + 1, chisq)


    # plt.imshow(mat, interpolation='none', cmap='Greys')
    # plt.colorbar()
    # plt.show()


if __name__ == '__main__':
    main()
