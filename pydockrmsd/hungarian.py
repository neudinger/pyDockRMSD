from math import sqrt
from munkres import Munkres, DISALLOWED


def readMol2(fname):
    f = open(fname)
    coords = []
    atoms = []
    readflag = False
    labels = []
    n = 0
    for line in f:
        if line.strip() == "@<TRIPOS>BOND":
            break
        elif line.strip() == "@<TRIPOS>ATOM":
            readflag = True
        elif readflag:
            n += 1
            parts = line.split()
            if parts[5] == "H":
                continue
            coords.append([float(i) for i in parts[2:5]])
            atoms.append(parts[5].split('.')[0])
            labels.append(parts[5]+str(n))
        else:
            continue
    return (coords, atoms, labels)


def hungarian(first_mol_path: str,
              second_mol_path: str):
    query = readMol2(first_mol_path)
    querycoords = query[0]
    queryatoms = query[1]
    temp = readMol2(second_mol_path)
    tempcoords = temp[0]
    tempatoms = temp[1]
    if len(querycoords) != len(tempcoords):
        print("identical atomcount needed")
        exit(1)

    distmat = [[0. for j in range(len(tempcoords))]
               for i in range(len(querycoords))]
    showmat = [[0. for j in range(len(tempcoords))]
               for i in range(len(querycoords))]
    for i in range(len(querycoords)):
        for j in range(len(tempcoords)):
            if queryatoms[i] == tempatoms[j]:
                distmat[i][j] = int(
                    sum([(querycoords[i][k]-tempcoords[j][k])**2 for k in
                         range(3)]*10000))
                showmat[i][j] = distmat[i][j]
            else:
                distmat[i][j] = DISALLOWED
                showmat[i][j] = 0.0

    m = Munkres()
    # try:
    indexes = m.compute(distmat)
    # except Exception as e:
    #     print(e)
    #     exit(1)

    squaredist = 0.
    for row, column in indexes:
        squaredist += (float(distmat[row][column])/10000)/len(querycoords)

    return sqrt(squaredist)
