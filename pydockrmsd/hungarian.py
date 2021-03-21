from math import sqrt
from munkres import Munkres, DISALLOWED
from typing import List, Tuple


def readMol2(fname: str) -> Tuple[List[float],
                                  List[str]]:
    coords = []
    atoms = []
    with open(fname) as file:
        n = 0
        readflag = False
        for line in file:
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
            else:
                continue
        file.close()
    return (coords, atoms)


def hungarian(first_mol_path: str,
              second_mol_path: str):
    query = readMol2(first_mol_path)
    querycoords: List[float] = query[0]
    queryatoms: List[str] = query[1]
    temp = readMol2(second_mol_path)
    tempcoords: List[float] = temp[0]
    tempatoms: List[str] = temp[1]
    if len(querycoords) != len(tempcoords):
        Exception("identical atomcount needed")
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
    squaredist = 0.
    for row, column in Munkres().compute(distmat):
        squaredist += (float(distmat[row][column])/10000)/len(querycoords)

    return sqrt(squaredist)
