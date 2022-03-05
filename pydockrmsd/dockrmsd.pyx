import os
import cython
from libc.stdio cimport *  # noqa: E999

cdef extern from "stdio.h":
    # FILE * fopen ( const char * filename, const char * mode )
    FILE * fopen(const char * , const char * )  # noqa: E203, E202
    # int fclose ( FILE * stream )
    int fclose(FILE * )  # noqa: E203, E202

cdef extern from "./DockRMSD_sources/DockRMSD.c":
    # int grabAtomCount(FILE * , size_t * )  # noqa: E203, E202
    ctypedef struct DockRMSD:
        float rmsd
        float total_of_possible_mappings
        char * optimal_mapping
        char * error
    DockRMSD dock_rmsd(FILE * , FILE * )  # noqa: E203, E202


@cython.embedsignature(True)
@cython.binding(True)
cdef class PyDockRMSD:
    """PyDockRMSD
    DockRMSD: an open-source tool for atom mapping and RMSD calculation of symmetric molecules through graph isomorphism

    Cython Binding of DockRMSD:
        - reproduce ./DockRMSD file1.mol2 file2.mol2

    Parameters
    ----------

        first_mol_path: str
            os.path to the mol2 file

        second_mol_path: str
            os.path to the mol2 file

    Returns
    -------

        PyDockRMSD
            property:
            - rmsd : float
            - total_of_possible_mappings : float
            - optimal_mapping : str
            - error : str

    C file Written by Eric Bell \

    v1.0 written 5/2/2019 \

    Latest update (v1.1) written 8/26/2019 \

    Some modification was made to be entirely compatible with Cython Library. \

    To be more precise all float are not rounded unlike normal binary.

    """  # noqa: E501
    cdef DockRMSD data

    def __init__(self,
                 first_mol_path: str,
                 second_mol_path: str):
        first_mol_path_byte_string: bytes = first_mol_path.encode("UTF-8")
        cdef char * firstmolpath = first_mol_path_byte_string
        cdef FILE * first_cfile
        first_cfile = fopen(firstmolpath, "r")
        if first_cfile == NULL:
            raise FileNotFoundError(
                2, "No such file or directory: '%s'", first_mol_path)

        second_mol_path_byte_string: bytes = second_mol_path.encode("UTF-8")
        cdef char * secondmolpath = second_mol_path_byte_string
        cdef FILE * second_cfile
        second_cfile = fopen(secondmolpath, "r")
        if second_cfile == NULL:
            raise FileNotFoundError(
                2, "No such file or directory: '%s'", second_mol_path)
        self.data = dock_rmsd(first_cfile, second_cfile)

    @property
    def rmsd(self) -> float:
        """Return root mean square deviation of atomic positions : float"""
        return self.data.rmsd

    @property
    def total_of_possible_mappings(self) -> float:
        return self.data.total_of_possible_mappings

    @property
    def optimal_mapping(self) -> str:
        """Find the deterministically optimal mapping between query and template
        atoms, an exhaustive assignment search reminiscent of the VF2 algorithm
        coupled with Dead-End Elimination (DEE) is implemented."""
        return self.data.optimal_mapping.decode("UTF-8")

    @property
    def error(self) -> str:
        """Return empty str if no error was found: str"""
        return self.data.error.decode("UTF-8")