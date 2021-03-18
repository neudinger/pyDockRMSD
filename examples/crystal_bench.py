# %%
import os
import glob
import rdkit
import pandas
import multiprocessing
import logging
import pathlib  # Read and close file one-liner
from functools import partial
from typing import Iterable
from rdkit import Chem
from rdkit import RDLogger
from rdkit.rdBase import DisableLog
import rdkit.RDLogger as rkl
from rdkit.rdBase import BlockLogs
# neudinger barre kevin library
from pydockrmsd.dockrmsd import PyDockRMSD
import pydockrmsd.hungarian as hungarian
script_dir = os.path.dirname(os.path.realpath(__file__))
logger = logging.getLogger(__name__)
# %%

logger = rkl.logger()
logger.setLevel(rkl.ERROR)
BlockLogs()
for level in RDLogger._levels:
    DisableLog(level)
# %%
help(PyDockRMSD)
# %%


def MolToJSON(mol: rdkit.Chem.Mol) -> str:
    return Chem.MolToJSON(mol) if mol else mol


def map_parallel(f, iter,
                 max_parallel=multiprocessing.cpu_count()) -> Iterable:
    """Just like map(f, iter) but each is done in a separate thread."""
    import sys
    import threading
    import traceback
    from queue import Queue, Empty
    total_items = 0
    queue = Queue()
    for i, arg in enumerate(iter):
        queue.put((i, arg))
        total_items += 1
    if max_parallel > total_items:
        max_parallel = total_items
    res = {}
    errors = {}

    class Worker(threading.Thread):
        def run(self):
            while not errors:
                try:
                    num, arg = queue.get(block=False)
                    try:
                        res[num] = f(arg)
                    except Exception:
                        errors[num] = sys.exc_info()
                except Empty:
                    break
    threads = [Worker() for _ in range(max_parallel)]
    [t.start() for t in threads]
    [t.join() for t in threads]
    if errors:
        if len(errors) > 1:
            logging.warning("map_parallel multiple errors: %d:\n%s" % (
                len(errors), errors))
        item_i = min(errors.keys())
        type, value, tb = errors[item_i]
        logging.info("map_parallel exception on item %s/%s:\n%s" % (
            item_i, total_items, "\n".join(traceback.format_tb(tb))))
        raise value
    return [res[i] for i in range(len(res))]

# %%


def rmsd_calculation(vinamol2_path: str,
                     crystal_mol2_path: str) -> pandas.DataFrame:
    vinamol: rdkit.Chem.Mol = Chem.rdmolfiles. \
        MolFromMol2File(vinamol2_path)
    DockRMSD: PyDockRMSD = PyDockRMSD(crystal_mol2_path, vinamol2_path)
    hungarian_rmsd: float = hungarian(crystal_mol2_path, vinamol2_path)
    result = {
        "crystal_mol2_path": crystal_mol2_path,
        "vinamol2_path": vinamol2_path,
        "docked_mol": vinamol,
        "hungarianRMSD": hungarian_rmsd,
        "DockRMSD": DockRMSD.rmsd,  # float
        "DockRMSD_optimal_mapping": DockRMSD.optimal_mapping,  # str
        "DockRMSD_total_of_possible_mappings":
        DockRMSD.total_of_possible_mappings,
        "DockRMSD_error": DockRMSD.error}
    return pandas.DataFrame(data=[result])


def compute(target: str):
    protein_path = glob.glob(f'{target}/*.pdb')
    crystal_mol2_path = glob.glob(f'{target}/crystal.mol2')
    target_df = pandas.DataFrame()
    if all([crystal_mol2_path, protein_path]):
        protein_path = protein_path[0]
        crystal_mol2_path = crystal_mol2_path[0]
        protein_id = protein_path.split("/")[-2]
        crystal_rmsd_calculation = partial(rmsd_calculation,
                                           crystal_mol2_path=crystal_mol2_path)
        mol2_files = glob.glob(f'{target}/vina[1-9]*.mol2')
        target_df = pandas.concat(map(crystal_rmsd_calculation, mol2_files))
        target_df["crystal_mol"] = Chem.rdmolfiles.MolFromMol2File(
            crystal_mol2_path)
        target_df["protein_id"] = protein_id
        target_df["protein"] = pathlib.Path(protein_path) \
            .read_text(encoding='ASCII')
    return target_df


# %%
def compute_example(targets: str) -> pandas.DataFrame:
    targets = glob.glob(f'{targets}*')
    result_df = pandas.concat(map_parallel(compute, targets))
    return result_df


# %%
result_df = compute_example(f"{script_dir}/data/targets/")
# %%
result_df["crystal_mol_json"] = result_df["crystal_mol"].map(MolToJSON)
result_df["docked_mol_json"] = result_df["docked_mol"].map(MolToJSON)
# %%
result_df
# %%
result_df.drop(columns=["crystal_mol", "docked_mol"], errors="ignore") \
    .to_parquet("dockrmsd_example.parquet")
# %%
result_df_2 = result_df.copy()
# %%
drop_col = ["crystal_mol", "docked_mol",
            "docked_mol_json", "crystal_mol_json",
            "crystal_mol2_path", "vinamol2_path", "protein"]
comparedf = (result_df.drop(columns=drop_col, errors="ignore") ==
             pandas.read_parquet("dockrmsd_example-sav.parquet")
             .drop(columns=drop_col, errors="ignore"))

# %%
if all(list(comparedf.all())):
    print("All RMSD was correctly computed")
# %%
