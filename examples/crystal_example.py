# %%
import logging
from rdkit import RDLogger
from rdkit.rdBase import DisableLog
import rdkit.RDLogger as rkl
from rdkit.rdBase import BlockLogs
# neudinger barre kevin library
from pydockrmsd.dockrmsd import PyDockRMSD
import pydockrmsd.hungarian as hungarian
logger = logging.getLogger(__name__)
# %%

logger = rkl.logger()
logger.setLevel(rkl.ERROR)
BlockLogs()
for level in RDLogger._levels:
    DisableLog(level)
# %%
# help(PyDockRMSD)

# %%
dockrmsd = PyDockRMSD("./data/targets/1a8i/crystal.mol2",
                      "./data/targets/1a8i/vina1.mol2")
# %%
print(dockrmsd.rmsd)
# %%
print(dockrmsd.total_of_possible_mappings)
# %%
print(dockrmsd.optimal_mapping)
# %%
print(dockrmsd.error)
# %%
print(hungarian("./data/targets/1a8i/crystal.mol2",
                "./data/targets/1a8i/vina1.mol2"))
