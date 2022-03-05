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
help(PyDockRMSD)

# %%
dockrmsd = PyDockRMSD("./data/targets/1a8i/crystal.mol2",
                      "./data/targets/1a8i/vina1.mol2")
if not dockrmsd.error:
    print(dockrmsd.rmsd)
    print(dockrmsd.total_of_possible_mappings)
    print(dockrmsd.optimal_mapping)
else:
    print(dockrmsd.error)
# %%
# dockrmsd = PyDockRMSD("./inputs_and_conversions/1ADQ_FV_deep.mol2",
#                       "./inputs_and_conversions/1adq_fixed_RF_short.mol2")
# if not dockrmsd.error:
#     print(dockrmsd.rmsd)
#     print(dockrmsd.total_of_possible_mappings)
#     print(dockrmsd.optimal_mapping)
# else:
#     print(dockrmsd.error)
# %%
print(hungarian("./data/targets/1a8i/crystal.mol2",
                "./data/targets/1a8i/vina1.mol2"))
