__author__ = "Barre Kevin"
__maintainer__ = "Barre kevin"
__credits__ = ["https://kevinbarre.fr/"]
__email__ = "kevin.barre@epitech.eu"
__status__ = "Production"

from .hungarian import hungarian

try:
    from .__version__ import __version__  # noqa: F401, F40
except Exception:
    pass

__all__ = [
    "hungarian",
    "dockrmsd"
]
