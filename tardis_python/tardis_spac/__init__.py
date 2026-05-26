from importlib.metadata import PackageNotFoundError, version

from . import stats, utils

try:
    from . import external
except ImportError:
    pass

__author__ = __maintainer__ = "zenglabPKU"
__email__ = "barry_2001@pku.edu.cn"

try:
    __version__ = version("tardis_spac")
except PackageNotFoundError:
    __version__ = "0.0.0"