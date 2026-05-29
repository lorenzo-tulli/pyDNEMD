from importlib.metadata import version, PackageNotFoundError

try:
    __version__ = version("pydnemd")
except PackageNotFoundError:
    __version__ = "unknown"
