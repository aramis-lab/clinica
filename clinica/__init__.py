try:
    from importlib.metadata import version
except ImportError:
    # Python 3.7 backport module.
    from importlib_metadata import version

__all__ = ["__version__"]

__version__ = version("clinica")
