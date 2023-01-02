"""This module contains utility functions which are used across Pydra workflows."""
import typing as ty


def sanitize_fwhm(
    fwhm: ty.Union[float, ty.Tuple[float], ty.List[float], ty.List[ty.List[float]]],
) -> ty.List[ty.Tuple[float, float, float]]:
    """Make sure the FWHM is in the right format for the Smooth SPM interface.

    Parameters
    ----------
    fwhm : Union[float, Tuple[float], List[float], List[List[float]]
        Smoothing kernel(s) that should get passed to the SPM Smooth() interface.
        There are three ways to specify fwhm:
            - A float
            - A list/tuple of floats
            - A list of lists of floats

    Returns
    -------
    fwhm : List[List[float]]
        The FWHM kernels as a list of lists of floats. All inner lists are of
        length 3 as they encode each physical dimension.

    Examples
    --------
    >>> _sanitize_fwhm(3.0)
    [[3.0, 3.0, 3.0]]
    >>> _sanitize_fwhm([3.0])
    [[3.0, 3.0, 3.0]]
    >>> _sanitize_fwhm((3,))
    [[3.0, 3.0, 3.0]]
    >>> _sanitize_fwhm([3.0, 2.0])
    [[3.0, 3.0, 3.0], [2.0, 2.0, 2.0]]
    >>> _sanitize_fwhm((3.0, 2.0))
    [[3.0, 3.0, 3.0], [2.0, 2.0, 2.0]]
    >>> _sanitize_fwhm([3.0, 2.0, 1.0])
    [[3.0, 3.0, 3.0], [2.0, 2.0, 2.0], [1.0, 1.0, 1.0]]
    >>> _sanitize_fwhm([[3.0, 2.0, 1.0], [2.0, 2.0, 1.0]])
    [[3.0, 2.0, 1.0], [2.0, 2.0, 1.0]]
    """
    if isinstance(fwhm, tuple):
        fwhm = list(fwhm)
    if isinstance(fwhm, (int, float)):
        fwhm = [[float(fwhm)] * 3]
    if isinstance(fwhm, list):
        if len(fwhm) == 0:
            raise ValueError("Empty FWHM list provided.")
        if isinstance(fwhm[0], list):
            if not all([isinstance(f, list) for f in fwhm]):
                raise ValueError(
                    "Expecting a list of lists of ints or a list of ints for FWHM."
                )
            if not all([len(f) == 3 for f in fwhm]):
                raise ValueError(
                    "When providing a list of lists of ints for FWHM, all inner lists must have length 3"
                )
            for f in fwhm:
                if not all([isinstance(ff, (int, float)) for ff in f]):
                    raise ValueError(
                        "Expecting a list of lists of ints or a list of ints for FWHM."
                    )
        else:
            if all([isinstance(f, (int, float)) for f in fwhm]):
                return [[float(f)] * 3 for f in fwhm]
            else:
                raise ValueError(
                    "Expecting a list of lists of floats or a list of floats for FWHM."
                )
    return fwhm
