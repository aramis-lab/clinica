from typing import Optional

import numpy as np

__all__ = ["Results"]


def _convert_arrays_to_lists(data: dict) -> dict:
    """If the input dictionary contains numpy arrays, this function will
    cast them to lists and return the same dictionary with the lists instead
    of the numpy arrays.

    Parameters
    ----------
    data : dict
        The dictionary to clean.

    Returns
    -------
    new_data : dict
        The dictionary with arrays casted to lists.
    """
    new_data = dict()
    for k, v in data.items():
        if isinstance(v, dict):
            new_data[k] = _convert_arrays_to_lists(v)
        elif isinstance(v, np.ndarray):
            new_data[k] = v.tolist()
        else:
            new_data[k] = v
    return new_data


class Results:
    """Common class for GLM results."""

    def to_dict(self) -> dict:
        """Returns the `Results` instance in dict format.

        Private attributes and all methods are not returned.

        This function does not perform any casting.

        Returns
        -------
        data : dict
            Resulting dictionary.
        """
        import inspect

        data = dict()
        for attribute in inspect.getmembers(self):
            name, value = attribute
            if not name.startswith("_"):
                if not inspect.ismethod(value):
                    if hasattr(value, "to_dict"):
                        data[name] = value.to_dict()
                    else:
                        data[name] = value
        return data

    def to_json(self, indent: Optional[int] = 4) -> str:
        """Returns the json of the `Results` instance.

        Parameters
        ----------
        indent : int, optional
            Indent to use. Default=4.

        Returns
        -------
        str :
            The JSON dumps of the results.
        """
        import json

        return json.dumps(_convert_arrays_to_lists(self.to_dict()), indent=indent)
