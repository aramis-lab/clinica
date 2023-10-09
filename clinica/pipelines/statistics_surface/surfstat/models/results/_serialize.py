from pathlib import Path
from typing import Callable

from ._statistics import StatisticsResults

__all__ = ["StatisticsResultsSerializer"]


class StatisticsResultsSerializer:
    """This class is responsible for writing instances of `StatisticsResults`
    to disk through different methods.

    Attributes
    ----------
    output_file : PathLike
        Path and filename root to be used.
    """

    def __init__(self, output_file: Path):
        self.output_file = output_file
        self.json_extension = "_results.json"
        self.json_indent = 4
        self.mat_extension = ".mat"

    def save(self, result: StatisticsResults, method: str) -> None:
        """Save provided `StatisticsResults` to disk with provided method.

        Parameters
        ----------
        result : StatisticsResults
            Results to be saved.

        method : str
            Name of the saving method to use.
        """
        writer = self._get_writer(method)
        writer(result)

    def _get_writer(self, method: str) -> Callable[[StatisticsResults], None]:
        """Returns a writer method from its name.

        Parameters
        ----------
        method : str
            The name of the writing method to use.

        Returns
        -------
        Callable :
            The writing method.
        """
        if method.lower() == "json":
            return self._write_to_json
        elif method.lower() == "mat":
            return self._write_to_mat
        raise NotImplementedError(f"Serializing method {method} is not implemented.")

    def _write_to_json(self, results: StatisticsResults) -> None:
        """Write the provided `StatisticsResults` to JSON format.

        Parameters
        ----------
        results : StatisticsResults
            The results to write to disk in JSON format.
        """
        import json
        import os

        from clinica.utils.stream import cprint

        out_json_file = Path(str(self.output_file) + self.json_extension)
        if not os.path.exists(out_json_file.parents[0]):
            os.makedirs(out_json_file.parents[0])
        cprint(
            msg=f"Writing results to JSON in {out_json_file}...",
            lvl="info",
        )
        with open(out_json_file, "w") as fp:
            json.dump(results.to_json(indent=self.json_indent), fp)

    def _write_to_mat(self, results: StatisticsResults) -> None:
        """Write the provided `StatisticsResults` to MAT format.

        Parameters
        ----------
        results : StatisticsResults
            The results to write to disk in MAT format.
        """
        from scipy.io import savemat

        from clinica.utils.stream import cprint

        # These labels are used for compatibility with the previous
        # MATLAB implementation of the Statistics Surface Pipeline
        # of Clinica.
        struct_labels = {
            "coefficients": "coef",
            "TStatistics": "tvaluewithmask",
            "uncorrectedPValue": "uncorrectedpvaluesstruct",
            "correctedPValue": "correctedpvaluesstruct",
            "FDR": "FDR",
        }
        for name, res in results.to_dict().items():
            if name in struct_labels:
                mat_filename = str(self.output_file) + "_" + name + self.mat_extension
                cprint(
                    msg=f"Writing {name} results to MAT in  {mat_filename}",
                    lvl="info",
                )
                savemat(mat_filename, {struct_labels[name]: res})
