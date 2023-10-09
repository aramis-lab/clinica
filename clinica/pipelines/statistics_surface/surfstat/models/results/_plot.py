from pathlib import Path
from typing import Callable

from nilearn.surface import Mesh

from ._statistics import StatisticsResults

__all__ = ["StatisticsResultsPlotter"]


class StatisticsResultsPlotter:
    """Class responsible to plotting results of GLM fit.

    Attributes
    ----------
    output_file : PathLike
        Path to the output file.

    mesh : nilearn.surface.Mesh
        The mesh to be used for plotting results.
    """

    def __init__(self, output_file: Path, mesh: Mesh):
        self.output_file = output_file
        self.mesh = mesh
        self.plotting_extension = ".png"
        self.no_plot = {"coefficients"}  # Elements which should not be plotted

    def plot(self, result: StatisticsResults, method: str) -> None:
        """Plot the results.

        Parameters
        ----------
        result : StatisticsResults
            The results to be plotted.

        method : str
            The plotting method to use.
        """
        plotter = self._get_plotter(method)
        plotter(result)

    def _get_plotter(self, method: str) -> Callable[[StatisticsResults], None]:
        """Returns the plotting method from its name.

        Parameters
        ----------
        method : str
            Name of the plotting method to use.

        Returns
        -------
        Callable :
            Plotting method.
        """
        if method == "nilearn_plot_surf_stat_map":
            return self._plot_stat_maps
        raise NotImplementedError(f"Plotting method {method} is not implemented.")

    def _plot_stat_maps(self, result: StatisticsResults) -> None:
        """Wrapper around the `nilearn.plotting.plot_surf_stat_map` method.

        Parameters
        ----------
        result : StatisticsResults
            The results to plot.
        """
        from nilearn.plotting import plot_surf_stat_map

        from clinica.utils.stream import cprint

        for name, res in result.to_dict().items():
            if name not in self.no_plot:
                texture = res
                threshold = None
                plot_filename = (
                    str(self.output_file) + "_" + name + self.plotting_extension
                )
                if isinstance(res, dict):
                    texture = res["P"]
                    threshold = res["thresh"]
                cprint(msg=f"Saving plot to {plot_filename}", lvl="info")
                plot_surf_stat_map(
                    self.mesh,
                    texture,
                    threshold=threshold,
                    output_file=plot_filename,
                    title=name,
                )
