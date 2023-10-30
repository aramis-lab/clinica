"""Common CLI options used by Clinica tools."""

from multiprocessing import cpu_count

import click
from click_option_group import optgroup

option = optgroup.option

global_option_group = optgroup.group(
    "Options common to all clinica tools", help="Options common to all clinica tools"
)


def validate_n_procs(ctx, param, value):
    if isinstance(value, int) and 1 <= value <= cpu_count():
        return value
    raise click.BadParameter(
        "The number of processes should be between 1 and the number "
        f"of CPU available on your system ({cpu_count()}). "
        f"A value of {value} was provided instead. "
    )


n_procs = option(
    "-np",
    "--n_procs",
    type=int,
    default=lambda: max(cpu_count() - 1, 1),
    show_default="Number of available CPU minus one",
    callback=validate_n_procs,
    help="Number of cores used to run in parallel.",
)
