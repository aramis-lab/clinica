import os
from pathlib import Path

import click

nipype_pipelines_folder = Path(os.path.dirname(__file__))
pydra_pipelines_folder = Path(os.path.dirname(__file__)).parent / "pydra"
pipeline_folders = [
    nipype_pipelines_folder,
    pydra_pipelines_folder,
]


class PipelineCommands(click.MultiCommand):
    """CLI MultiCommand for executing pipelines."""

    def list_commands(self, ctx):
        ns = {}
        rv = []
        for folder in pipeline_folders:
            for filename in folder.glob("**/*_cli.py"):
                with open(filename, "r") as fn:
                    code = compile(fn.read(), filename, "exec")
                    eval(code, ns, ns)
                rv.append(ns["pipeline_name"])
        rv.sort()
        return rv

    def get_command(self, ctx, name):
        ns = {}
        for folder in pipeline_folders:
            for filename in folder.glob("**/*_cli.py"):
                with open(filename, "r") as fn:
                    code = compile(fn.read(), filename, "exec")
                    eval(code, ns, ns)
                if ns["pipeline_name"] == name:
                    return ns["cli"]


@click.command(cls=PipelineCommands, name="run")
def cli() -> None:
    """Run pipelines on BIDS and CAPS datasets."""
    pass


if __name__ == "__main__":
    cli()
