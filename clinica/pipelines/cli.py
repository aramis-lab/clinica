import click


class RegistrationOrderGroup(click.Group):
    """CLI group which lists commands by order or registration."""

    def list_commands(self, ctx):
        return self.commands.keys()


@click.group(cls=RegistrationOrderGroup, name="run")
def cli() -> None:
    """Run pipelines on BIDS and CAPS datasets."""
    pass


if __name__ == "__main__":
    cli()
