"""Script that generate the skeleton for a new Pipeline."""
from os import getcwd

import click


@click.group("generate")
def cli():
    """Instantiate a new pipeline from available templates."""
    pass


@cli.command("template")
@click.argument("name")
@click.option(
    "-d",
    "--output-dir",
    type=click.Path(exists=True, writable=True),
    default=getcwd(),
    help="Target path to the generated template.",
)
def generate_template(name: str, output_dir: str):
    """Generate the skeleton for a new pipeline."""
    from os import mkdir
    from os.path import isdir, join, splitext

    from jinja2 import Environment, PackageLoader, select_autoescape

    pipeline = dict(
        title=name,
        module_name=name.replace(" ", "_").lower(),
        class_name=name.replace(" ", "_"),
        command_name=name.replace(" ", "-").lower(),
        dir=join(output_dir, name.replace(" ", "_").lower()),
    )

    env = Environment(
        loader=PackageLoader("clinica", "resources/templates/pipeline_template"),
        autoescape=select_autoescape(["py"]),
    )

    if not isdir(pipeline["dir"]):
        mkdir(pipeline["dir"])

    for template_file in env.list_templates(extensions="j2"):
        template = env.get_template(template_file)
        rendered_file, template_ext = splitext(template_file)
        rendered_filename, rendered_file_ext = splitext(rendered_file)
        if rendered_file_ext == ".py" and rendered_filename != "__init__":
            path_to_write = join(
                pipeline["dir"], pipeline["module_name"] + "_" + rendered_file
            )
        else:
            path_to_write = join(pipeline["dir"], rendered_file)
        with open(path_to_write, "w+") as file_to_write:
            print(f"Generating template {file_to_write.name}")
            file_to_write.write(template.render(pipeline=pipeline))


if __name__ == "__main__":
    cli()
