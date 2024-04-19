from __future__ import annotations

from pathlib import PurePath

from pydra.mark import annotate, task


def download_file(url: str, to: str) -> PurePath:
    from shutil import copyfileobj
    from ssl import SSLContext
    from urllib.request import urlopen

    print(f"Downloading {url} to {to}...")

    response = urlopen(url=url, context=SSLContext())
    with open(to, mode="wb") as f:
        copyfileobj(response, f)

    return PurePath(to)


@task
@annotate({"return": {"mni_template_file": PurePath}})
def download_mni_template_2009a() -> PurePath:
    return _download_mni_template(version="09a")


@task
@annotate({"return": {"mni_template_file": PurePath}})
def download_mni_template_2009c() -> PurePath:
    return _download_mni_template(version="09c")


def _download_mni_template(version: str) -> PurePath:
    from pathlib import Path

    return download_file(
        url=f"https://aramislab.paris.inria.fr/files/data/img_t1_linear/mni_icbm152_t1_tal_nlin_sym_{version}.nii.gz",
        to=str(Path.cwd() / f"mni_icbm152_t1_tal_nlin_sym_{version}.nii.gz"),
    )


@task
@annotate({"return": {"ref_template_file": PurePath}})
def download_ref_template() -> PurePath:
    from pathlib import Path

    return download_file(
        url="https://aramislab.paris.inria.fr/files/data/img_t1_linear/ref_cropped_template.nii.gz",
        to=str(Path.cwd() / "ref_cropped_template.nii.gz"),
    )
