from pathlib import Path
from typing import Iterable, Optional

import click
from pandas import DataFrame

from clinica.iotools.converters import cli_param


def parse_description(description: str) -> Optional[str]:
    from pandas import NA

    # Normalize description.
    description = description.lower().replace("-", "")

    if "mprage" in description:
        return "T1w"
    elif "flair" in description:
        return "FLAIR"
    elif "t2" in description:
        return "T2w"
    elif "asl" in description:
        return "PDw"
    elif "pib" in description:
        return "pet"
    elif "fdg" in description:
        return "pet"
    else:
        return NA


def compute_filename(
    subject_label: str,
    session_label: str,
    datatype: str,
    suffix: str,
    trc_label: str,
    rec_label: str,
) -> str:
    from pandas import notna

    return (
        f"{datatype}/sub-{subject_label}_ses-{session_label}"
        f"{'_trc-'+trc_label if notna(trc_label) else ''}"
        f"{'_rec-'+rec_label if notna(rec_label) else ''}"
        f"_{suffix}.nii.gz"
    )


def write_tsv(dataframe: DataFrame, path: Path, *args, **kwargs) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    dataframe.to_csv(
        path, sep="\t", na_rep="n/a", date_format="%Y-%m-%d", *args, **kwargs
    )


def find_image_data_dir(sourcedata_dir: Path, image_data_id: str) -> Path:
    dirs = set(f.parent for f in sourcedata_dir.rglob(f"*_{image_data_id}.*"))
    return dirs.pop()


def do_conversion(
    image_data_dir: Path,
    bids_filename: Path,
):
    bids_filename.parent.mkdir(parents=True, exist_ok=True)

    run_dcm2niix(
        source_directory=image_data_dir,
        output_directory=bids_filename.parent,
        output_format=bids_filename.name.replace(".nii.gz", ""),
    )


def run_dcm2niix(
    source_directory: Path,
    output_directory: Path,
    output_format: str,
) -> Iterable[Path]:
    import subprocess

    subprocess.run(
        f"dcm2niix -9 -b y -ba y -f {output_format} -o {output_directory} -z i {source_directory}",
        shell=True,
    )

    return output_directory.glob(f"{output_format}.nii*")


@click.command(name="nifd-to-bids")
@cli_param.dataset_directory
@cli_param.clinical_data_directory
@cli_param.bids_directory
def cli(
    dataset_directory: str,
    clinical_data_directory: str,
    bids_directory: str,
) -> None:
    """NIFD to BIDS converter.

    Convert the imaging and clinical data of NIFD (http://4rtni-ftldni.ini.usc.edu/), located in DATASET_DIRECTORY and
    CLINICAL_DATA_DIRECTORY respectively, to a BIDS dataset in the target BIDS_DIRECTORY.
    """
    from pathlib import Path
    from shutil import copytree
    from tempfile import TemporaryDirectory

    from pandas import NA, read_csv

    from clinica.utils.check_dependency import check_dcm2niix
    from clinica.utils.stream import cprint

    check_dcm2niix()

    try:
        image_data_file = next(Path(dataset_directory).glob("*.csv"))
    except StopIteration:
        raise FileNotFoundError("Imaging collection file not found.")

    image_data = (
        read_csv(
            image_data_file,
            index_col=0,
            usecols=[
                "Image Data ID",
                "Subject",
                "Group",
                "Sex",
                "Age",
                "Visit",
                "Modality",
                "Description",
                "Acq Date",
            ],
            parse_dates=["Acq Date"],
        )
        .rename(columns=lambda x: x.lower().replace(" ", "_"))
        .rename_axis(index=lambda x: x.lower().replace(" ", "_"))
        .convert_dtypes()
    )

    image_data = (
        image_data.assign(
            subject_label=lambda d: d.subject.apply(
                lambda x: f"NIFD{x.replace('_', '')}"
            ),
            session_label=lambda d: d.visit.apply(lambda x: f"M{(6 * (x - 1)):02d}"),
        )
        .drop(columns=["subject", "visit"])
        .sort_values(by=["subject_label", "session_label"])
    )

    image_data = (
        image_data.assign(
            datatype=lambda d: d.modality.replace({"PET": "pet", "MRI": "anat"}),
            suffix=lambda d: d.description.apply(parse_description),
            trc_label=lambda d: d.description.apply(
                lambda x: "11CPIB" if ":PIB" in x else ("18FFDG" if ":FDG" in x else NA)
            ),
            rec_label=lambda d: d.description.apply(
                lambda x: "IR" if ":IR" in x else ("RP" if ":RP" in x else NA)
            ),
        )
        .drop(columns=["description"])
        .dropna(subset=["suffix"])
    )

    image_data = image_data.assign(
        filename=lambda d: d.apply(
            lambda x: compute_filename(
                subject_label=x.subject_label,
                session_label=x.session_label,
                datatype=x.datatype,
                suffix=x.suffix,
                trc_label=x.trc_label,
                rec_label=x.rec_label,
            ),
            axis=1,
        )
    ).drop(columns=["datatype", "suffix", "trc_label", "rec_label"])

    with TemporaryDirectory() as tmpdir:
        staging_directory = Path(tmpdir)

        subject_data = (
            image_data[["subject_label", "group", "sex", "age"]]
            .drop_duplicates(keep="first")
            .assign(participant_id=lambda d: "sub-" + d.subject_label)
            .drop(columns=["subject_label"])
            .set_index("participant_id")
        )

        participants_file = staging_directory / "participants.tsv"

        write_tsv(dataframe=subject_data, path=participants_file)

        for subject_label, per_subject in image_data.groupby(by="subject_label"):
            session_data = (
                per_subject[["session_label", "acq_date"]]
                .drop_duplicates(keep="first")
                .assign(session_id=lambda d: "ses-" + d.session_label)
                .drop(columns=["session_label"])
                .set_index("session_id")
            )

            sessions_file = (
                staging_directory
                / f"sub-{subject_label}"
                / f"sub-{subject_label}_sessions.tsv"
            )

            write_tsv(dataframe=session_data, path=sessions_file)

            for session_label, per_session in per_subject.groupby(by="session_label"):
                scans_data = (
                    per_session[["filename", "acq_date"]]
                    .drop_duplicates(keep="last")
                    .set_index("filename")
                )

                scans_file = (
                    staging_directory
                    / f"sub-{subject_label}"
                    / f"ses-{session_label}"
                    / f"sub-{subject_label}_ses-{session_label}_scans.tsv"
                )

                write_tsv(dataframe=scans_data, path=scans_file)

                for image_data_id, filename in per_session.filename.items():
                    bids_filename = (
                        staging_directory
                        / f"sub-{subject_label}"
                        / f"ses-{session_label}"
                        / Path(filename)
                    )

                    image_data_dir = find_image_data_dir(
                        sourcedata_dir=Path(dataset_directory),
                        image_data_id=image_data_id,
                    )

                    do_conversion(
                        image_data_dir=image_data_dir,
                        bids_filename=bids_filename,
                    )

        copytree(staging_directory, bids_directory, dirs_exist_ok=True)

    cprint("Conversion to BIDS succeeded.")


if __name__ == "__main__":
    cli()
