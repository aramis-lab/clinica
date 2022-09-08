from enum import Enum
from typing import IO

from attrs import define, fields
from cattr.gen import make_dict_unstructure_fn, override
from cattr.preconf.json import make_converter

BIDS_VERSION = "1.7.0"


@define
class BIDSReadme:
    """Model representing a BIDS ReadMe.

    See
    """

    class DatasetType(str, Enum):
        raw = "raw"
        derivative = "derivative"

    name: str
    bids_version: str = BIDS_VERSION
    dataset_type: DatasetType = DatasetType.raw

    def write(self, to: IO[str], readme_dict):
        import clinica

        # datadict_link = {
        #     "OASIS-1": "https://www.oasis-brains.org/#access",
        #     "OASIS-3": "https://www.oasis-brains.org/#access",
        #     "ADNI": "http://adni.loni.usc.edu",
        #     "AIBL": "http://adni.loni.usc.edu/study-design/collaborative-studies/aibl/",
        #     "UKB": "https://www.ukbiobank.ac.uk/",
        #     "NIFD": "https://ida.loni.usc.edu/home/projectPage.jsp?project=NIFD&page=HOME&subPage=OVERVIEW_PR#",
        #     "HABS": "https://habs.mgh.harvard.edu",
        # }
        # datadict_description = {
        #     "OASIS-3": "OASIS-3 is a retrospective compilation of data for 1378 participants that were collected across several ongoing projects through the WUSTL Knight ADRC over the course of 30years. Participants include 755 cognitively normal adults and 622 individuals at various stages of cognitive decline ranging in age from 42-95yrs. All participants were assigned a new random identifier and all dates were removed and normalized to reflect days from entry into study. The dataset contains 2842 MR sessions which include T1w, T2w, FLAIR, ASL, SWI, time of flight, resting-state BOLD, and DTI sequences. Many of the MR sessions are accompanied by volumetric segmentation files produced through FreeSurfer processing. PET imaging from different tracers, PIB, AV45, and FDG, totaling over 2157 raw imaging scans and the accompanying post-processed files from the Pet Unified Pipeline (PUP) are also available in OASIS-3. ",
        #     "OASIS-1": "This set consists of a cross-sectional collection of 416 subjects aged 18 to 96. For each subject, 3 or 4 individual T1-weighted MRI scans obtained in single scan sessions are included. The subjects are all right-handed and include both men and women. 100 of the included subjects over the age of 60 have been clinically diagnosed with very mild to moderate Alzheimer’s disease (AD). Additionally, a reliability data set is included containing 20 nondemented subjects imaged on a subsequent visit within 90 days of their initial session.",
        #     "ADNI": "ADNI is a global research effort that actively supports the investigation and development of treatments that slow or stop the progression of Alzheimer's disease (AD).This multisite, longitudinal study assesses clinical, imaging, genetic and biospecimen biomarkers through the process of normal aging to mild cognitive impairment (MCI) and AD dementia.With established, standardized methods for imaging and biomarker collection and analysis, ADNI facilitates a way for scientists to conduct cohesive research and share compatible data with other researchers around the world.",
        #     "AIBL": "The Australian Imaging, Biomarker & Lifestyle Flagship Study of Ageing (AIBL) seeks to discover which biomarkers, cognitive characteristics, and health and lifestyle factors determine the development of AD.Although AIBL and ADNI have many of the same goals, there are differences between the two projects.",
        #     "UKB": "UK Biobank is a large-scale biomedical database and research resource, containing in-depth genetic and health information from half a million UK participants. The database is regularly augmented with additional data and is globally accessible to approved researchers undertaking vital research into the most common and life-threatening diseases. It is a major contributor to the advancement of modern medicine it and has led to the discovery of several scientific advances and numerous treatments to improve human health.",
        #     "NIFD": "NIFD is the nickname for the frontotemporal lobar degeneration neuroimaging initiative (FTLDNI, AG032306), which was funded by the NIA and NINDS to characterize longitudinal clinical and imaging changes in FTLD.The imaging and clinical methods are the same for NIFD and for the 4-Repeat Tauopathy Neuroimaging Initiative (4RTNI), which is also available for download from LONI. Controls for NIFD are the same controls as those collected for 4RTNI.",
        #     "HABS": "The overall goal of the Harvard Aging Brain Study (HABS) is to elucidate the earliest changes in molecular, functional and structural imaging markers that signal the transition from normal cognition to progressive cognitive decline along the trajectory of preclinical Alzheimer’s Disease.",
        # }

        to.write(
            f"This BIDS directory was generated with Clinica v{clinica.__version__}.\n"
            f"More information on https://www.clinica.run\n"
            f"\n"
            f"Study: {self.name}\n"
            f"\n"
            f"{readme_dict['desc']}\n\n"
            f"Find more about it and about the data user agreement: {readme_dict['link']}"
        )


def _rename(name: str) -> str:
    """Rename attributes following the specification for the JSON file.

    Basically pascal case with known acronyms such as BIDS fully capitalized.
    """
    return "".join(
        word.upper() if word == "bids" else word.capitalize()
        for word in name.split("_")
    )


# Register a JSON converter for the BIDS dataset description model.
converter = make_converter()

converter.register_unstructure_hook(
    BIDSReadme,
    make_dict_unstructure_fn(
        BIDSReadme,
        converter,
        **{a.name: override(rename=_rename(a.name)) for a in fields(BIDSReadme)},
    ),
)
