# coding: utf8

import clinica.engine as ce


class T1VolumeExistingTemplateCLI(ce.CmdParser):
    def define_name(self):
        """Define the sub-command name to run this pipeline."""
        self._name = "t1-volume-existing-template"

    def define_description(self):
        """Define a description of this pipeline."""
        self._description = (
            "Volume-based processing of T1-weighted MR images using an existing DARTEL template:\n"
            "https://aramislab.paris.inria.fr/clinica/docs/public/latest/Pipelines/T1_Volume/"
        )

    def define_options(self):
        """Define the sub-command arguments."""
        from clinica.engine.cmdparser import PIPELINE_CATEGORIES

        # Clinica compulsory arguments (e.g. BIDS, CAPS, group_label)
        clinica_comp = self._args.add_argument_group(
            PIPELINE_CATEGORIES["CLINICA_COMPULSORY"]
        )
        clinica_comp.add_argument("bids_directory", help="Path to the BIDS directory.")
        clinica_comp.add_argument("caps_directory", help="Path to the CAPS directory.")
        clinica_comp.add_argument(
            "group_label",
            help="User-defined identifier for the provided group of subjects.",
        )
        # Optional arguments (e.g. FWHM)
        optional = self._args.add_argument_group(PIPELINE_CATEGORIES["OPTIONAL"])
        optional.add_argument(
            "-s",
            "--smooth",
            nargs="+",
            type=int,
            default=[8],
            help="A list of integers specifying the different isomorphic FWHM in millimeters "
            "to smooth the image (default: --smooth 8).",
        )
        # Clinica standard arguments (e.g. --n_procs)
        self.add_clinica_standard_arguments()
        # Advanced arguments (i.e. tricky parameters)
        advanced = self._args.add_argument_group(PIPELINE_CATEGORIES["ADVANCED"])
        # t1-volume-tissue-segmentation
        advanced.add_argument(
            "-tc",
            "--tissue_classes",
            metavar="",
            nargs="+",
            type=int,
            default=[1, 2, 3],
            choices=range(1, 7),
            help="Tissue classes (1: gray matter (GM), 2: white matter (WM), "
            "3: cerebrospinal fluid (CSF), 4: bone, 5: soft-tissue, 6: background) to save "
            "(default: GM, WM and CSF i.e. --tissue_classes 1 2 3).",
        )
        advanced.add_argument(
            "-tpm",
            "--tissue_probability_maps",
            metavar="TissueProbabilityMap.nii",
            default=None,
            help="Tissue probability maps to use for segmentation "
            "(default: TPM.nii from SPM software).",
        )
        advanced.add_argument(
            "-dswu",
            "--dont_save_warped_unmodulated",
            action="store_true",
            default=False,
            help="Do not save warped unmodulated images for tissues specified "
            "in --tissue_classes flag.",
        )
        advanced.add_argument(
            "-swm",
            "--save_warped_modulated",
            action="store_true",
            default=False,
            help="Save warped modulated images for tissues specified in --tissue_classes flag.",
        )
        # t1-volume-tissue-segmentation / t1-volume-create-dartel
        advanced.add_argument(
            "-dt",
            "--dartel_tissues",
            metavar="",
            nargs="+",
            type=int,
            default=[1, 2, 3],
            choices=range(1, 7),
            help="Tissues to use for DARTEL template calculation "
            "(default: GM, WM and CSF i.e. --dartel_tissues 1 2 3).",
        )
        # t1-volume-dartel2mni
        advanced.add_argument(
            "-t",
            "--tissues",
            metavar="",
            nargs="+",
            type=int,
            default=[1, 2, 3],
            choices=range(1, 7),
            help="Tissues to create flow fields to DARTEL template "
            "(default: GM, WM and CSF i.e. --tissues 1 2 3).",
        )
        advanced.add_argument(
            "-m",
            "--modulate",
            type=bool,
            default=True,
            metavar=("True/False"),
            help="A boolean. Modulate output images - no modulation preserves concentrations "
            "(default: --modulate True).",
        )
        advanced.add_argument(
            "-vs",
            "--voxel_size",
            metavar=("float"),
            nargs=3,
            type=float,
            help="A list of 3 floats specifying the voxel sizeof the output image "
            "(default: --voxel_size 1.5 1.5 1.5).",
        )

    def run_command(self, args):
        """Run the pipeline with defined args."""
        from colorama import Fore

        from clinica.utils.stream import cprint

        from ..t1_volume_dartel2mni.t1_volume_dartel2mni_cli import (
            T1VolumeDartel2MNICLI,
        )
        from ..t1_volume_parcellation.t1_volume_parcellation_cli import (
            T1VolumeParcellationCLI,
        )
        from ..t1_volume_register_dartel.t1_volume_register_dartel_cli import (
            T1VolumeRegisterDartelCLI,
        )
        from ..t1_volume_tissue_segmentation.t1_volume_tissue_segmentation_cli import (
            T1VolumeTissueSegmentationCLI,
        )

        cprint(
            f"The t1-volume-existing-template pipeline is divided into 4 parts:\n"
            f"\t{Fore.BLUE}t1-volume-tissue-segmentation pipeline{Fore.RESET}: "
            f"Tissue segmentation, bias correction and spatial normalization to MNI space\n"
            f"\t{Fore.BLUE}t1-volume-register-dartel pipeline{Fore.RESET}: "
            f"Inter-subject registration using an existing DARTEL template\n"
            f"\t{Fore.BLUE}t1-volume-dartel2mni pipeline{Fore.RESET}: "
            f"DARTEL template to MNI\n"
            f"\t{Fore.BLUE}t1-volume-parcellation pipeline{Fore.RESET}: "
            f"Atlas statistics"
        )

        cprint(
            f"{Fore.BLUE}\nPart 1/4: Running t1-volume-segmentation pipeline{Fore.RESET}"
        )
        tissue_segmentation_cli = T1VolumeTissueSegmentationCLI()
        tissue_segmentation_cli.run_command(args)

        cprint(
            f"{Fore.BLUE}\nPart 2/4: Running t1-volume-register-dartel pipeline{Fore.RESET}"
        )
        register_dartel_cli = T1VolumeRegisterDartelCLI()
        register_dartel_cli.run_command(args)

        cprint(
            f"{Fore.BLUE}\nPart 3/4: Running t1-volume-dartel2mni pipeline{Fore.RESET}"
        )
        dartel2mni_cli = T1VolumeDartel2MNICLI()
        dartel2mni_cli.run_command(args)

        cprint(
            f"{Fore.BLUE}\nPart 4/4: Running t1-volume-parcellation pipeline{Fore.RESET}"
        )
        parcellation_cli = T1VolumeParcellationCLI()
        parcellation_cli.run_command(args)
