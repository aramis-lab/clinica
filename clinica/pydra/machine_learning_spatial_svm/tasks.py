from typing import List

import numpy as np
from pydra.mark import annotate, task


@task
@annotate({"return": {"g_fisher_tensor": np.array, "output_fisher_tensor": str}})
def obtain_g_fisher_tensor_task(dartel_input, FWHM):
    """Pydra task for obtaining g fisher tensor.

    ..note::
        Please refer to the documentation of function
        `clinica.pipelines.machine_learning_spatial_svm.spatial_svm_utils.obtain_g_fisher_tensor`.
    """
    from clinica.pipelines.machine_learning_spatial_svm.spatial_svm_utils import (
        obtain_g_fisher_tensor,
    )

    return obtain_g_fisher_tensor(dartel_input, FWHM)


@task
@annotate({"return": {"t_step": np.array, "output_data": str}})
def obtain_time_step_estimation_task(dartel_input, FWHM, g):
    """Pydra task for computing time step estimation.

    ..note::
        Please refer to the documentation of function
        `clinica.pipelines.machine_learning_spatial_svm.spatial_svm_utils.obtain_time_step_estimation`.
    """
    from clinica.pipelines.machine_learning_spatial_svm.spatial_svm_utils import (
        obtain_time_step_estimation,
    )

    return obtain_time_step_estimation(dartel_input, FWHM, g)


@task
@annotate({"return": {"regularized": str}})
def heat_solver_equation_task(input_image, g, FWHM, t_step, dartel_input):
    """Pydra task for computing heat solver equation.

    ..note::
        Please refer to the documentation of function
        `clinica.pipelines.machine_learning_spatial_svm.spatial_svm_utils.heat_solver_equation`.
    """
    from clinica.pipelines.machine_learning_spatial_svm.spatial_svm_utils import (
        heat_solver_equation,
    )

    return heat_solver_equation(input_image, g, FWHM, t_step, dartel_input)
