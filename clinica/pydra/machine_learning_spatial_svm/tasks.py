from typing import List

import numpy as np
from pydra.mark import annotate, task


@task
@annotate({"return": {"g_fisher_tensor": np.array, "output_fisher_tensor": str}})
def obtain_g_fisher_tensor_task(dartel_input, FWHM):
    from clinica.pipelines.machine_learning_spatial_svm.spatial_svm_utils import (
        obtain_g_fisher_tensor,
    )

    return obtain_g_fisher_tensor(dartel_input, FWHM)


@task
@annotate({"return": {"t_step": np.array, "output_data": str}})
def obtain_time_step_estimation_task(dartel_input, FWHM, g):
    from clinica.pipelines.machine_learning_spatial_svm.spatial_svm_utils import (
        obtain_time_step_estimation,
    )

    return obtain_time_step_estimation(dartel_input, FWHM, g)


@task
@annotate({"return": {"regularized": str}})
def heat_solver_equation_task(input_image, g, FWHM, t_step, dartel_input):
    from clinica.pipelines.machine_learning_spatial_svm.spatial_svm_utils import (
        heat_solver_equation,
    )

    return heat_solver_equation(input_image, g, FWHM, t_step, dartel_input)
