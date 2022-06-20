
import os
import pytest
import numpy as np
import pandas as pd
from pathlib import Path
from numpy.testing import assert_array_equal
from brainstat.stats.terms import FixedEffect

CURRENT_DIR = os.path.dirname(os.path.realpath(__file__))


@pytest.fixture
def df():
    return pd.read_csv(Path(CURRENT_DIR) / "data/subjects.tsv", sep="\t")


def test_is_categorical(df):
    from clinica.pipelines.statistics_surface._model import _is_categorical
    with pytest.raises(
        ValueError,
        match=(
            "Term foo from the design matrix is not in the columns of the "
            "provided TSV file. Please make sure that there is no typo"
        )
    ):
        _is_categorical(df, "foo")
    assert _is_categorical(df, "sex")
    assert not _is_categorical(df, "age")


def test_build_model_term_error(df):
    from clinica.pipelines.statistics_surface._model import _build_model_term
    with pytest.raises(
            ValueError,
            match=(
                    "Term foo from the design matrix is not in the columns of the "
                    "provided TSV file. Please make sure that there is no typo"
            )
    ):
        _build_model_term("foo", df)


@pytest.mark.parametrize(
        "design",
        ["1 + age", "1+age", "age +1", "age"])
def test_build_model_intercept(design, df):
    """Test that we get the same results with equivalent designs.
    Especially, the fact that adding explicitly the intercept doesn't change the results.
    Test also that spaces in the design expression have no effect.
    """
    from clinica.pipelines.statistics_surface._model import _build_model
    model = _build_model(design, df)
    assert isinstance(model, FixedEffect)
    assert len(model.m.columns) == 2
    assert_array_equal(
        model.intercept,
        np.array([1, 1, 1, 1, 1, 1, 1])
    )
    assert_array_equal(
        model.age,
        np.array([78., 73.4, 70.8, 82.3, 60.6, 72.1, 74.2])
    )


def test_build_model(df):
    from clinica.pipelines.statistics_surface._model import _build_model
    model = _build_model("1 + age + sex", df)
    assert isinstance(model, FixedEffect)
    assert len(model.m.columns) == 4
    assert_array_equal(
        model.intercept,
        np.array([1, 1, 1, 1, 1, 1, 1])
    )
    assert_array_equal(
        model.age,
        np.array([78., 73.4, 70.8, 82.3, 60.6, 72.1, 74.2])
    )
    assert_array_equal(
        model.sex_Female,
        np.array([1, 0, 1, 1, 1, 0, 1])
    )
    assert_array_equal(
        model.sex_Male,
        np.array([0, 1, 0, 0, 0, 1, 0])
    )
    model = _build_model("1 + age + sex + age * sex", df)
    assert isinstance(model, FixedEffect)
    assert len(model.m.columns) == 6
    assert_array_equal(
        model.intercept,
        np.array([1, 1, 1, 1, 1, 1, 1])
    )
    assert_array_equal(
        model.age,
        np.array([78., 73.4, 70.8, 82.3, 60.6, 72.1, 74.2])
    )
    assert_array_equal(
        model.sex_Female,
        np.array([1, 0, 1, 1, 1, 0, 1])
    )
    assert_array_equal(
        model.sex_Male,
        np.array([0, 1, 0, 0, 0, 1, 0])
    )
    assert_array_equal(
        getattr(model, "age*sex_Female"),
        np.array([78.,  0., 70.8, 82.3, 60.6,  0., 74.2])
    )
    assert_array_equal(
        getattr(model, "age*sex_Male"),
        np.array([0., 73.4,  0.,  0.,  0., 72.1,  0.])
    )


@pytest.mark.parametrize("parameters",
                         [
                             {"sizeoffwhm": 3.0},
                             {
                                 "thresholduncorrectedpvalue": 0.6,
                                 "thresholdcorrectedpvalue": 0.8,
                                 "clusterthreshold": 0.44,
                                 "sizeoffwhm": 2.66,
                             },
                         ]
)
def test_glm_instantiation(df, parameters):
    from clinica.pipelines.statistics_surface._model import(
        GLM,
        DEFAULT_THRESHOLD_UNCORRECTED_P_VALUE,
        DEFAULT_CLUSTER_THRESHOLD,
        DEFAULT_THRESHOLD_CORRECTED_P_VALUE,
    )
    from clinica.pipelines.statistics_surface.clinica_surfstat import DEFAULT_FWHM
    design = "1 + age"
    model = GLM(design, df, "feature_label", "age", **parameters)
    assert model._two_tailed
    assert model._correction == ["fdr", "rft"]
    assert model.feature_label == "feature_label"
    assert model.fwhm == parameters["sizeoffwhm"]
    threshold_uncorrected_pvalue = parameters.pop(
        "thresholduncorrectedpvalue", DEFAULT_THRESHOLD_UNCORRECTED_P_VALUE
    )
    assert model.threshold_uncorrected_pvalue == threshold_uncorrected_pvalue
    threshold_corrected_pvalue = parameters.pop(
        "thresholdcorrectedpvalue", DEFAULT_THRESHOLD_CORRECTED_P_VALUE
    )
    assert model.threshold_corrected_pvalue == threshold_corrected_pvalue
    cluster_threshold = parameters.pop(
        "clusterthreshold", DEFAULT_CLUSTER_THRESHOLD
    )
    assert model.cluster_threshold == cluster_threshold
    assert model.contrasts == dict()
    assert model.filenames == dict()
    assert model.contrast_names == list()
    assert isinstance(model.model, FixedEffect)

