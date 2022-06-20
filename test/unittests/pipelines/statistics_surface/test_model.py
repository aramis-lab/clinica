
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


@pytest.mark.parametrize("contrast", ["age", "-age"])
def test_correlation_glm_instantiation(df, contrast):
    from clinica.pipelines.statistics_surface._model import CorrelationGLM
    design = "1 + age"
    parameters = {
        "sizeoffwhm": 2.0,
        "group_label": "group_label",
    }
    model = CorrelationGLM(design, df, "feature_label", contrast, **parameters)
    assert not model.with_interaction
    assert model.group_label == "group_label"
    assert model.feature_label == "feature_label"
    assert model.fwhm == 2.0
    sign = "positive" if contrast == "age" else "negative"
    assert model.contrast_sign == sign
    assert model.absolute_contrast_name == "age"
    assert isinstance(model.contrasts, dict)
    assert len(model.contrasts) == 1
    mult = 1 if sign == "positive" else -1
    assert_array_equal(
        model.contrasts[contrast].values,
        mult * np.array([78., 73.4, 70.8, 82.3, 60.6, 72.1, 74.2])
    )
    with pytest.raises(
        ValueError,
        match="Unknown contrast foo"
    ):
        model.filename_root("foo")
    expected = f"group-group_label_correlation-age-{sign}_measure-feature_label_fwhm-2.0"
    assert model.filename_root(contrast) == expected


def test_group_glm_instantiation(df):
    from clinica.pipelines.statistics_surface._model import GroupGLM
    design = "1 + age"
    parameters = {
        "sizeoffwhm": 2.0,
        "group_label": "group_label",
    }
    with pytest.raises(
        ValueError,
        match="Contrast should refer to a categorical variable for group comparison."
    ):
        GroupGLM(design, df, "feature_label", "age", **parameters)
    model = GroupGLM(design, df, "feature_label", "sex", **parameters)
    assert not model.with_interaction
    assert model.group_label == "group_label"
    assert model.fwhm == 2.0
    assert isinstance(model.contrasts, dict)
    contrast_names = ['Female-lt-Male', 'Male-lt-Female']
    assert set(model.contrasts.keys()) == set(contrast_names)
    for contrast_name, sign in zip(contrast_names, [1, -1]):
        assert_array_equal(
            model.contrasts[contrast_name].values,
            sign * np.array([1, -1, 1, 1, 1, -1, 1])
        )
    with pytest.raises(
            ValueError,
            match="Unknown contrast foo"
    ):
        model.filename_root("foo")
    for contrast_name in contrast_names:
        expected = f'group-group_label_{contrast_name}_measure-feature_label_fwhm-2.0'
        assert model.filename_root(contrast_name) == expected


def test_group_glm_with_interaction_instantiation(df):
    from clinica.pipelines.statistics_surface._model import GroupGLMWithInteraction
    design = "1 + age"
    parameters = {
        "sizeoffwhm": 2.0,
        "group_label": "group_label",
    }
    with pytest.raises(
            ValueError,
            match=(
                "The contrast must be an interaction between one continuous "
                "variable and one categorical variable."
            )
    ):
        GroupGLMWithInteraction(design, df, "feature_label", "age", **parameters)
    model = GroupGLMWithInteraction(design, df, "feature_label", "age * sex", **parameters)
    assert model.with_interaction
    assert model.group_label == "group_label"
    assert model.fwhm == 2.0
    assert isinstance(model.contrasts, dict)
    assert len(model.contrasts) == 1
    assert_array_equal(
        model.contrasts["age * sex"].values,
        np.array([78., -73.4, 70.8, 82.3, 60.6, -72.1, 74.2])
    )
    with pytest.raises(
            ValueError,
            match="Unknown contrast foo"
    ):
        model.filename_root("foo")
    assert model.filename_root("age * sex") == "interaction-age * sex_measure-feature_label_fwhm-2.0"
