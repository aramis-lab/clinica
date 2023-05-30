import os
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from brainstat.stats.terms import FixedEffect
from numpy.testing import assert_array_almost_equal, assert_array_equal

CURRENT_DIR = os.path.dirname(os.path.realpath(__file__))


@pytest.fixture
def df():
    return pd.read_csv(Path(CURRENT_DIR) / "data/subjects.tsv", sep="\t")


def test_missing_column_error(df):
    from clinica.pipelines.statistics_surface._model import _check_column_in_df

    with pytest.raises(
        ValueError,
        match=(
            "Term foo from the design matrix is not in the columns of the "
            "provided TSV file. Please make sure that there is no typo"
        ),
    ):
        _check_column_in_df(df, "foo")


def test_is_categorical(df):
    from clinica.pipelines.statistics_surface._model import _categorical_column

    assert _categorical_column(df, "sex")
    assert not _categorical_column(df, "age")


def test_build_model_term_error(df):
    from clinica.pipelines.statistics_surface._model import _build_model_term

    assert isinstance(_build_model_term("sex", df), FixedEffect)


@pytest.mark.parametrize("design", ["1 + age", "1+age", "age +1", "age"])
def test_build_model_intercept(design, df):
    """Test that we get the same results with equivalent designs.
    Especially, the fact that adding explicitly the intercept doesn't change the results.
    Test also that spaces in the design expression have no effect.
    """
    from clinica.pipelines.statistics_surface._model import _build_model

    model = _build_model(design, df)
    assert isinstance(model, FixedEffect)
    assert len(model.m.columns) == 2
    assert_array_equal(model.intercept, np.array([1, 1, 1, 1, 1, 1, 1]))
    assert_array_equal(model.age, np.array([78.0, 73.4, 70.8, 82.3, 60.6, 72.1, 74.2]))


def test_build_model(df):
    from clinica.pipelines.statistics_surface._model import _build_model

    model = _build_model("1 + age + sex", df)
    assert isinstance(model, FixedEffect)
    assert len(model.m.columns) == 4
    assert_array_equal(model.intercept, np.array([1, 1, 1, 1, 1, 1, 1]))
    assert_array_equal(model.age, np.array([78.0, 73.4, 70.8, 82.3, 60.6, 72.1, 74.2]))
    assert_array_equal(model.sex_Female, np.array([1, 0, 1, 1, 1, 0, 1]))
    assert_array_equal(model.sex_Male, np.array([0, 1, 0, 0, 0, 1, 0]))
    model = _build_model("1 + age + sex + age * sex", df)
    assert isinstance(model, FixedEffect)
    assert len(model.m.columns) == 6
    assert_array_equal(model.intercept, np.array([1, 1, 1, 1, 1, 1, 1]))
    assert_array_equal(model.age, np.array([78.0, 73.4, 70.8, 82.3, 60.6, 72.1, 74.2]))
    assert_array_equal(model.sex_Female, np.array([1, 0, 1, 1, 1, 0, 1]))
    assert_array_equal(model.sex_Male, np.array([0, 1, 0, 0, 0, 1, 0]))
    assert_array_equal(
        getattr(model, "age*sex_Female"),
        np.array([78.0, 0.0, 70.8, 82.3, 60.6, 0.0, 74.2]),
    )
    assert_array_equal(
        getattr(model, "age*sex_Male"), np.array([0.0, 73.4, 0.0, 0.0, 0.0, 72.1, 0.0])
    )


@pytest.mark.parametrize(
    "parameters",
    [
        {"sizeoffwhm": 3.0},
        {
            "thresholduncorrectedpvalue": 0.6,
            "thresholdcorrectedpvalue": 0.8,
            "clusterthreshold": 0.44,
            "sizeoffwhm": 2.66,
        },
    ],
)
def test_glm_instantiation(df, parameters):
    from clinica.pipelines.statistics_surface._model import GLM

    model = GLM("1 + age", df, "feature_label", "age")
    assert not model._two_tailed
    assert model._correction == ["fdr", "rft"]
    assert model.feature_label == "feature_label"
    assert model.fwhm == 20
    assert model.threshold_uncorrected_pvalue == 0.001
    assert model.threshold_corrected_pvalue == 0.05
    assert model.cluster_threshold == 0.001
    assert model.contrasts == dict()
    assert model.filenames == dict()
    assert model.contrast_names == list()
    assert isinstance(model.model, FixedEffect)


@pytest.mark.parametrize("contrast", ["age", "-age"])
def test_correlation_glm_instantiation(df, contrast):
    from clinica.pipelines.statistics_surface._model import CorrelationGLM

    model = CorrelationGLM("1 + age", df, "feature_label", contrast, "group_label")
    assert not model.with_interaction
    assert model.group_label == "group_label"
    assert model.feature_label == "feature_label"
    assert model.fwhm == 20
    sign = "positive" if contrast == "age" else "negative"
    assert model.contrast_sign == sign
    assert model.absolute_contrast_name == "age"
    assert isinstance(model.contrasts, dict)
    assert len(model.contrasts) == 1
    mult = 1 if sign == "positive" else -1
    assert_array_equal(
        model.contrasts[contrast].values,
        mult * np.array([78.0, 73.4, 70.8, 82.3, 60.6, 72.1, 74.2]),
    )
    with pytest.raises(ValueError, match="Unknown contrast foo"):
        model.filename_root("foo")
    expected = f"group-group_label_correlation-age-{sign}_measure-feature_label_fwhm-20"
    assert model.filename_root(contrast) == expected


def test_group_glm_instantiation(df):
    from clinica.pipelines.statistics_surface._model import GroupGLM

    with pytest.raises(
        ValueError,
        match="Contrast should refer to a categorical variable for group comparison.",
    ):
        GroupGLM("1 + age", df, "feature_label", "age", "group_label")
    model = GroupGLM("1 + age", df, "feature_label", "sex", "group_label")
    assert not model.with_interaction
    assert model.group_label == "group_label"
    assert model.fwhm == 20
    assert isinstance(model.contrasts, dict)
    contrast_names = ["Female-lt-Male", "Male-lt-Female"]
    assert set(model.contrasts.keys()) == set(contrast_names)
    for contrast_name, sign in zip(contrast_names, [-1, 1]):
        assert_array_equal(
            model.contrasts[contrast_name].values,
            sign * np.array([1, -1, 1, 1, 1, -1, 1]),
        )
    with pytest.raises(ValueError, match="Unknown contrast foo"):
        model.filename_root("foo")
    for contrast_name in contrast_names:
        expected = f"group-group_label_{contrast_name}_measure-feature_label_fwhm-20"
        assert model.filename_root(contrast_name) == expected


def test_group_glm_with_interaction_instantiation(df):
    from clinica.pipelines.statistics_surface._model import GroupGLMWithInteraction

    with pytest.raises(
        ValueError,
        match=(
            "The contrast must be an interaction between one continuous "
            "variable and one categorical variable."
        ),
    ):
        GroupGLMWithInteraction("1 + age", df, "feature_label", "age", "group_label")
    model = GroupGLMWithInteraction(
        "1 + age", df, "feature_label", "age * sex", "group_label"
    )
    assert model.with_interaction
    assert model.group_label == "group_label"
    assert model.fwhm == 20
    assert isinstance(model.contrasts, dict)
    assert len(model.contrasts) == 1
    assert_array_equal(
        model.contrasts["age * sex"].values,
        np.array([78.0, -73.4, 70.8, 82.3, 60.6, -72.1, 74.2]),
    )
    with pytest.raises(ValueError, match="Unknown contrast foo"):
        model.filename_root("foo")
    assert (
        model.filename_root("age * sex")
        == "interaction-age * sex_measure-feature_label_fwhm-20"
    )


def test_create_glm_model(df):
    from clinica.pipelines.statistics_surface._model import (
        CorrelationGLM,
        GroupGLM,
        GroupGLMWithInteraction,
        create_glm_model,
    )

    model = create_glm_model(
        "correlation", "age", df, "age", feature_label="feature_label"
    )
    assert isinstance(model, CorrelationGLM)
    model = create_glm_model(
        "group_comparison",
        "age",
        df,
        "sex",
        feature_label="feature_label",
        group_label="group_label",
    )
    assert isinstance(model, GroupGLM)
    model = create_glm_model(
        "group_comparison",
        "age * sex",
        df,
        "age * sex",
        feature_label="feature_label",
        group_label="group_label",
    )
    assert isinstance(model, GroupGLMWithInteraction)


def test_p_value_results():
    from clinica.pipelines.statistics_surface._model import PValueResults

    pvalues = np.random.random((5, 10))
    threshold = 0.3
    mask = pvalues >= threshold
    results = PValueResults(pvalues, mask, threshold)
    d = results.to_dict()
    assert_array_equal(d["P"], pvalues)
    assert_array_equal(d["mask"], mask)
    assert d["thresh"] == threshold
    d = results.to_json(indent=2)
    assert isinstance(d, str)
    for sub in ["P", "mask", "thresh"]:
        assert sub in d


def test_statistics_results_serializer(tmp_path):
    import json

    from scipy.io import loadmat

    from clinica.pipelines.statistics_surface._model import (
        CorrectedPValueResults,
        PValueResults,
        StatisticsResults,
        StatisticsResultsSerializer,
    )

    dummy_input = np.empty([3, 6])
    uncorrected = PValueResults(*[dummy_input] * 2, 0.1)
    corrected = CorrectedPValueResults(*[dummy_input] * 3, 0.2)
    results = StatisticsResults(*[dummy_input] * 2, uncorrected, dummy_input, corrected)
    serializer = StatisticsResultsSerializer(str(tmp_path / Path("out/dummy")))
    assert serializer.output_file == str(tmp_path / Path("out/dummy"))
    assert serializer.json_extension == "_results.json"
    assert serializer.json_indent == 4
    assert serializer.mat_extension == ".mat"
    with pytest.raises(
        NotImplementedError, match="Serializing method foo is not implemented."
    ):
        serializer.save(results, "foo")
    serializer.save(results, "json")
    assert os.path.exists(tmp_path / Path("out/dummy_results.json"))
    with open(tmp_path / Path("out/dummy_results.json"), "r") as fp:
        serialized = json.load(fp)
    serializer.save(results, "mat")
    names = [
        "coefficients",
        "TStatistics",
        "uncorrectedPValue",
        "correctedPValue",
        "FDR",
    ]
    keys = [
        "coef",
        "tvaluewithmask",
        "uncorrectedpvaluesstruct",
        "correctedpvaluesstruct",
        "FDR",
    ]
    for name, key in zip(names, keys):
        assert os.path.exists(tmp_path / Path(f"out/dummy_{name}.mat"))
        mat = loadmat(tmp_path / Path(f"out/dummy_{name}.mat"))
        if key in ["uncorrectedpvaluesstruct", "correctedpvaluesstruct"]:
            assert_array_almost_equal(mat[key]["P"][0, 0], dummy_input)
            assert_array_almost_equal(mat[key]["mask"][0, 0], dummy_input)
        else:
            assert_array_almost_equal(mat[key], dummy_input)
