import os
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from numpy.testing import assert_array_almost_equal, assert_array_equal

CURRENT_DIR = os.path.dirname(os.path.realpath(__file__))


def is_brainstat_missing() -> bool:
    try:
        import brainstat
    except ImportError:
        return True
    return False


@pytest.fixture
def df():
    return pd.read_csv(Path(CURRENT_DIR) / "data/subjects.tsv", sep="\t")


@pytest.mark.skipif(is_brainstat_missing(), reason="Brainstat is not installed.")
def test_missing_column_error(df):
    from clinica.pipelines.statistics_surface.surfstat.models._utils import (
        check_column_in_df,
    )

    with pytest.raises(
        ValueError,
        match=(
            "Term foo from the design matrix is not in the columns of the "
            "provided TSV file. Please make sure that there is no typo"
        ),
    ):
        check_column_in_df(df, "foo")


@pytest.mark.skipif(is_brainstat_missing(), reason="Brainstat is not installed.")
def test_is_categorical(df):
    from clinica.pipelines.statistics_surface.surfstat.models._utils import (
        is_categorical,
    )

    assert is_categorical(df, "sex")
    assert not is_categorical(df, "age")


@pytest.mark.skipif(is_brainstat_missing(), reason="Brainstat is not installed.")
def test_build_model_term_error(df):
    from brainstat.stats.terms import FixedEffect

    from clinica.pipelines.statistics_surface.surfstat.models._utils import (
        _build_model_term,
    )

    assert isinstance(_build_model_term("sex", df), FixedEffect)


@pytest.mark.skipif(is_brainstat_missing(), reason="Brainstat is not installed.")
@pytest.mark.parametrize("design", ["1 + age", "1+age", "age +1", "age"])
def test_build_model_intercept(design, df):
    """Test that we get the same results with equivalent designs.
    Especially, the fact that adding explicitly the intercept doesn't change the results.
    Test also that spaces in the design expression have no effect.
    """
    from brainstat.stats.terms import FixedEffect

    from clinica.pipelines.statistics_surface.surfstat.models._utils import build_model

    model = build_model(design, df)
    assert isinstance(model, FixedEffect)
    assert len(model.m.columns) == 2
    assert_array_equal(model.intercept, np.array([1, 1, 1, 1, 1, 1, 1]))
    assert_array_equal(model.age, np.array([78.0, 73.4, 70.8, 82.3, 60.6, 72.1, 74.2]))


@pytest.mark.skipif(is_brainstat_missing(), reason="Brainstat is not installed.")
def test_build_model(df):
    from brainstat.stats.terms import FixedEffect

    from clinica.pipelines.statistics_surface.surfstat.models._utils import build_model

    model = build_model("1 + age + sex", df)
    assert isinstance(model, FixedEffect)
    assert len(model.m.columns) == 4
    assert_array_equal(model.intercept, np.array([1, 1, 1, 1, 1, 1, 1]))
    assert_array_equal(model.age, np.array([78.0, 73.4, 70.8, 82.3, 60.6, 72.1, 74.2]))
    assert_array_equal(model.sex_Female, np.array([1, 0, 1, 1, 1, 0, 1]))
    assert_array_equal(model.sex_Male, np.array([0, 1, 0, 0, 0, 1, 0]))
    model = build_model("1 + age + sex + age * sex", df)
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


@pytest.mark.skipif(is_brainstat_missing(), reason="Brainstat is not installed.")
def test_base_glm_instantiation_error(df):
    """Test that the base abstract GLM class cannot be instantiated."""
    from clinica.pipelines.statistics_surface.surfstat.models._base import GLM

    with pytest.raises(NotImplementedError):
        GLM("1 + age", df, "feature_label", "age")


@pytest.mark.skipif(is_brainstat_missing(), reason="Brainstat is not installed.")
@pytest.mark.parametrize(
    "model,contrast",
    [
        ("CorrelationGLM", "age"),
        ("GroupGLM", "sex"),
        ("GroupGLMWithInteraction", "age * sex"),
    ],
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
def test_common_parameters_glm_instantiation(df, model, contrast, parameters):
    from brainstat.stats.terms import FixedEffect

    from clinica.pipelines.statistics_surface.surfstat.models._correlation import (
        CorrelationGLM,
    )
    from clinica.pipelines.statistics_surface.surfstat.models._group import (
        GroupGLM,
        GroupGLMWithInteraction,
    )

    if model == "CorrelationGLM":
        model = CorrelationGLM
    elif model == "GroupGLM":
        model = GroupGLM
    elif model == "GroupGLMWithInteraction":
        model = GroupGLMWithInteraction
    else:
        raise ValueError(f"Unknown model {model}.")

    model_instance = model("1 + age", df, "feature_label", contrast, "group")

    assert not model_instance._two_tailed
    assert model_instance._correction == ["fdr", "rft"]
    assert model_instance.feature_label == "feature_label"
    assert model_instance.fwhm == 20
    assert model_instance.threshold_uncorrected_pvalue == 0.001
    assert model_instance.threshold_corrected_pvalue == 0.05
    assert model_instance.cluster_threshold == 0.001
    assert model_instance.filenames == dict()
    assert isinstance(model_instance.model, FixedEffect)


@pytest.mark.skipif(is_brainstat_missing(), reason="Brainstat is not installed.")
@pytest.mark.parametrize("contrast", ["age", "-age"])
def test_correlation_glm_instantiation(df, contrast):
    from clinica.pipelines.statistics_surface.surfstat.models._contrast import (
        CorrelationContrast,
    )
    from clinica.pipelines.statistics_surface.surfstat.models._correlation import (
        CorrelationGLM,
    )

    model = CorrelationGLM("1 + age", df, "feature_label", contrast, "group")
    assert not model.with_interaction
    assert model.group_label == "group"
    assert model.feature_label == "feature_label"
    assert model.fwhm == 20
    sign = "positive" if contrast == "age" else "negative"
    assert len(model.contrasts) == 1
    assert isinstance(model.contrasts[0], CorrelationContrast)
    assert model.contrasts[0].sign == sign
    assert model.contrasts[0].name == contrast
    mult = 1 if sign == "positive" else -1
    assert_array_equal(
        model.contrasts[0].built_contrast.values,
        mult * np.array([78.0, 73.4, 70.8, 82.3, 60.6, 72.1, 74.2]),
    )
    with pytest.raises(ValueError, match="Unknown contrast foo"):
        model.get_output_filename("foo")
    expected = f"group-group_correlation-age-{sign}_measure-feature_label_fwhm-20"
    assert model.get_output_filename(contrast) == expected


@pytest.mark.skipif(is_brainstat_missing(), reason="Brainstat is not installed.")
def test_group_glm_instantiation(df):
    from clinica.pipelines.statistics_surface.surfstat.models._group import GroupGLM

    with pytest.raises(
        ValueError,
        match="Contrast should refer to a categorical variable for group comparison.",
    ):
        GroupGLM("1 + age", df, "feature_label", "age", "group_label")
    model = GroupGLM("1 + age", df, "feature_label", "sex", "group_label")
    assert not model.with_interaction
    assert model.group_label == "group_label"
    assert model.fwhm == 20
    assert isinstance(model.contrasts, list)
    contrast_names = ["Female-lt-Male", "Male-lt-Female"]
    assert set(model.contrast_names) == set(contrast_names)
    for contrast_name, sign in zip(contrast_names, [-1, 1]):
        assert_array_equal(
            model.get_contrast_by_name(contrast_name).built_contrast,
            sign * np.array([1, -1, 1, 1, 1, -1, 1]),
        )
    with pytest.raises(ValueError, match="Unknown contrast foo"):
        model.get_output_filename("foo")
    for contrast_name in contrast_names:
        expected = f"group-group_label_{contrast_name}_measure-feature_label_fwhm-20"
        assert model.get_output_filename(contrast_name) == expected


@pytest.mark.skipif(is_brainstat_missing(), reason="Brainstat is not installed.")
def test_group_glm_with_interaction_instantiation(df):
    from clinica.pipelines.statistics_surface.surfstat.models._group import (
        GroupGLMWithInteraction,
    )

    with pytest.raises(
        ValueError,
        match=(
            "The contrast must be an interaction between one continuous "
            "variable and one categorical variable."
        ),
    ):
        GroupGLMWithInteraction("1 + age", df, "feature_label", "age", "group")
    model = GroupGLMWithInteraction(
        "1 + age", df, "feature_label", "age * sex", "group"
    )
    assert model.with_interaction
    assert model.group_label == "group"
    assert model.fwhm == 20
    assert isinstance(model.contrasts, list)
    assert len(model.contrasts) == 1
    assert_array_equal(
        model.get_contrast_by_name("age * sex").built_contrast,
        np.array([78.0, -73.4, 70.8, 82.3, 60.6, -72.1, 74.2]),
    )
    with pytest.raises(ValueError, match="Unknown contrast foo"):
        model.get_output_filename("foo")
    assert (
        model.get_output_filename("age * sex")
        == "interaction-age * sex_measure-feature_label_fwhm-20"
    )


@pytest.mark.skipif(is_brainstat_missing(), reason="Brainstat is not installed.")
def test_create_glm_model(df):
    from clinica.pipelines.statistics_surface.surfstat.models import create_glm_model
    from clinica.pipelines.statistics_surface.surfstat.models._group import (
        GroupGLM,
        GroupGLMWithInteraction,
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


@pytest.mark.skipif(is_brainstat_missing(), reason="Brainstat is not installed.")
def test_p_value_results():
    from clinica.pipelines.statistics_surface.surfstat.models.results._statistics import (
        PValueResults,
    )

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


@pytest.mark.skipif(is_brainstat_missing(), reason="Brainstat is not installed.")
def test_statistics_results_serializer(tmp_path):
    import json

    from scipy.io import loadmat

    from clinica.pipelines.statistics_surface.surfstat.models.results import (
        StatisticsResults,
        StatisticsResultsSerializer,
    )
    from clinica.pipelines.statistics_surface.surfstat.models.results._statistics import (
        CorrectedPValueResults,
        PValueResults,
    )

    dummy_input = np.empty([3, 6])
    uncorrected = PValueResults(*[dummy_input] * 2, 0.1)
    corrected = CorrectedPValueResults(*[dummy_input] * 3, 0.2)
    results = StatisticsResults(*[dummy_input] * 2, uncorrected, dummy_input, corrected)

    out_file = tmp_path / "out" / "dummy"
    serializer = StatisticsResultsSerializer(out_file)
    assert serializer.output_file == out_file
    assert serializer.json_extension == "_results.json"
    assert serializer.json_indent == 4
    assert serializer.mat_extension == ".mat"

    with pytest.raises(
        NotImplementedError, match="Serializing method foo is not implemented."
    ):
        serializer.save(results, "foo")
    serializer.save(results, "json")
    assert (tmp_path / "out" / "dummy_results.json").exists()
    with open(tmp_path / "out" / "dummy_results.json", "r") as fp:
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
        assert (tmp_path / "out" / f"dummy_{name}.mat").exists()
        mat = loadmat(tmp_path / "out" / f"dummy_{name}.mat")
        if key in ("uncorrectedpvaluesstruct", "correctedpvaluesstruct"):
            assert_array_almost_equal(mat[key]["P"][0, 0], dummy_input)
            assert_array_almost_equal(mat[key]["mask"][0, 0], dummy_input)
        else:
            assert_array_almost_equal(mat[key], dummy_input)
