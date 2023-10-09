import numpy as np
import pandas as pd

from ._utils import check_column_in_df, is_categorical

__all__ = [
    "Contrast",
    "CorrelationContrast",
    "GroupContrast",
    "GroupContrastWithInteraction",
]


class Contrast:
    def __init__(self, name: str, built_contrast: pd.Series):
        self.name = name
        self.built_contrast = built_contrast


class CorrelationContrast(Contrast):
    def __init__(
        self, name: str, absolute_name: str, built_contrast: pd.Series, sign: str
    ):
        super().__init__(name, built_contrast)
        self.sign = sign
        self.absolute_name = absolute_name

    @classmethod
    def from_string(cls, contrast: str, df: pd.DataFrame):
        absolute_contrast_name = contrast
        contrast_sign = "positive"
        if contrast.startswith("-"):
            absolute_contrast_name = contrast[1:].lstrip()
            contrast_sign = "negative"
        check_column_in_df(df, absolute_contrast_name)
        built_contrast = df[absolute_contrast_name]
        if contrast_sign == "negative":
            built_contrast *= -1
        return cls(contrast, absolute_contrast_name, built_contrast, contrast_sign)


class GroupContrastWithInteraction(Contrast):
    def __init__(self, name: str, built_contrast: pd.Series):
        super().__init__(name, built_contrast)

    @classmethod
    def from_string(cls, contrast: str, df: pd.DataFrame):
        contrast_elements = [_.strip() for _ in contrast.split("*")]
        for contrast_element in contrast_elements:
            check_column_in_df(df, contrast_element)
        categorical = [is_categorical(df, _) for _ in contrast_elements]
        if len(contrast_elements) != 2 or sum(categorical) != 1:
            raise ValueError(
                "The contrast must be an interaction between one continuous "
                "variable and one categorical variable. Your contrast contains "
                f"the following variables : {contrast_elements}"
            )
        idx = 0 if categorical[0] else 1
        categorical_contrast = contrast_elements[idx]
        continue_contrast = contrast_elements[(idx + 1) % 2]
        group_values = np.unique(df[categorical_contrast])
        built_contrast = df[continue_contrast].where(
            df[categorical_contrast] == group_values[0], 0
        ) - df[continue_contrast].where(df[categorical_contrast] == group_values[1], 0)
        return cls(contrast, built_contrast)


class GroupContrast(Contrast):
    def __int__(self, name: str, built_contrast: pd.Series):
        super().__init__(name, built_contrast)

    @classmethod
    def from_string(cls, contrast: str, df: pd.DataFrame, sign: str):
        check_column_in_df(df, contrast)
        if not is_categorical(df, contrast):
            raise ValueError(
                "Contrast should refer to a categorical variable for group comparison. "
                "Please select 'correlation' for 'glm_type' otherwise."
            )
        group_values = np.unique(df[contrast])
        i, j = (0, 1) if sign == "positive" else (1, 0)
        contrast_name = f"{group_values[j]}-lt-{group_values[i]}"
        built_contrast = (df[contrast] == group_values[i]).astype(int) - (
            df[contrast] == group_values[j]
        ).astype(int)

        return cls(contrast_name, built_contrast)
