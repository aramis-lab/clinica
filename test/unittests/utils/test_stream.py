from typing import Type

import pytest

from clinica.utils.exceptions import ClinicaMissingDependencyError
from clinica.utils.stream import LoggingLevel


@pytest.mark.parametrize(
    "level,expected",
    [
        ("debug", LoggingLevel.DEBUG),
        ("info", LoggingLevel.INFO),
        ("warning", LoggingLevel.WARNING),
        ("error", LoggingLevel.ERROR),
        ("critical", LoggingLevel.CRITICAL),
        (1, LoggingLevel.DEBUG),
        (2, LoggingLevel.INFO),
        (3, LoggingLevel.WARNING),
        (4, LoggingLevel.ERROR),
        (5, LoggingLevel.CRITICAL),
    ],
)
def test_get_logging_level(level, expected):
    from clinica.utils.stream import get_logging_level

    assert get_logging_level(level) == expected


def test_get_logging_level_error():
    from clinica.utils.stream import get_logging_level

    with pytest.raises(
        ValueError,
        match="Logging level foo is not valid.",
    ):
        get_logging_level("foo")


@pytest.mark.parametrize(
    "message,error_type",
    [
        ("foo bar error", FileNotFoundError),
        ("bar baz error", RuntimeError),
        ("missing dependency", ClinicaMissingDependencyError),
    ],
)
def test_log_and_raise(message: str, error_type: Type[Exception]):
    from clinica.utils.stream import log_and_raise

    with pytest.raises(
        error_type,
        match=message,
    ):
        log_and_raise(message, error_type)


@pytest.mark.parametrize(
    "message,warning_type",
    [
        ("foo bar warning", UserWarning),
        ("bar baz warning", DeprecationWarning),
        ("warning", FutureWarning),
    ],
)
def test_log_and_warn(message: str, warning_type: Type[Warning]):
    from clinica.utils.stream import log_and_warn

    with pytest.warns(
        warning_type,
        match=message,
    ):
        log_and_warn(message, warning_type)
