"""This module handles stream and log redirection."""
import warnings
from enum import IntEnum
from typing import Optional, Type, Union

__all__ = [
    "LoggingLevel",
    "get_logging_level",
    "cprint",
    "log_and_raise",
    "log_and_warn",
]


class LoggingLevel(IntEnum):
    DEBUG = 1
    INFO = 2
    WARNING = 3
    ERROR = 4
    CRITICAL = 5


def get_logging_level(level: Union[str, int, LoggingLevel]) -> LoggingLevel:
    if isinstance(level, str):
        if level == "debug":
            return LoggingLevel.DEBUG
        if level == "info":
            return LoggingLevel.INFO
        if level == "warning":
            return LoggingLevel.WARNING
        if level == "error":
            return LoggingLevel.ERROR
        if level == "critical":
            return LoggingLevel.CRITICAL
        raise ValueError(f"Logging level {level} is not valid.")
    return LoggingLevel(level)


def cprint(msg: str, lvl: Optional[Union[str, int, LoggingLevel]] = None) -> None:
    """Print message to the console at the desired logging level.

    Parameters
    ----------
    msg : str
        The message to log.

    lvl : str or int, or LoggingLevel, optional
        The logging level:
            - 1 = "debug"
            - 2 = "info"
            - 3 = "warning"
            - 4 = "error"
            - 5 = "critical"
        The default value is "info".
    """
    from logging import getLogger

    lvl = lvl or LoggingLevel.INFO
    # Use the package level logger.
    logger = getLogger("clinica")

    # Log message as info level.
    if lvl == LoggingLevel.DEBUG:
        logger.debug(msg=msg)
    elif lvl == LoggingLevel.INFO:
        logger.info(msg=msg)
    elif lvl == LoggingLevel.WARNING:
        logger.warning(msg=msg)
    elif lvl == LoggingLevel.ERROR:
        logger.error(msg=msg)
    elif lvl == LoggingLevel.CRITICAL:
        logger.critical(msg=msg)
    else:
        pass


def log_and_raise(message: str, error_type: Type[Exception]):
    """Log the error message using cprint and raise.

    Parameters
    ----------
    message : str
        The error message.

    error_type : Exception
        The error type to raise.
    """
    cprint(message, lvl=LoggingLevel.ERROR)
    raise error_type(message)


def log_and_warn(message: str, warning_type: Type[Warning]):
    """Log the warning message using cprint and give it to the user.

    Parameters
    ----------
    message : str
        The warning message.

    warning_type : Warning
        The warning type to use.
    """
    cprint(message, lvl=LoggingLevel.WARNING)
    warnings.warn(message, warning_type)
