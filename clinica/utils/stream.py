"""This module handles stream and log redirection."""
import warnings
from enum import Enum
from typing import Type


class LoggingLevel(str, Enum):
    debug = "debug"
    info = "info"
    warning = "warning"
    error = "error"
    critical = "critical"


def cprint(msg: str, lvl: str = "info") -> None:
    """
    Print message to the console at the desired logging level.

    Args:
        msg (str): Message to print.
        lvl (str): Logging level between "debug", "info", "warning", "error" and "critical".
                   The default value is "info".
    """
    from logging import getLogger

    # Use the package level logger.
    logger = getLogger("clinica")

    # Log message as info level.
    if lvl == LoggingLevel.debug:
        logger.debug(msg=msg)
    elif lvl == LoggingLevel.info:
        logger.info(msg=msg)
    elif lvl == LoggingLevel.warning:
        logger.warning(msg=msg)
    elif lvl == LoggingLevel.error:
        logger.error(msg=msg)
    elif lvl == LoggingLevel.critical:
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
    cprint(message, lvl="error")
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
    cprint(message, lvl="warning")
    warnings.warn(message, warning_type)
