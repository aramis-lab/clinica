"""This module handles stream and log redirection."""

from enum import Enum


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
