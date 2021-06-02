# coding: utf8

"""This module handles stream and log redirection."""


def cprint(msg):
    from logging import getLogger

    # Use the package level logger.
    logger = getLogger("clinica")

    # Log message as info level.
    logger.info(msg=msg)
