"""Logging helpers for the :mod:`mcal` package.

Warnings and other diagnostics are emitted through the standard
:mod:`logging` module to ``stderr`` so that they stay separate from the
program results printed to ``stdout``. Colour is applied only when
``stderr`` is an interactive terminal.
"""

from __future__ import annotations

import logging
import sys

_RESET = "\033[0m"
_LEVEL_COLORS: dict[int, str] = {
    logging.WARNING: "\033[33m",  # yellow
    logging.ERROR: "\033[31m",  # red
    logging.CRITICAL: "\033[31m",  # red
}


class _ColorFormatter(logging.Formatter):
    """Formatter that prefixes the level name, colourised on a TTY.

    Parameters
    ----------
    color : bool
        If ``True``, wrap the formatted record in an ANSI colour code that
        depends on the record's level. If ``False``, emit plain text with no
        escape sequences (suitable for redirected or piped output).
    """

    def __init__(self, *, color: bool) -> None:
        super().__init__("%(levelname)s: %(message)s")
        self._color = color

    def format(self, record: logging.LogRecord) -> str:
        """Format a log record, optionally wrapping it in an ANSI colour code.

        Parameters
        ----------
        record : logging.LogRecord
            The record to format.

        Returns
        -------
        str
            The formatted message, colourised when ``color`` is enabled and a
            colour is defined for the record's level.
        """
        text = super().format(record)
        color = _LEVEL_COLORS.get(record.levelno)
        if self._color and color is not None:
            return f"{color}{text}{_RESET}"
        return text


def configure_logging(level: int = logging.WARNING) -> None:
    """Attach a ``stderr`` handler to the ``mcal`` package logger.

    Idempotent: calling it more than once leaves the existing handler in
    place and only updates the level. Intended to be called once from the
    CLI entry point. When never called (e.g. library use), logging falls
    back to its last-resort handler, which still emits ``WARNING`` and above
    to ``stderr``.

    Parameters
    ----------
    level : int, optional
        Minimum level the package logger emits, by default
        :data:`logging.WARNING`.

    Returns
    -------
    None
    """
    logger = logging.getLogger("mcal")
    logger.setLevel(level)
    if logger.handlers:
        return
    handler = logging.StreamHandler(sys.stderr)
    handler.setFormatter(_ColorFormatter(color=sys.stderr.isatty()))
    logger.addHandler(handler)
    logger.propagate = False


def get_logger(name: str) -> logging.Logger:
    """Return a logger for the given name.

    Parameters
    ----------
    name : str
        Logger name, conventionally the module ``__name__`` so the logger
        becomes a child of the ``mcal`` package logger.

    Returns
    -------
    logging.Logger
        The requested logger.
    """
    return logging.getLogger(name)
