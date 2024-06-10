from pathlib import Path

from clinica.utils.filemanip import UserProvidedPath

__all__ = ["validate_input_path"]


def validate_input_path(input_path: UserProvidedPath, check_exist: bool = True) -> Path:
    """Take a user provided path as input and returns a resolved, existing if check_exist is True, Path."""
    input_path = Path(input_path).resolve()
    if check_exist and not input_path.exists():
        raise FileNotFoundError(f"The provided input path {input_path} does not exist.")
    return input_path
