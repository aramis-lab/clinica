from pathlib import Path
from typing import Union

__all__ = ["validate_input_path"]


def validate_input_path(input_path: Union[str, Path], check_exist: bool = True) -> Path:
    input_path = Path(input_path).resolve()
    if check_exist and not input_path.exists():
        raise FileNotFoundError(f"The provided input path {input_path} does not exist.")
    return input_path
