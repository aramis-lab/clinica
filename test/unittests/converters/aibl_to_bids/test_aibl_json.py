import pytest
from pydicom.dataset import Dataset

test = Dataset()
test.add(0x00100020, "LO", "12345")
