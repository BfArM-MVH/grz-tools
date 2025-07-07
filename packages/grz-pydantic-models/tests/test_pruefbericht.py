import importlib.resources

import pytest
from grz_pydantic_models.pruefbericht import Pruefbericht

from . import resources


@pytest.mark.parametrize("example_name", ("submission_example", "test_submission_example"))
def test_example(example_name: str):
    metadata_str = (
        importlib.resources.files(resources).joinpath("example_pruefberichte", f"{example_name}.json").read_text()
    )
    Pruefbericht.model_validate_json(metadata_str)
