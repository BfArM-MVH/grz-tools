import importlib.resources

from grz_pydantic_models.pruefbericht import Pruefbericht
from grz_pydantic_models_testing import example_pruefberichte


def test_example():
    metadata_str = importlib.resources.files(example_pruefberichte).joinpath("submission_example.json").read_text()
    Pruefbericht.model_validate_json(metadata_str)