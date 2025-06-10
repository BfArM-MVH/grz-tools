import importlib.resources

from grz_pydantic_models.pruefbericht import Pruefbericht

from . import resources


def test_example():
    metadata_str = (
        importlib.resources.files(resources).joinpath("example_pruefberichte", "submission_example.json").read_text()
    )
    Pruefbericht.model_validate_json(metadata_str)
