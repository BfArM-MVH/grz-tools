"""
grz-pydantic-models-testing
============================
Installable test fixtures and example data for the GRZ pydantic models.

Usage with ``importlib.resources``::

    import importlib.resources
    import grz_pydantic_models_testing

    # Access a specific example file
    json_text = (
        importlib.resources.files(grz_pydantic_models_testing)
        .joinpath("example_metadata", "wes_tumor_germline", "v1.2.1.json")
        .read_text()
    )

    # Or via the sub-package directly
    from grz_pydantic_models_testing import example_metadata

    json_text = (
        importlib.resources.files(example_metadata)
        .joinpath("wes_tumor_germline", "v1.2.1.json")
        .read_text()
    )
"""

__version__ = "0.1.0"
