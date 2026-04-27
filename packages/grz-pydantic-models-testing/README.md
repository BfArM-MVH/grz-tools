# grz-pydantic-models-testing

Installable test fixtures and example metadata for [`grz-pydantic-models`](../grz-pydantic-models).

This package ships the example and failing metadata JSON files as proper package data so
that any downstream package can access them via `importlib.resources` without relying on
filesystem layout or editable-install assumptions.

## Installation

```bash
pip install grz-pydantic-models-testing
```

## Usage

```python
import importlib.resources
import grz_pydantic_models_testing
from grz_pydantic_models_testing import example_metadata, failing_metadata

# Read an example metadata file
json_text = (
    importlib.resources.files(example_metadata)
    .joinpath("wes_tumor_germline", "v1.2.1.json")
    .read_text()
)

# Parse it directly
from grz_pydantic_models.submission.metadata import GrzSubmissionMetadata
metadata = GrzSubmissionMetadata.model_validate_json(json_text)
```

## Contents

| Sub-package | Description |
|---|---|
| `example_metadata` | Valid example submissions for all supported study types and schema versions |
| `example_pruefberichte` | Example Prüfbericht payloads |
| `example_research_consent` | Example research consent FHIR resources |
| `failing_metadata` | Invalid metadata fixtures for negative testing |
