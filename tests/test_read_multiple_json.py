"""Tests for `read_multiple_json` from the file_operations module."""

import io

import pytest

from grz_cli.utils.io import read_multiple_json


# Basic JSON objects (single-line JSON objects)
def test_read_multiple_json_basic():
    json_data = """{"name": "Alice", "age": 30}
                   {"name": "Bob", "age": 25}
                   {"name": "Charlie", "age": 35}"""
    input_stream = io.StringIO(json_data)

    result = list(read_multiple_json(input_stream))

    assert result == [
        {"name": "Alice", "age": 30},
        {"name": "Bob", "age": 25},
        {"name": "Charlie", "age": 35},
    ]


# Basic JSON objects without white spaces
def test_read_multiple_json_nowhitespace():
    json_data = """{"name": "Alice", "age": 30}{"name": "Bob", "age": 25}{"name": "Charlie", "age": 35}"""
    input_stream = io.StringIO(json_data)

    result = list(read_multiple_json(input_stream))

    assert result == [
        {"name": "Alice", "age": 30},
        {"name": "Bob", "age": 25},
        {"name": "Charlie", "age": 35},
    ]


# Multiline JSON objects (JSON objects span multiple lines)
def test_read_multiple_json_multiline():
    json_data = """{
                    "name": "Alice",
                    "age": 30
                   }
                   {
                    "name": "Bob",
                    "age": 25
                   }"""
    input_stream = io.StringIO(json_data)

    result = list(read_multiple_json(input_stream))

    assert result == [{"name": "Alice", "age": 30}, {"name": "Bob", "age": 25}]


# JSON objects that span across chunks (larger buffer)
def test_read_multiple_json_large_buffer():
    json_data = """{"name": "Alice", "age": 30}
                   {"name": "Bob", "age": 25}
                   {"name": "Charlie", "age": 35}"""
    input_stream = io.StringIO(json_data)

    # Simulate reading in smaller buffer chunks to ensure that multiple reads work properly
    result = list(read_multiple_json(input_stream, buffer_size=10))

    assert result == [
        {"name": "Alice", "age": 30},
        {"name": "Bob", "age": 25},
        {"name": "Charlie", "age": 35},
    ]


# Incomplete JSON objects (incomplete at the end)
def test_read_multiple_json_incomplete():
    json_data = """{"name": "Alice", "age": 30}
                   {"name": "Bob", "age": 25}
                   {"name": "Charlie", "age": 35
                   """
    input_stream = io.StringIO(json_data)

    with pytest.raises(ValueError):
        # The last object is incomplete
        _result = list(read_multiple_json(input_stream))


# Empty input
def test_read_multiple_json_empty():
    json_data = ""
    input_stream = io.StringIO(json_data)

    result = list(read_multiple_json(input_stream))

    assert result == []


# Invalid JSON (no valid objects)
def test_read_multiple_json_invalid():
    json_data = """invalid data"""
    input_stream = io.StringIO(json_data)

    with pytest.raises(ValueError):
        _result = list(read_multiple_json(input_stream))


# Mixed valid and invalid JSON objects
def test_read_multiple_json_mixed_valid_invalid():
    json_data = """{"name": "Alice", "age": 30}
                   invalid json here
                   {"name": "Bob", "age": 25}"""
    input_stream = io.StringIO(json_data)

    with pytest.raises(ValueError):
        _result = list(read_multiple_json(input_stream))
