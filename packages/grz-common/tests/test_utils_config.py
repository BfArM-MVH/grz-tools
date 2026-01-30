import copy

import pytest
from grz_common.utils.config import merge_config_dicts


def test_merge_config_dicts_no_conflict():
    a = {
        "key1": "value1",
        "key2": {
            "subkey1": [1, 2],
            "subkey2": "subvalue2",
        },
        "key3": [1, 2, 3],
    }
    b = {
        "key2": {
            "subkey1": [3, 4],
            "subkey3": "subvalue3",
        },
        "key4": "value4",
    }
    expected = {
        "key1": "value1",
        "key2": {
            "subkey1": [3, 4],  # from b, replaces a's list
            "subkey2": "subvalue2",
            "subkey3": "subvalue3",
        },
        "key3": [1, 2, 3],
        "key4": "value4",
    }
    result = merge_config_dicts(a, b)
    assert result == expected


def test_merge_config_dicts_missing_values():
    a = {
        "key2": {
            "subkey1": None,
        },
        "key3": [1, 2, 3],
    }
    b = {
        "key2": {
            "subkey1": [1, 2],
        },
        "key3": None,
    }
    expected = {
        "key2": {
            "subkey1": [1, 2],  # from b, as a's value is None
        },
        "key3": [1, 2, 3],  # remains unchanged from a as b's value is None
    }
    result = merge_config_dicts(a, b)
    assert result == expected


def test_merge_config_dicts_conflict():
    a = {
        "key1": {
            "subkey1": "value1",
        }
    }
    b = {
        "key1": {
            "subkey1": [1, 2, 3],
        }
    }
    with pytest.raises(ValueError) as excinfo:
        merge_config_dicts(a, b)
    assert "Conflict at key1.subkey1" in str(excinfo.value)


def test_merge_config_dicts_unmodified_inputs():
    a = {
        "key1": "value1",
        "key2": {
            "subkey1": [1, 2],
        },
    }
    b = {
        "key2": {
            "subkey1": [3, 4],
            "subkey2": "subvalue2",
        },
        "key3": "value3",
    }

    # Create deep copies for later comparison
    a_original = copy.deepcopy(a)
    b_original = copy.deepcopy(b)

    result = merge_config_dicts(a, b)

    # Modify the result
    result["key1"] = "new_value"
    result["key2"]["subkey1"].append(5)
    result["key3"] = "new_value3"

    # Check that original dicts are unmodified
    assert a == a_original
    assert b == b_original
