# Copyright (c) SHRY Development Team.
# Distributed under the terms of the MIT License.

"""Some helper functions related to testing."""

import os
from functools import wraps


def chdir(directory):
    """
    Wrapper for running a function inside `directory', and exit immediately afterwards.

    Args:
        directory: Target directory."""

    def wrapper(function):
        @wraps(function)
        def wrapped(*args, **kwargs):
            cwd = os.getcwd()
            try:
                os.chdir(testing_dir(directory))
                return function(*args, **kwargs)
            finally:
                os.chdir(cwd)

        return wrapped

    return wrapper


def testing_dir(path="."):
    """Get relative directory from the current testing directory."""
    return os.path.join(os.path.dirname(os.path.realpath(__file__)), path)
