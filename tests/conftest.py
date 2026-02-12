import pytest


def pytest_addoption(parser):
    """Add custom CLI options for redundancy test selection."""
    parser.addoption(
        "--SG-redundancy",
        action="store",
        default="0",
        help="Run SG redundancy test with N cases (e.g. 1, 5, all). 0 disables.",
    )
    parser.addoption(
        "--MSG-redundancy",
        action="store",
        default="0",
        help="Run MSG redundancy test with N cases (e.g. 1, 5, all). 0 disables.",
    )
