"""
Configuration loader for the revised TCO2 study.
Loads YAML configs and provides convenient access.
"""

import yaml
import os
from pathlib import Path


def find_project_root():
    """Find the project root directory (contains config/)."""
    current = Path(__file__).resolve().parent.parent
    if (current / "config").exists():
        return current
    raise FileNotFoundError("Cannot find project root (config/ directory)")


def load_config(config_name="global_config.yaml"):
    """Load a YAML configuration file from the config directory."""
    root = find_project_root()
    config_path = root / "config" / config_name
    with open(config_path, "r") as f:
        return yaml.safe_load(f)


def get_path(key, config=None):
    """Resolve a path from the global config, relative to project root."""
    if config is None:
        config = load_config()
    root = find_project_root()
    return root / config["paths"][key]


def get_data_path(config=None):
    """Get the raw GLODAP data file path."""
    return get_path("raw_data", config)


def get_basin_config(config=None):
    """Get basin definitions from config."""
    if config is None:
        config = load_config()
    return config["basins"]


def get_qc_config(config=None):
    """Get quality control thresholds from config."""
    if config is None:
        config = load_config()
    return config["quality_control"]
