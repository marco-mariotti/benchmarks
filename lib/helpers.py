import json
from pathlib import Path


JSON = str | int | float | bool | None | dict[str, "JSON"] | list["JSON"]


def get_config() -> JSON:
    return json.load(Path("config.json").open())


def library_to_language(library: str) -> str:
    return get_config()["library_to_language"].get(library)
