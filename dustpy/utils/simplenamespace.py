import types


class SimpleNamespace(types.SimpleNamespace):
    """Class for SimpleNamespace that does not allow to add new attributes."""

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)

    def __setattr__(self, key, val):
        if key not in self.__dict__:
            raise RuntimeError(f"'{key}' is not a valid option.")
        super().__setattr__(key, val)
