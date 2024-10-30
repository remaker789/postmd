from collections import UserDict
import inspect
import logging

import numpy as np

logger = logging.getLogger(__name__)

class Results(UserDict):
    r"""Container object for storing results.

    :class:`Results` are dictionaries that provide two ways by which values
    can be accessed: by dictionary key ``results["value_key"]`` or by object
    attribute, ``results.value_key``. :class:`Results` stores all results
    obtained from an analysis after calling :meth:`~AnalysisBase.run()`.

    The implementation is copied from MDAnalysis, similar to the :class:`sklearn.utils.Bunch`
    class in `scikit-learn`_.

    .. _`scikit-learn`: https://scikit-learn.org/

    Raises
    ------
    AttributeError
        If an assigned attribute has the same name as a default attribute.

    ValueError
        If a key is not of type ``str`` and therefore is not able to be
        accessed by attribute.

    Examples
    --------
    >>> from postmd.analysis.base import Results
    >>> results = Results(a=1, b=2)
    >>> results['b']
    2
    >>> results.b
    2
    >>> results.a = 3
    >>> results['a']
    3
    >>> results.c = [1, 2, 3, 4]
    >>> results['c']
    [1, 2, 3, 4]


    .. versionadded:: 2.0.0
    """

    def _validate_key(self, key):
        if key in dir(self):
            raise AttributeError(f"'{key}' is a protected dictionary "
                                 "attribute")
        elif isinstance(key, str) and not key.isidentifier():
            raise ValueError(f"'{key}' is not a valid attribute")

    def __init__(self, *args, **kwargs):
        kwargs = dict(*args, **kwargs)
        if "data" in kwargs.keys():
            raise AttributeError(f"'data' is a protected dictionary attribute")
        self.__dict__["data"] = {}
        self.update(kwargs)

    def __setitem__(self, key, item):
        self._validate_key(key)
        super().__setitem__(key, item)

    def __setattr__(self, attr, val):
        if attr == 'data':
            super().__setattr__(attr, val)
        else:
            self.__setitem__(attr, val)

    def __getattr__(self, attr):
        try:
            return self[attr]
        except KeyError as err:
            raise AttributeError("'Results' object has no "
                                 f"attribute '{attr}'") from err

    def __delattr__(self, attr):
        try:
            del self[attr]
        except KeyError as err:
            raise AttributeError("'Results' object has no "
                                 f"attribute '{attr}'") from err

    def __getstate__(self):
        return self.data

    def __setstate__(self, state):
        self.data = state