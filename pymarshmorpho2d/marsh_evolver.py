#! /usr/bin/env python
"""pyMarshMorpho2D model."""

import numpy as np
from landlab import Component


class MarshEvolver(Component):
    """Simulate tidal marsh evolution."""

    _name = "MarshEvolver"

    _cite_as = """
    """

    _info = {
        "topographic__elevation": {
            "dtype": float,
            "intent": "inout",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Land surface topographic elevation",
        },
    }

    def __init__(self, grid):
        """Initialize the MarshEvolver."""
        super(MarshEvolver, self).__init__(grid)
