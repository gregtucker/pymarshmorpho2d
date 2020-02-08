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
        "water__depth": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Water depth",
        },
        "fully_wet__depth": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "m",
            "mapping": "node",
            "doc": "Water depth, with depth > 0 everywhere",
        },
    }

    def __init__(self, grid,
                 rel_sl_rise_rate=2.74e-6,
                 tidal_range=3.1,
                 ):
        """Initialize the MarshEvolver.

        Parameters
        ----------
        grid : ModelGrid object
            Landlab model grid
        rel_sl_rise_rate : float
            Rate of relative sea-level rise, m/day
        tidal_range : float
            Tidal range, m
        """

        super(MarshEvolver, self).__init__(grid)

        self._elev = self._grid.at_node['topographic__elevation']
        self._water_depth = self._grid.at_node['water__depth']
        self._fully_wet_depth = self._grid.at_node['fully_wet__depth']

        self._mean_sea_level = 0.0
        self._rel_sl_rise_rate = rel_sl_rise_rate
        self._tidal_range = tidal_range
        self._tidal_half_range = tidal_range / 2.0

    def get_water_depth(self, min_depth=0.01):
        """Calculate the water depth field."""

        depth_at_mean_high_water = np.amax(0, - self._elev
                                              + self._mean_sea_level
                                              + self._tidal_half_range)
        self._fully_wet_depth = (0.5 * (depth_at_mean_high_water
                                        + np.amax(0.0, depth_at_mean_high_water
                                                       - self._tidal_range)))
        self._fully_wet_depth[self._fully_wet_depth < min_depth] = min_depth
        self._water_depth[:] = self._fully_wet_depth
        self._water_depth[z > (self._mean_sea_level + self._tidal_half_range)] = 0.0

    def run_one_step(self, dt):
        """Advance in time."""

        # Update sea level
        self._mean_sea_level += self._rel_sl_rise_rate * dt

        # water depth

        # vegetation

        # roughness
