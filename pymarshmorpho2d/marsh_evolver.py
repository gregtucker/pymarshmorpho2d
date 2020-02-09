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
        "water_depth": {
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
        "veg_is_present": {
            "dtype": bool,
            "intent": "out",
            "optional": False,
            "units": "-",
            "mapping": "node",
            "doc": "True where marsh vegetation is present",
        },
        "vegetation": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "?",
            "mapping": "node",
            "doc": "Some measure of vegetation...?",
        },
        "roughness": {
            "dtype": float,
            "intent": "out",
            "optional": False,
            "units": "s/m^1/3",
            "mapping": "node",
            "doc": "Manning roughness coefficient",
        },
    }

    def __init__(self, grid,
                 rel_sl_rise_rate=2.74e-6,
                 tidal_range=3.1,
                 tidal_range_for_veg=3.1,
                 roughness_with_veg=0.1,
                 roughness_without_veg=0.02
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
        tidal_range_for_veg : float
            Tidal range for vegetation model, m (normally same as tidal range)
        """

        super(MarshEvolver, self).__init__(grid)
        self.initialize_output_fields()

        # Get references to fields
        self._elev = self._grid.at_node['topographic__elevation']
        self._water_depth = self._grid.at_node['water_depth']
        self._fully_wet_depth = self._grid.at_node['fully_wet__depth']
        self._veg_is_present = self._grid.at_node['veg_is_present']
        self._vegetation = self._grid.at_node['vegetation']
        self._roughness = self._grid.at_node['roughness']

        # Set parameter values
        self._mean_sea_level = 0.0
        self._rel_sl_rise_rate = rel_sl_rise_rate
        self._tidal_range = tidal_range
        self._tidal_range_for_veg = tidal_range_for_veg
        self._tidal_half_range = tidal_range / 2.0
        self._roughness_with_veg = roughness_with_veg
        self._roughness_without_veg = roughness_without_veg

        # lower and upper limits for veg growth [m]
        # see McKee, K.L., Patrick, W.H., Jr., 1988.
        # these are "dBlo" and "dBup" in matlab original
        self._min_elev_for_veg_growth = -(0.237 * self._tidal_range_for_veg
                                          - 0.092)
        self._max_elev_for_veg_growth = self._tidal_range_for_veg / 2.0

    def get_water_depth(self, min_depth=0.01):
        """Calculate the water depth field."""

        depth_at_mean_high_water = np.maximum((-self._elev) 
                                              + self._mean_sea_level
                                              + self._tidal_half_range, 0.0)
        self._fully_wet_depth = (0.5 * (depth_at_mean_high_water
                                        + np.maximum(depth_at_mean_high_water
                                                     - self._tidal_range, 0.0)))
        self._fully_wet_depth[self._fully_wet_depth < min_depth] = min_depth
        self._water_depth[:] = self._fully_wet_depth
        self._water_depth[self._elev 
                          > (self._mean_sea_level 
                             + self._tidal_half_range)] = 0.0

    def update_vegetation(self):
        """Update vegetation."""
        height_above_msl = self._elev - self._mean_sea_level
        self._veg_is_present[:] = (height_above_msl
                                   > self._min_elev_for_veg_growth)
        self._vegetation = (4 * (height_above_msl
                                 - self._max_elev_for_veg_growth)
                              * (self._min_elev_for_veg_growth
                                 - height_above_msl)
                              / (self._min_elev_for_veg_growth
                                 - self._max_elev_for_veg_growth)**2)
        self._vegetation[height_above_msl > self._max_elev_for_veg_growth] = 0.0
        self._vegetation[height_above_msl < self._min_elev_for_veg_growth] = 0.0
        
    def update_roughness(self):
        """Update Manning's n values."""
        self._roughness[:] = self._roughness_without_veg
        self._roughness[self._veg_is_present] = self._roughness_with_veg

    def run_one_step(self, dt):
        """Advance in time."""

        # Update sea level
        self._mean_sea_level += self._rel_sl_rise_rate * dt

        # water depth
        self.get_water_depth()

        # vegetation
        self.update_vegetation()

        # roughness
        self.update_roughness()