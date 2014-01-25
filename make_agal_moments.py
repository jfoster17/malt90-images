#!/usr/bin/env python
# encoding: utf-8
"""
make_agal_moments

Make moment maps using the catalog from Jill/Scott

Make a moment map for each ATLASGAL source. Simpler
than 

"""

import sys,os
import malt_params as malt
import numpy as np
from astropy.table import Table
