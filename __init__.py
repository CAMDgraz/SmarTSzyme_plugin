# -*- coding: utf-8 -*-
"""
Created on Wed Nov 20 16:46:31 2024

@author: aliaa
"""

from pymol.plugins import addmenuitemqt
from . import test_load_plugin


# Function to launch the plugin
def launch_plugin():
    global plugin_window # need this or window crashes immediately after it is launched

    plugin_window = test_load_plugin.MainWindow()
    plugin_window.show()

# PyMOL plugin initialization
def __init_plugin__(app=None):
    # Add the plugin to the PyMOL menu
    addmenuitemqt("Aliaa GUI Plugin", launch_plugin)
