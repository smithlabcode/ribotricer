#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Part of ribotricer software
#
# Copyright (C) 2020 Saket Choudhary, Wenzheng Li, and Andrew D Smith
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

"""
Minimal setup.py for backward compatibility.

All package configuration is now in pyproject.toml.
This file is kept for compatibility with older pip versions
and editable installs.
"""

import setuptools

if __name__ == "__main__":
    setuptools.setup()
