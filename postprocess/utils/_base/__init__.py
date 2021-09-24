#!/usr/bin/env python
# -*- coding: utf-8 -*-
import json
import os

from .spin import Spin
from .latt import Latt

_fn = os.path.join(
    os.path.dirname(__file__), "elements.json"
)
with open(_fn, "r") as ele:
    periodic_table = json.load(ele)

order_table = {v: k for k, v in periodic_table.items()}

if __name__ == "__main__":
    pass
