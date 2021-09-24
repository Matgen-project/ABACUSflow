#!/usr/bin/env python
# -*- coding: utf-8 -*-

from enum import unique, Enum


@unique
class Spin(Enum):
    """
    Enum type for Spin.  Only up and down.
    Usage: Spin.up, Spin.down.
    """

    up, down = (1, -1)

    def __int__(self):
        return self.value

    def __float__(self):
        return float(self.value)

    def __str__(self):
        return str(self.value)


if __name__ == "__main__":
    pass
