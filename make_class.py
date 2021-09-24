#!/usr/env/python

import os
import os.path as osp
from monty.os import cd
from pathlib import Path
import sys


def _yield(path):
   yield from (osp.join(path, i) for i in os.listdir(path))


def _mk_cls(path, cls_num):
    lst = list(_yield(path))
    clst = [lst[i:i + cls_num] for i in range(0, len(lst), cls_num)]
    with cd(path):
        for idx, ist in enumerate(clst):
            cls_name = f"class_{idx + 1}"
            print(cls_name)
            os.makedirs(cls_name)
            for stru in ist:
                _, name = osp.split(stru)
                print(name)
                os.system("mv {} {}/".format(name, cls_name))

if __name__ == "__main__":
    args = sys.argv

    root = args[1]
    _mk_cls(root, 70)
