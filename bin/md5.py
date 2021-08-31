#!/usr/bin/env python
# Jianhong Ou @ 2021 @ duke

import hashlib
import sys

inputArgs = sys.argv

for i in inputArgs[1:]:
    with open(i, 'rb') as f:
        fh = hashlib.md5()
        while chunk := f.read(8192):
            fh.update(chunk)
        print(fh.hexdigest(), "\t", i, "\n")
