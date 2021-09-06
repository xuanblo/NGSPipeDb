#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
2021
"""

import os
import sys

if not __package__:
    path = os.path.join(os.path.dirname(__file__), os.pardir)
    sys.path.insert(0, path)

def main():
    from ngspipedbcli.app_click import ngspipedb_cli_main
    ngspipedb_cli_main()

if __name__ == '__main__':
    sys.exit(main())