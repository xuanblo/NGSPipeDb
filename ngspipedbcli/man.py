#!/usr/bin/env python3
# -*- coding: utf-8 -*-

################################
# rich ngspipedb description
################################

from rich import print
from rich.console import Console
from rich.table import Column, Table
from rich.markdown import Markdown

console = Console()
table = Table(show_header=True, header_style="bold magenta")
table.add_column("Module", style="dim", width=12)
table.add_column("Description")
table.add_column("document", justify="right")
table.add_column("Version", justify="right")
table.add_row(
    "Dev 20, 2019", "www.liu-lab.com", "$275,000,000", "$375,126,118"
)
table.add_row(
    "May 25, 2018",
    "[red]Solo[/red]: A Star Wars Story",
    "$275,000,000",
    "$393,151,347",
)
table.add_row(
    "Dec 15, 2017",
    "Star Wars Ep. VIII: The Last Jedi",
    "$262,000,000",
    "[bold]$1,332,539,889[/bold]",
)

def man(args):
    print('man, {0}!'.format(args.name))

def doc_download(args):
    MARKDOWN = """
    # Download

    Rich can do a pretty *decent* job of rendering markdown.

    1. This is a list [item](./xx.fa)
    2. This is another list item
    """
    md = Markdown(MARKDOWN)
    console.print(md)