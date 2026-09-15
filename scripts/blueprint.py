#!/usr/bin/env python3

# Standard library imports:
import sys
from pathlib import Path
import datetime
import argparse

def createParser():
    parser = argparse.ArgumentParser(prog='<prog_name>',
                                     description='',
                                     epilog='Created by mimitair.')

    #parser.add_argument('path_to_cif_collection', action='store', type=Path, help='Path to the cif collection')
    # add more args here
    
    return parser


def main():
    args = createParser().parse_args()

    # logic

if __name__ == "__main__":

    start = datetime.datetime.now()
    print(f"--- SCRIPT: {__file__} ---")
    print(f"--- COMMAND: {' '.join(sys.argv)} ---")
    print(f"--- START: {start} ---")

    main()
    
    end = datetime.datetime.now()
    print(f"--- END: {end} ---")
    print(f"--- TIME ELAPSED: {end - start} ---")
