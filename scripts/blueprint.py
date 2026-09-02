#!/usr/bin/env python3

"""

Usage:
"""

# Standard library imports:
import sys
from pathlib import Path
import datetime

def main():
    pass

if __name__ == "__main__":

    print(f"--- SCRIPT: {__file__} ---")
    print(f"--- COMMAND: {' '.join(sys.argv)} ---")
    print(f"--- START: {datetime.datetime.now()} ---")

    main()
    
    print(f"--- END: {datetime.datetime.now()} ---")
