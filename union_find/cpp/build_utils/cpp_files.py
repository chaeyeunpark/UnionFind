#!/usr/bin/env python3
from pathlib import Path
import json
import sys
import os

FILE_EXTENSIONS = ['.hpp', '.cpp']
SOURCE_DIRS = ['binding',  'include', 'src']
EXCLUDE_DIRS = []
PROJECT_SOURCE_DIR = Path(__file__).parent.parent.resolve()

if __name__ == '__main__':
    include_examples = False
    include_tests = False

    for arg in sys.argv[1:]:
        if arg == '--include-examples':
            SOURCE_DIRS += ['examples']
        if arg == '--include-tests':
            SOURCE_DIRS += ['tests']


    file_list = []
    for source_dir in SOURCE_DIRS:
        source_dir_full = PROJECT_SOURCE_DIR.joinpath(source_dir)
        for ext in FILE_EXTENSIONS:
            file_list += [p for p in source_dir_full.rglob('*' + ext)]

    file_list = set(file_list)

    for exclude_dir in EXCLUDE_DIRS:
        exclude_dir_full = PROJECT_SOURCE_DIR.joinpath(exclude_dir)
        for ext in FILE_EXTENSIONS:
            for exclude_file in exclude_dir_full.rglob('*' + ext):
                file_list.discard(exclude_file)

    file_list = [{'name': str(v)} for v in file_list]
    json.dump(file_list, sys.stdout)
