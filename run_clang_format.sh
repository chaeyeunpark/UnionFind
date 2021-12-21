#!/usr/bin/env bash

SCRIPT_DIR="$( cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd )"
FILE_ARR=( $("${SCRIPT_DIR}/build_utils/yavque_files.py" "include_tests" | jq -r '.[] | .name') )

clang-format-12 -i ${FILE_ARR[@]/%/}
