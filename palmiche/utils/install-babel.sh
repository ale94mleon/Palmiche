#!/usr/bin/env bash

set -e

SCRIPT_DIR=$(cd -- "$( dirname -- "${BASH_SOURCE[0]}" )" &> /dev/null && pwd)

NAME="babel"

(
  cd "${SCRIPT_DIR}"

  git clone --depth 1 --branch openbabel-3-1-1 https://github.com/openbabel/openbabel.git "${NAME}-src"

  PATH="${SCRIPT_DIR}/cmake/bin:${PATH}"

  cmake -S "${NAME}-src" -B "${NAME}-build" -DCMAKE_INSTALL_PREFIX="$(pwd)/${NAME}"
  
  ( cd "${NAME}-build" && make -j 8 && make install )

  rm -rf "${NAME}-build"
  rm -rf "${NAME}-src"
)
