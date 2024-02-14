#!/bin/sh

set -e

INSTALL_DIR=$1
MOD_DIR=$2

cd ${INSTALL_DIR}

echo "Building mod files"
rm -rf x86_64
nrnivmodl ${MOD_DIR} >nrnivmodl.log 2>&1
