#!/bin/bash

set -e
set -x

cd efel/tests
nosetests -s -v -x
