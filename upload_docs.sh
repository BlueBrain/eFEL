#!/bin/bash

set -e
set -x

SOURCES='doc/source efel doc/Makefile'

git checkout gh-pages
rm -rf build _sources _static _modules
git checkout master $(SOURCES)
git reset HEAD
make -f doc/Makefile html
mv -fv build/html/* ./
rm -rf $(SOURCES) build
git add -A
git commit -m "Generated github docs for `git log master -1 --pretty=short \
--abbrev-commit`" && git push origin gh-pages ; git checkout master
