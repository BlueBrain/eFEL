TEST_REQUIREMENTS=nose coverage virtualenv

all: install
install: clean
	python setup.py sdist
	pip install `ls dist/efel-*.tar.gz`[neo] --upgrade
virtualenv: clean
	virtualenv pyenv
	. ./pyenv/bin/activate	
doc_efeatures:
	rm -rf docs/build_efeatures && \
	mkdir docs/build_efeatures && \
	cd docs/source/tex && \
	ls -al ../../build_efeatures && \
	pdflatex -output-directory=../../build_efeatures efeature-documentation.tex
doc: install doc_efeatures
	pip install sphinx sphinx-autobuild sphinx_rtd_theme
	cd docs; $(MAKE) clean; $(MAKE) html
doc_upload: doc
	cd docs/build/html && \
	cp ../../build_efeatures/efeature-documentation.pdf . && \
	touch .nojekyll && \
	git init . && \
	git add . && \
	git commit -m "Updating docs" && \
	git push "git@github.com:BlueBrain/eFEL.git" master:gh-pages --force && \
	rm -rf .git
install_test_requirements:
	pip install -q $(TEST_REQUIREMENTS) --upgrade
update_version:
	cd efel && \
	python -c 'import version; version._get_version_number()' && \
	git add GITHASH.txt && \
	git add VERSION.txt && \
	git commit -m 'Updated version number'
test: virtualenv install install_test_requirements
	cd efel/tests; nosetests -s -v -x --with-coverage --cover-xml \
	   --cover-package efel 
debugtest: virtualenv install install_test_requirements
	cd efel/tests; nosetests -a debugtest -s -v -x --with-coverage --cover-xml \
	   --cover-package efel 
pypi: test
	pip install twine --upgrade
	rm -rf dist
	python setup.py sdist bdist
	twine upload dist/*
clean:
	rm -rf efel/tests/log/fllog.txt
	rm -rf build_cmake
	rm -rf build
	rm -rf docs/build
	rm -rf dist
	rm -rf pyenv
cpp:
	mkdir -p build_cmake && \
	cd build_cmake && \
	cmake .. && \
	make -j
push: clean install test doc doc_upload
	git push
	git push --tags
format:	
	clang-format -i -style="google" efel/cppcore/*.cpp
graph:
	pip install pygraphviz==1.3.1
	utils/efel_graph_dependency -i efel/DependencyV5.txt --graph dependencies.png --graph-deps
