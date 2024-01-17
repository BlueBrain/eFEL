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
	pip install -r requirements_docs.txt
	cd docs; $(MAKE) clean; $(MAKE) html SPHINXOPTS=-W
doc_upload: doc
	cd docs/build/html && \
	cp ../../build_efeatures/efeature-documentation.pdf . && \
	touch .nojekyll && \
	git init . && \
	git add . && \
	git commit -m "Updating docs" && \
	git push "git@github.com:BlueBrain/eFEL.git" master:gh-pages --force && \
	rm -rf .git
update_version:
	cd efel && \
	python -c 'import version; version._get_version_number()' && \
	git add GITHASH.txt && \
	git add VERSION.txt && \
	git commit -m 'Updated version number'
pypi: test
	pip install twine --upgrade
	rm -rf dist
	python setup.py sdist bdist
	twine upload dist/*
clean:
	rm -f tests/log/fllog.txt
	rm -f fllog.txt
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
	pip install pygraphviz==1.11
	python utils/efel_graph_dependency.py -i efel/DependencyV5.txt --graph dependencies.png --graph-deps
