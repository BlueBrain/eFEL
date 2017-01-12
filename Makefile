all: install
install: clean
	python setup.py sdist
	pip install dist/*.tar.gz --upgrade
install3: clean
	python3 setup.py sdist
	pip3 install dist/*.tar.gz --upgrade
doc_efeatures:
	rm -rf docs/build_efeatures && \
	mkdir docs/build_efeatures && \
	cd docs/source/tex && \
	ls -al ../../build_efeatures && \
	pdflatex -output-directory=../../build_efeatures efeature-documentation.tex
doc: install doc_efeatures
	pip install sphinx sphinx-autobuild
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
update_version:
	cd efel && \
	python -c 'import version; version._get_version_number()' && \
	git add GITHASH.txt && \
	git add VERSION.txt && \
	git commit -m 'Updated version number'
test: install
	pip install nose coverage --upgrade
	cd efel/tests; nosetests -s -v -x --with-coverage --cover-xml \
		--cover-package efel || cat log/fllog.txt
	cat efel/tests/log/fllog.txt
test3: install3
	pip3 install nose coverage --upgrade
	cd efel/tests; nosetests-3.4 -s -v -x --with-coverage --cover-xml \
		--cover-package efel
pypi: test
	pip install twine --upgrade
	rm -rf dist
	python setup.py sdist bdist
	twine upload dist/*
clean:
	rm -rf build_cmake
	rm -rf build
	rm -rf docs/build
	rm -rf dist
cpp:
	mkdir -p build_cmake && \
	cd build_cmake && \
	cmake .. && \
	make -j
push: clean install test doc doc_upload
	git push
	git push --tags
graph:
	pip install pygraphviz==1.3.1
	utils/efel_graph_dependency -i efel/DependencyV5.txt --graph dependencies.png --graph-deps
