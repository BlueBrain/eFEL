all: install	
install:
	pip install . --upgrade
doc_efeatures:
	rm -rf docs/build_efeatures && \
	mkdir docs/build_efeatures && \
	cd docs/source/tex && \
	ls -al ../../build_efeatures && \
	pdflatex -output-directory=../../build_efeatures efeature-documentation.tex
doc: install doc_efeatures
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
test: install
	cd efel/tests; nosetests -s -v -x
pypi: test
	rm -rf dist
	python setup.py sdist bdist
	twine upload dist/*
clean:
	rm -rf build_cmake
	rm -rf build
	rm -rf docs/build
cpp:
	mkdir -p build_cmake && \
	cd build_cmake && \
	cmake .. && \
	make -j
push: clean install test doc doc_upload
	git push
