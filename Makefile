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
	cp ../../build_efeatures/efeature-documentation.pdf . \
	touch .nojekyll && \
	git init . && \
	git add . && \
	git commit -m "Updating docs" && \
	git push "git@github.com:BlueBrain/eFEL.git" master:gh-pages --force && \
	rm -rf .git
test: install
	cd efel/tests; nosetests -s -v -x
clean:
	rm -rf build
	rm -rf docs/build
