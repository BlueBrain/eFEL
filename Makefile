all: install	
install:
	pip install . --upgrade
doc: install
	cd docs; $(MAKE) clean; $(MAKE) html
doc_upload: doc
	cd docs/build/html && \
	git init . && \
	git add . && \
	git commit -m "Updating docs" && \
	git push "git@github.com:BlueBrain/eFEL.git" master:gh-pages --force && \
	rm -rf .git
test:
	cd efel/tests; nosetests -s -v -x
clean:
	rm -rf build
	rm -rf docs/build
