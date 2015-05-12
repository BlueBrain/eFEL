all: install	
install:
	pip install . --upgrade
doc: install
	cd docs; $(MAKE) clean; $(MAKE) html
doc_upload: doc
	git checkout gh-pages                                                            
	rm -rf build _sources _static _modules efel
	git reset HEAD
	mv -fv docs/build/html/* ./
	rm -rf docs/build
	git add -A                                                                       
	git commit -m "Generated github docs for `git log master -1 --pretty=short --abbrev-commit`" && git push origin gh-pages ; git checkout master
test:
	cd efel/tests; nosetests -s -v -x
clean:
	rm -rf build
	rm -rf docs/build
