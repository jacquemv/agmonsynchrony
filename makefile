all:
	python -m build #--no-isolation

local:
	python setup.py build_ext -i

sdist:
	python setup.py sdist

check: sdist
	twine check dist/*

upload: sdist
	twine upload dist/*

install:
	ln -s $(PWD)/agmonsynchrony $(HOME)/python

clean:
	rm -rf build dist *.egg-info