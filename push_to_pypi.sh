#!/usr/bin/env bash

rm -rf dist
rm -rf build
rm -rf redpatch.egg-info
python setup.py sdist bdist_wheel
python -m twine upload --repository-url https://upload.pypi.org/legacy/ dist/*