python setup.py bdist_wheel --universal
twine upload --repository-url https://test.pypi.org/legacy/ dist/*