python setup.py bdist_wheel --universal
twine upload --repository-url https://test.pypi.org/legacy/ dist/*
pip uninstall -y memdof
pip install -i https://test.pypi.org/simple/ memdof