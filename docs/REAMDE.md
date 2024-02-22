## Sphinx
In the root dir, run the `sphinx-apidoc` command to produce the rst files of packages
```shell
sphinx-apidoc -o docs .
```
Remember to add `modules.rst` in the `index.rst`, like:
```index.rst
.. toctree::
   :maxdepth: 2
   :caption: Contents:

   modules
```

Then, enter `docs` dir, and make the format you want.
```shell
cd docs
make html # here make html
```