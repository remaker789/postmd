# sphinx-apidoc -M -e -f -t docs/_templates/apidoc  -o docs .
# sphinx-apidoc -M -e -f  -o docs .
sphinx-apidoc -M -f  -o docs .


# -e, --separate: Put documentation for each module on its own page.
# -T, --no-toc: Do not create a table of contents file. Ignored when --full is provided.
# -E, --no-headings: Do not create headings for the modules/packages. This is useful, for example, when docstrings already contain headings.
cd docs 
./make.bat clean
./make.bat html
cd ../