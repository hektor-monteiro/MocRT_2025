#!/bin/bash
set -e

# Reset conf.py and index.rst
echo "Resetting docs configuration..."
cat > docs/conf.py << EOL
# Configuration file for the Sphinx documentation builder.
#
# For the full list of built-in configuration values, see the documentation:
# https://www.sphinx-doc.org/en/master/usage/configuration.html

# -- Project information -----------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#project-information

project = 'MOCASSIN'
copyright = '2025, Jules'
author = 'Jules'
release = '1.0'

# -- General configuration ---------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#general-configuration

extensions = [ 'sphinxfortran.fortran_domain', 'sphinxfortran.fortran_autodoc' ]
fortran_src = []

templates_path = ['_templates']
exclude_patterns = ['_build', 'Thumbs.db', '.DS_Store']

language = 'en'

# -- Options for HTML output -------------------------------------------------
# https://www.sphinx-doc.org/en/master/usage/configuration.html#options-for-html-output

html_theme = 'alabaster'
html_static_path = ['_static']
EOL

cat > docs/index.rst << EOL
.. MOCASSIN documentation master file, created by
   sphinx-quickstart on Thu Jul 25 10:33:20 2024.
   You can adapt this file completely to your liking, but it should at least
   contain the root \`toctree\` directive.

Welcome to MOCASSIN's documentation!
===================================

.. toctree::
   :maxdepth: 2
   :caption: Contents:

EOL

# Find all .f90 files in the source directory
FILES=$(find source -name "*.f90")

# Iterate over each file, adding it to the documentation and trying to build
for file in $FILES; do
  module_name=$(basename "$file" .f90)
  echo "Testing $module_name..."

  # Add the file to conf.py
  sed -i "s|fortran_src = .*|fortran_src = [ '../source/$module_name.f90' ]|" docs/conf.py

  # Add the file to index.rst
  echo "   $module_name" >> docs/index.rst

  # Create the rst file
  echo "$module_name" > "docs/$module_name.rst"
  echo "==========" >> "docs/$module_name.rst"
  echo "" >> "docs/$module_name.rst"
  echo ".. f:automodule:: $module_name" >> "docs/$module_name.rst"


  # Try to build the documentation
  if ! (cd docs && make html > /dev/null 2>&1); then
    echo "Error building with $module_name"
    exit 1
  fi
done

echo "All files built successfully!"
