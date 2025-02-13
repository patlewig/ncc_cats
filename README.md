ncc_cats
==============================

A python package that codifies the current EPA New Chemicals Categories (NCC). These are a set of 56 chemical categories that EPA uses within its assessments under the New Chemicals Program assessment. More information about these categories can be found here (https://www.epa.gov/reviewing-new-chemicals-under-toxic-substances-control-act-tsca/chemical-categories-used-review-new).
The package herein allows for chemicals is be assigned in accordance with the NCC based on their chemical structure and selected physicochemical parameters. The NCC implemented herein used the NCC profiler from the OECD QSAR Toolbox v4.6 as a foundation. 

Install locally by using `pip install -e .` when in the root directory.

Project Organisation
------------

    ├── LICENSE
    ├── Makefile           <- Makefile with commands like `make data` or `make train`
    ├── README.md          <- The top-level README for developers using this project.
    ├── data
    │   ├── external       <- Data from third party sources.
    │   ├── interim        <- Intermediate data that has been transformed.
    │   ├── processed      <- The final, canonical data sets for modeling.
    │   └── raw            <- example data files
    │
    ├── docs               <- A default mkdocs project; see mkdocs.org for details
    │
    ├── notebooks          <- Jupyter notebooks. Naming convention is a number (for ordering),
    │                         and a short `-` delimited description, e.g.
    │                         `01-walkthrough`.
    │
    ├── references         <- Data dictionaries, manuals, and all other explanatory materials.
    │
    ├── reports            <- Generated analysis as HTML, PDF, LaTeX, etc.
    │   └── figures        <- Generated graphics and figures to be used in reporting
    │
    ├── requirements.txt   <- The requirements file for reproducing the analysis environment, e.g.
    │                         generated with `pip freeze > requirements.txt`
    │
    ├── setup.py           <- makes project pip installable (pip install -e .) so epa_ncc can be imported
    ├── epa_ncc                <- Source code for use in this project.
    │   ├── __init__.py    <- Makes src a Python module
    │   │
    │   ├── ncc_categories <- Scripts to process NCC
    │   │   └── ncc_categories.py
    │   │
    │   ├── data           <- Contains the xml file which encode some of the categories
    │   
    │
    └── tox.ini            <- tox file with settings for running tox; see tox.readthedocs.io


--------

<p><small>Project based on the <a target="_blank" href="https://drivendata.github.io/cookiecutter-data-science/">cookiecutter data science project template</a>. #cookiecutterdatascience</small></p>
