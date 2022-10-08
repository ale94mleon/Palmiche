Installation
============

This package is not release on PyPi, but could be easily installed via pip directly from the `GitHub Repo <https://github.com/ale94mleon/palmiche/>`_.

Via pip (standard)
------------------

In this case you must have a correct installation
ofopenbable and autodock-vina. If you already have it, skip the creation of the conda environment and the conda dependencies installation:

Create conda environment and install conda dependencies
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. code-block:: bash

    conda create -n palmiche
    conda activate palmiche

Then install the dependencies libraries:

.. code-block:: bash

    conda install -y -c conda-forge vina openbabel

..  In the future we will consider to use the python modules `vina on pypi <https://pypi.org/project/vina/>`_. Finally:

pip install
~~~~~~~~~~~

.. code-block:: bash

    # To get the version on development (only posable alternative for now)
    pip install git+https://github.com/ale94mleon/palmiche.git@main
