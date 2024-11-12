Installation
=========================

CompactObject is an open-source tool designed for comprehensive neutron star equation inference. It is built to be easy to install and use. Follow the step-by-step installation guide below to get started.

Using Python Virtual Environment (Recommended)
------------------------------------------

If you are not using Anaconda, you can create a virtual environment using Python's `venv` module:

1. **Create a Virtual Environment**

   Run the following command to create a virtual environment named `CompactObject`:

   .. code-block:: bash

       python3 -m venv CompactObject

   *You can specify a different path by replacing `CompactObject` with your desired directory name.*

2. **Activate the Environment**

   Activate the virtual environment with:

   .. code-block:: bash

       source CompactObject/bin/activate

3. **Install the CompactObject package**

    Once the environment is activated, you can install the package from PyPI, which is called `CompactObject-TOV`:

    .. code-block:: bash
    
         pip install CompactObject-TOV
   
    The dependencies should be automatically installed for you. 
   
    while if you have trouble to use pip install, maybe you could git clone this repository, Then

   .. code-block:: bash
    
         pip install -r requirements.txt

**Alternative**:

    You can install the CompactObject package using pip inside the repository directory:
    
   .. code-block:: bash
    
         pip install -e .

    To upgrade to the latest version on PyPI, run:

   .. code-block:: bash

       pip install CompactObject-TOV --upgrade


5. **Using the Package**

   You are now ready to use CompactObject. Each time you want to use the package, ensure you activate the environment:

   .. code-block:: bash

       source CompactObject/bin/activate


Alternative: Using Anaconda
---------------------------

1. **Create a Virtual Environment**

   You can also create a virtual environment for CompactObject using Anaconda:

   .. code-block:: bash

       conda create -n CompactObject

   When prompted to proceed, type `y` and press Enter:

   .. code-block:: none

       Proceed ([y]/n)?

2. **Activate the Environment**

   Activate your newly created environment with the following command:

   .. code-block:: bash

       conda activate CompactObject

   **Note:** Once you create this environment, you don't need to create it again. Simply activate it whenever you want to use CompactObject.

3. **Install the CompactObject package**

    Once the environment is activated, you can install the package from PyPI, which is called `CompactObject-TOV`:

    .. code-block:: bash
    
         pip install CompactObject-TOV
   
    The dependencies should be automatically installed for you. 
   
    while if you have trouble to use pip install, maybe you could git clone this repository, Then

   .. code-block:: bash
    
         pip install -r requirements.txt

**Alternative**:

    You can install the CompactObject package using pip inside the repository directory:
    
   .. code-block:: bash
    
         pip install -e .

    To upgrade to the latest version on PyPI, run:

   .. code-block:: bash

       pip install CompactObject-TOV --upgrade

5. **You're Ready to Use CompactObject!**

   You are now ready to use CompactObject. Whenever you want to use this package, remember to activate the environment first:

   .. code-block:: bash

       conda activate CompactObject

..    Our package automatically installs all necessary dependencies for you. The dependencies include:

..    - `corner`
..    - `csv`
..    - `itertools`
..    - `math`
..    - `matplotlib`
..    - `numba`
..    - `numbaminpack`
..    - `numpy`
..    - `os`
..    - `pandas`
..    - `scipy`
..    - `sys`
..    - `ultranest`

`CompactObject-TOV` optionally depends on `numbaminpack`. However, it may hard to install if don't have a fortran complier. Please
refer to this page `NumbaMinpack documentation <https://pypi.org/project/NumbaMinpack/>`_, and you can skip
this dependency if you are not using "fastRMF_EoS" and "pQCD"

.. Summary
.. -------

.. - **Using Anaconda:**
..   1. Create and activate the `CompactObject` environment.
..   2. Install CompactObject with `pip`.
..   3. Activate the environment whenever you use the package.

.. - **Using Python Virtual Environment:**
..   1. Create and activate the `CompactObject` virtual environment.
..   2. Install CompactObject with `pip`.
..   3. Activate the environment whenever you use the package.

If you encounter any issues or have questions, feel free to reach out for support. Happy computing!
