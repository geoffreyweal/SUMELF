# Installation

In this article, we will look at how to install the SUMELF program and all requisites required for this program.

## Pre-requisites

### Python 3 and ``pip3``

This program is designed to work with **Python 3**. This program can only be used with Python 3.7. This is because the CSD Python API can only run using Python 3.7. 

To find out if you have Python 3 on your computer and what version you have, type into the terminal

```bash
python --version
```

If you have Python 3 on your computer, you will get the version of python you have on your computer. E.g.

```bash
user@computer_name path % python --version
Python 3.7.9
```

If you have Python 3, you may have ``pip`` installed on your computer as well. ``pip`` is a python package installation tool that is recommended by Python for installing Python packages. To see if you have ``pip`` installed, type into the terminal

```bash
pip list
```
If you get back a list of python packages install on your computer, you have ``pip`` installed. E.g.

```bash
user@computer_name Documentation % pip3 list
Package                       Version
----------------------------- ---------
alabaster                     0.7.12
asap3                         3.11.10
ase                           3.20.1
Babel                         2.8.0
certifi                       2020.6.20
chardet                       3.0.4
click                         7.1.2
cycler                        0.10.0
docutils                      0.16
Flask                         1.1.2
idna                          2.10
imagesize                     1.2.0
itsdangerous                  1.1.0
Jinja2                        2.11.2
kiwisolver                    1.2.0
MarkupSafe                    1.1.1
matplotlib                    3.3.1
numpy                         1.19.1
packaging                     20.4
Pillow                        7.2.0
pip                           20.2.4
Pygments                      2.7.1
pyparsing                     2.4.7
python-dateutil               2.8.1
pytz                          2020.1
requests                      2.24.0
scipy                         1.5.2
setuptools                    41.2.0
six                           1.15.0
snowballstemmer               2.0.0
Sphinx                        3.2.1
sphinx-pyreverse              0.0.13
sphinx-rtd-theme              0.5.0
sphinx-tabs                   1.3.0
sphinxcontrib-applehelp       1.0.2
sphinxcontrib-devhelp         1.0.2
sphinxcontrib-htmlhelp        1.0.3
sphinxcontrib-jsmath          1.0.1
sphinxcontrib-plantuml        0.18.1
sphinxcontrib-qthelp          1.0.3
sphinxcontrib-serializinghtml 1.1.4
sphinxcontrib-websupport      1.2.4
urllib3                       1.25.10
Werkzeug                      1.0.1
wheel                         0.33.1
xlrd                          1.2.0
```

If you do not see this, you probably do not have ``pip`` installed on your computer. If this is the case, check out [PIP Installation](https://pip.pypa.io/en/stable/installing). 

!!! note

	In most cases, ``pip`` and ``pip3`` are synonymous for the Python Installation Package for Python 3. **However in some cases,** ``pip`` **will be directed to the Python Installation Package for Python 2 rather than Python 3.** To check this, run in the terminal:

	```bash
	pip --version
	```

	If the output indicates you this Python Installation Package is for Python 2 and not Python 3, only install packages using the ``pip3`` name. 

	For the rest of this documentation, ``pip`` will be used, however if your computer's ``pip``  refers to Python 2 and not Python 3, use ``pip3``  instead of ``pip``. 


### Atomic Simulation Environment (ASE)

The SUMELF program uses the Atomic Simulation Environment (ASE) to read in the crystal data from various file format, to process the crystals, and to save the the crystals to disk. Read more about [ASE here](https://wiki.fysik.dtu.dk/ase). 

The installation of ASE can be found on the [ASE installation page](https://wiki.fysik.dtu.dk/ase/install.html), however from experience if you are using ASE for the first time, it is best to install ASE using ``pip``, the package manager that is an extension of python to keep all your program easily managed and easy to import into your python. 

To install ASE using ``pip``, perform the following in your terminal.

```bash
pip install --upgrade --user ase
```

Installing using ``pip`` ensures that ASE is being installed to be used by Python 3, and not Python 2. Installing ASE like this will also install all the requisite program needed for ASE. This installation includes the use of features such as viewing the xyz files of structure and looking at ase databases through a website. These should be already assessible, which you can test by entering into the terminal:

```bash
ase gui
```

This should show a GUI with nothing in it, as shown below.

<figure markdown="span">
	<img alt="ase_gui_blank.png" src="Shared_Images/ase_gui_blank.png?raw=true" width="500" />
  <!-- ![ase_gui_blank.png](Images/ase_gui_blank.png){ width="500" } -->
  <figcaption>This is a blank ase gui screen that you would see if enter ase gui into the terminal.</figcaption>
</figure>

However, **in the case that this does not work**, we need to manually add a path to your ``~/.bashrc`` so you can use the ASE features externally outside python. Do the following; first enter the following into the terminal:

```bash
pip show ase
```

This will give a bunch of information, including the location of ase on your computer. For example, when I do this I get:

```bash
user@computer_name docs % pip show ase
Name: ase
Version: 3.20.1
Summary: Atomic Simulation Environment
Home-page: https://wiki.fysik.dtu.dk/ase
Author: None
Author-email: None
License: LGPLv2.1+
Location: /Users/geoffreyweal/Library/Python/3.7/lib/python/site-packages
Requires: matplotlib, scipy, numpy
Required-by: 
```

Copy the ``Location`` line. If we remove the ``lib/python/site-packages`` bit and replace it with ``bin``, this gives us the location of useful ASE programs. The example below is for Python 3.7. 

```bash
/Users/geoffreyweal/Library/Python/3.7/bin
```

Next, add this to your ``~/.bashrc`` file as below:

```bash
############################################################
# For ASE
export PATH=/Users/geoffreyweal/Library/Python/3.7/bin:$PATH
############################################################
```

Write ``source ~/.bashrc`` in the terminal and press enter. Once you have done this, try to run ``ase gui`` in the terminal. This will hopefully show the ase gui and allow you to access the useful ASE programs through the terminal. 


### Networkx

``Networkx`` is a python program that is used in the SUMELF program to describe the bonding structure between atoms in the crystal structure. The easiest way to install ``Networkx`` is though ``pip``. Type the following into the terminal:

```bash
pip3 install --upgrade --user networkx
```

### Packaging

The ``packaging`` program is also used in this program to check the versions of ASE that you are using for compatibility issues. The easiest way to install ``packaging`` is though ``pip``. Type the following into the terminal:

```bash
pip3 install --upgrade --user packaging
```

### Pymatgen

``Pymatgen`` is a python program that is used in SUMELF to determine symmetric molecules within a crystal structure. The easiest way to install ``Pymatgen`` is though ``pip``. Type the following into the terminal:

```bash
pip3 install --upgrade --user pymatgen
```

This package and other required packages may take a bit of time to install. 

### TQDM

The ``tqdm`` program is used by this program to provide progress bars that are useful for easily monitoring progress during this program. The easiest way to install ``tqdm`` is though ``pip``. Type the following into the terminal:

```bash
pip3 install --upgrade --user tqdm
```

### Xlsxwriter

The ``xlsxwriter`` program is used by this program to write the output data from Gaussian jobs to an excel file(s). The easiest way to install ``xlsxwriter`` is though ``pip``. Type the following into the terminal:

```bash
pip3 install --upgrade --user xlsxwriter
```

### Memory Profiler

The ``memory_profiler`` program is a really useful program for determining how much memory is being used by a python program or any program using python. It is used here to determine how much memory is being used to process matrix data from ``output.log`` files. The easiest way to install ``memory_profiler`` is though ``pip``. Type the following into the terminal:

```bash
pip3 install --upgrade --user memory_profiler
```


## Setting up the SUMELF Program

There are three ways to install SUMELF on your system. These ways are described below:


### Install SUMELF through ``pip3``

To install the SUMELF program using ``pip3``, perform the following in your terminal:

```bash
pip3 install --upgrade --user git+https://github.com/geoffreyweal/SUMELF.git
```

To install all related programs at the same time, performing the following in your terminal: 

```bash
pip3 install --upgrade --user git+https://github.com/geoffreyweal/SUMELF.git git+https://github.com/geoffreyweal/ACSD.git git+https://github.com/geoffreyweal/ReCrystals.git git+https://github.com/geoffreyweal/RSGC.git git+https://github.com/geoffreyweal/ReJig.git git+https://github.com/geoffreyweal/ECCP.git git+https://github.com/geoffreyweal/EKMC.git git+https://github.com/geoffreyweal/SORE.git
```


### Install SUMELF through ``conda``

You can install the SUMELF program on ``conda`` through ``pip``. [Click here](https://conda.io/projects/conda/en/latest/user-guide/tasks/manage-pkgs.html#installing-non-conda-packages) to see more information about installing SUMELF through ``conda``.


### Manual installation

First, download the SUMELF program to your computer. You can do this by cloning a version of this from Github, or obtaining a version of the program from the authors. If you are obtaining this program via Github, you want to ``cd`` to the directory that you want to place this program in on the terminal, and then clone the program from Github through the terminal as well: 

```bash
cd PATH/TO/WHERE_YOU_WANT_TO_PLACE_SUMELF_ON_YOUR_COMPUTER
git clone https://github.com/geoffreyweal/SUMELF
```

Second, you need to change permissions to use this program by using ``chmod``. In the terminal write:

```bash
chmod -R 777 SUMELF
```

Third, you will want to add a python path for the ASCD program to the ``~/.bashrc`` to indicate its location. You can do this by entering the following into the terminal and clicking enter:

```bash
echo '
###################################
# Used for the SUMELF Program
export PATH_TO_SUMELF="'$PWD'/SUMELF" 
export PYTHONPATH="$PATH_TO_SUMELF":$PYTHONPATH
export PATH="$PATH_TO_SUMELF"/bin:$PATH
###################################
' >> ~/.bashrc
```

You can check that this has been entered into your ``~/.bashrc`` file by typing ``vim ~/.bashrc`` into the terminal, and scrolling down to the bottom of the terminal. 

!!! tip

	Make sure that the path given to ``PATH_TO_SUMELF`` is the correct path to the SUMELF folder. 

Finally, source your ``~/.bashrc`` file by typing the following into the terminal and pressing the enter button:

```bash
source ~/.bashrc
```

Once you have run ``source ~/.bashrc``, the SUMELF program should be all ready to go! You can check this by typing the following into the terminal:

```bash
which SUMELF
```

This should give you the path to the SUMELF program. If the terminal tells you it can not find this program, check that the path you gave for ``PATH_TO_SUMELF`` is the correct path to the SUMELF folder. 

#### Summary of ``~/.bashrc`` input

You want to have the following in your ``~/.bashrc``:

```bash
###################################
# Used for the SUMELF Program
export PATH_TO_SUMELF="<Path_to_SUMELF>" 
export PYTHONPATH="$PATH_TO_SUMELF":$PYTHONPATH
export PATH="$PATH_TO_SUMELF"/bin:$PATH
###################################
```

where ``"<Path_to_SUMELF>"`` is the directory path that you place the SUMELF program. You can find this by changing directory (``cd``) into the SUMELF folder and typing ``pwd`` into the terminal. This will give you the full path to the SUMELF program. 


## Other Useful Commands

There are several commands that are useful for running all the programs in the grand scheme for simulating/calculating exciton and charge diffusion. 

### The ``qme`` command

You may use ``squeue`` to figure out what jobs are running in slurm. For monitoring what slurm jobs are running, I have found the following alias useful. To include this in your ``~/.bashrc`` file, type the following into the terminal and click enter: 

```bash
echo "
alias qme='squeue -o \"%.18i %.7P %.5Q %.70j %.8u %.8T %.10M %.11l %.6D %.4C %.6m %.15S %.10R %.8q\" -u $USER --sort=+i'
" >> ~/.bashrc
```

To run this, type ``qme`` into your computer. 

### The ``no_of_jobs_running_or_queued`` command

The ``no_of_jobs_running_or_queued`` command is designed to determine the number of jobs that are either running or in the queue in slurm. 

Add this to your ``~/.bashrc`` file by typing the following into the terminal and pressing enter. 

```bash
echo "
alias no_of_jobs_running_or_queued=\"squeue -u $USER | wc -l\"
" >> ~/.bashrc
```



