# University of Michigan Astro 461 -- Python Cookbook

Hello, and welcome to Astro 461! I've created this repository to serve as a
tutorial about some useful tools for astronomy analysis in Python for those of
you who have not used Python extensively thus far. I hope it proves to be
helpful! 

## Prerequisites

For the remainder of this file, I will assume you already have `anaconda` and
 `git` installed. If this is not the case, please follow the steps in the 
 `INSTALL.md` file before continuing!

## Setup

To start using the code in this repository, use your terminal to navigate to the
directory you'd like all of your Astro 461 work to be in. Then type 

```bash
git clone https://sjhamilton820@bitbucket.org/sjhamilton820/astro-461-python-cookbook.git
```

Then navigate to the directory that was just created:

```bash
cd astro-461-python-cookbook
```

You should see several Jupyter notebooks dealing with various important
concepts, such as image calibration and
spectroscopy. In order to make sure these work out of the box, I've also
included a `requirements.txt` file which lists all of the python modules you 
will need
for these notebooks. We will use python's `pip` to install any of the ones 
that may not have been included with Anaconda for some reason. 

```bash
pip install -r requirements.txt
```

You should now be ready to go!

## Usage

You should see a number of iPython(Jupyter) notebooks in this new directory, in
addition to a subdirectory called `Sampledata`. The notebooks make use of the
FITS files located within `Sampledata`, which were taking using instruments
you'll be using later in this class! 
