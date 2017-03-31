# University of Michigan Astro 461 -- Python Cookbook

Hello, and welcome to Astro 461! I've created this repository to serve as a
tutorial about some useful tools for astronomy analysis in Python for those of
you who have not used Python extensively thus far. I hope it proves to be
helpful! 

## Installation

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
concepts, such as image calibration, aperture photometry, astrometry, and
spectroscopy. In order to make sure these work out of the box, I've also
included a `requirements.txt` file which lists all of the modules you will need
for these notebooks. We will use python's `pip` to install them. First, upgrade
to the latest version of `pip`:

```bash
pip3 install --upgrade pip3
```

(if you prefer Python 2 instead, just replace all `pip3` occurrences with `pip`)

Then, install any missing packages with

```bash
pip3 install -r requirements.txt
```

You should now be ready to go!

## Usage

You should see a number of iPython(Jupyter) notebooks in this new directory, in
addition to a subdirectory called `Sampledata`. The notebooks make use of the
FITS files located within `Sampledata`, which were taking using instruments
you'll be using later in this class! 
