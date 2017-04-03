# Installing Prerequisite Software

You will find both Anaconda and Git to be very useful for this tutorial 
repository and for Astro 461 in general. Below, I detail how to install each 
them if you still need to do so.

## Anaconda

1. Download anaconda from [here](https://www.continuum.io/downloads). Whether
   you use Python 2 or 3 is up to you, though we recommend Python 3.
    * Mac OSX users should choose the command line installer
2. AstroConda is a package repository maintained by the Space 
   Telescope Institute, and contains all of the packages and software (not just 
   python!) that one may find useful for astronomy analysis. Add the 
   AstroConda channel to anaconda by typing at your terminal:
   
   ```bash
   conda config --add channels http://ssb.stsci.edu/astroconda
   ```
   
3. You now have a few choices in terms of which AstroConda software to 
   install. The options are detailed 
   [here](https://astroconda.readthedocs.io/en/latest/installation.html). I've 
   opted for the "Legacy Software Stack," which includes most of the useful 
   analysis packages, in addition to IRAF and PyRAF. Install these by typing 
   at your terminal:
   
   ```bash
   conda install stsci pyraf iraf python=3
   ```
   
   The `python=3` at the end simply ensures that Anaconda installs the 
   compatible versions of everything.
    * You may see a wall of text after entering this command about new 
      packages being installed, packages being updated/downgraded, etc. This is
      all okay, and is Anaconda ensuring everything is compatible with 
      everything else. 
      
4. The last piece of software that will be useful for us is Source Extractor 
   (hopefully the name is self-explanatory!). Install it with:
   
   ```bash
   conda install sextractor
   ```
   
## Git

Since there are different sets of directions for Linux/Mac/Windows users, I 
will simply direct you to the 
[Git installation page](https://git-scm.com/book/en/v2/Getting-Started-Installing-Git), where installation for 
each platform is explicitly detailed.