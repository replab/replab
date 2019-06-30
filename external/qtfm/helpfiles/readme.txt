Description of help files stored in this directory and how to process them
--------------------------------------------------------------------------

The files in this directory include help text files and the master XML
files used to create them. This file describes how this system of files
works and how to edit/expand the help information for QTFM.

1. How does Matlab access the help information to display help about QTFM?

In Matlab's own documentation, see the page called 'Display Custom
Documentation' (R2016a), under 'Toolbox Distribution'. This describes how a
custom help system can be added to Matlab. There are two key files, both of
which are implemented in QTFM:

a. info.xml     - This files resides in the main toolbox folder, and Matlab
                  finds it because that folder is on the path. The file
                  tells Matlab the name of the toolbox and where to find
                  the help information. In the case of QTFM this is the
                  folder called 'helpfiles' which is a sub-folder of the
                  main folder 'qtfm'.

b. helptoc.xml  - This file resides in the helpfiles folder, and it tells
                  Matlab the main entry points to the help information. It
                  is used to construct a table of contents which appears in
                  the help browser when the QTFM help information is
                  selected.

Both files are fairly simple XML format, and should it should be obvious
how they work, but if not, the Matlab documentation 'Display Custom
Documentation' explains the detailed format of both files.

2. How are the help files created/edited/processed?

The online help is stored in HTML files, and it is these files that are
displayed by the Matlab help browser (including graphics files which are
linked in the usual HTML manner). Obviously, hyperlinks can be included to
other help files, or external URLs.

Most of the HTML files are created from master XML files (i.e. the HTML
files are not themselves edited, it is the XML files that are edited). The
exceptions to this are mainly listed in the helptoc.xml file, but a full
list is given here:

a. qtfm.html         - This file is the top-level help page.
b. overview.html     - This page provides an overview of quaternions.
c. alphabetical.html - This is a list of toolbox functions, alphabetically
                       ordered.
d. categorical.html  - This is a list of toolbox functions, ordered into
                       categories.
e. matlist.html      - This is a list of built-in Matlab functions that
                       work with the toolbox. There is some text in the
                       the file that explains what 'work' means.
f. contrib.html      - This file lists a small number of people that have
                       contributed ideas or code to QTFM.
g. license.html      - This file is linked from the foot of almost all of
                       the HTML files and it describes the license terms in
                       some detail, with links to the GNU GPL and FDL.

The master XML files for the various QTFM functions are the main part of
the help documentation. The XML files are stored in a sub-folder of the
'helpfiles' folder called 'xmlfiles'. Each XML file must be processed into
its corresponding HTML file. This is done sporadically and before each new
release to ensure that the HTML files in a release are up-to-date. The XML
files are distributed with the toolbox so that anyone who wishes to adapt
or extend the toolbox, as is permitted under the GNU GPL licence may do so.

The conversion from XML to HTML is done by Matlab (what else!). There is a
Matlab script file called 'process.m' in the 'helpfiles' folder that
performs this processing, which takes only a few seconds (it processes all
the files every time, from a list contained in the file 'process.m'. An
easier way to run this script is to run the function qtfm_helpup which is
in the main QTFM folder, and therefore on the search path.

The conversion is controlled by an XSLT style file invoked by the process
script. This style file is the most complex part of the whole system. The
file name of the style file is 'qtfmfunction.xsl'. The style sheet uses a
DTD file (document type definition) stored in the file 'qtfmfunction.dtd'.
This DTD file sets out the definitions of the XML markups used, which are
largely self-explanatory when studied alongside an example XML file and the
corresponding HTML file when displayed in the browser. Not all elements are
compulsory: the 'See also' section, for example, may be omitted.

3. Building a searchable database from the help information.

Matlab provides a way to make a searchable database of the help 
information. This can be done easily by running the function
qtfm_helpdb, which is located in the main QTFM folder, and is therefore on
the search path. The database is not distributed with the toolbox, because
each release of Matlab may create a new version.

4. PDF documentation

It is intended to use the same XML files as master documents for production
of LaTeX files which can be processed into PDF documentation for printing,
but this has not been done, although a tentative XML to LaTeX XSLT style
sheet has been written and distributed with the toolbox for some time.

The small number of HTML files which are not based on XML masters would
have to be edited manually into the LaTeX source or otherwise converted to
some source form from which the HTML and LaTeX could be produced. There are
no imminent plans to do this.


Steve Sangwine
May/June 2008
Revised March 2013
Rewritten June 2016

$Id: readme.txt 1004 2017-11-15 17:14:09Z sangwine $