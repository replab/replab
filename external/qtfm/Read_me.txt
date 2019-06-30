Quaternion toolbox for Matlab - installation instructions
---------------------------------------------------------

To install this toolbox:

1. Unzip the distribution file and move the directory/folder qtfm
   to a convenient location. (On Windows some unzip utilities may
   put this folder inside another of the same name - you can use
   just the innermost qtfm folder.)
   The location does not need to be in the same location as toolboxes
   supplied with Matlab, although it could be if your Matlab
   installation is to be used by multiple users on the same machine.

2. Set the Matlab path to include the directory/folder qtfm.
   This is done from the Matlab File -> Set Path menu up to
   version 7 of Matlab, or the Set Path button on the toolbar in
   later versions.
   [The QTFM folder must be near the top of the path, i.e.
   higher than the standard Matlab folders, otherwise the
   overloading of Matlab functions will not work.]

3. Help information is available under Supplemental Software in the
   Help window. (Older versions of Matlab displayed the help for this
   toolbox in different ways.) You can also access help information in
   the command window: try help qtfm. The command qtfm_helpdb will build
   a searchable database from the help files so that searches in the
   help window can return results from the QTFM documentation.

4. help quaternion/xxx shows help text for each function xxx
   and more detailed help is available in the help documentation
   (see 3 above).

5. To run a test of the toolbox, type the command 'qtfm_test'
   in the Matlab command window. This runs a test of many parts
   of the toolbox and will allow you to confirm that installation
   is correct. It is also possible to invoke the test code from
   the Matlab start menu prior to Matlab 8
   (Toolboxes -> Quaternion -> Run test code).

6. If you find the toolbox useful and would like to be kept
   informed of updates and future releases, subscribe to the
   mailing list qtfm-announce@lists.sourceforge.net. You can
   subscribe to this list at:
   https://lists.sourceforge.net/lists/listinfo/qtfm-announce

7. LaTeX/BiBTeX users please note the file qtfm.bib provided
   in the directory tex/bibtex if you wish to cite QTFM itself
   or one of the published papers cited within the source code.
   There is a separate entry in the file for QFTM version 2.
   Please use this if you cite QTFM for work on octonions.

8. We would welcome contributed code, ideas for new functions,
   or expanded examples to include in the help documentation.
   (Contributed code must be made available under the same
   license terms as QTFM itself). Please contact us by email
   as in below.

9. To send feedback, code, documentation corrections or bug
   reports, please use the reporting mechanisms at the project
   webpage at Sourceforge

   https://sourceforge.net/projects/qtfm/

   or send email to:

   sangwine@users.sourceforge.net
   n-le_bihan@users.sourceforge.net

Steve Sangwine
Nicolas Le Bihan

27 March 2006
Updated 5 June 2006 and 23 May 2008 and 2 April 2013 and 27 January 2016.
$Id: Read_me.txt 1004 2017-11-15 17:14:09Z sangwine $
