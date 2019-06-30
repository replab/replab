<?xml version="1.0" encoding="UTF-8"?>
<!-- This XML stylesheet converts an XML description of a QTFM function
     into a LaTeX file. See the file readme.txt for further details.

     (c) Stephen J. Sangwine and Nicolas Le Bihan, 2008. -->
     
<xsl:stylesheet version="1.0"
                xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
                
<xsl:output method="text"
            encoding="UTF-8" 
            omit-xml-declaration="yes" 
            indent="no"/>
                
<xsl:template match="function">
% Generated automatically from an XML file of the same name.
% Copyright: Stephen J. Sangwine and Nicolas Le Bihan, 2008-2014.
\function{<xsl:value-of select="@name"/>}
<xsl:apply-templates/>
</xsl:template>

<xsl:template match="short">
\short{<xsl:value-of select="text()"/>}
<xsl:if test = "/function/@overload = 'true' ">
\overload
</xsl:if>
</xsl:template>

<xsl:template match="syntax">
\syntax{<xsl:value-of select="text()"/>}
</xsl:template>

<xsl:template match="long">
\description
</xsl:template>

<xsl:template match="para">
\par
<xsl:apply-templates/>
</xsl:template>

<xsl:template match="code">
\code{<xsl:apply-templates/>}
</xsl:template>

<xsl:template match="precode">
\precode{<xsl:value-of select="text()"/>}
</xsl:template>

<xsl:template match="examples">
\examples
<xsl:apply-templates/>
</xsl:template>

<xsl:template match="seealso">
\seealso
<!-- Here we output a list of comma separated hyperlinks, one for QTFM functions and one for Matlab
     functions. What we output differs according to whether there are 0, 1 or more functions to be
     listed.
  -->
<xsl:if test="count(qtfmfun)=1">QTFM function: <xsl:apply-templates select="qtfmfun"/><br/></xsl:if>
<xsl:if test="count(matlabfun)=1">MATLAB® function: <xsl:apply-templates select="matlabfun"/></xsl:if>
<xsl:if test="count(qtfmfun)>1">QTFM functions: <xsl:apply-templates select="qtfmfun"/><br/></xsl:if>
<xsl:if test="count(matlabfun)>1">MATLAB® functions: <xsl:apply-templates select="matlabfun"/></xsl:if>
</xsl:template>

<xsl:template match="matlabfun">
<xsl:value-of select="@name"/>
<xsl:if test="not(position()=last())">, </xsl:if>
</xsl:template>

<xsl:template match="qtfmfun">
<xsl:value-of select="@name"/>
<xsl:if test="not(position()=last())">, </xsl:if>
</xsl:template>

<xsl:template match="references">
\references
\begin{enumerate}
<xsl:apply-templates select="reference"/>
\end{enumerate}
</xsl:template>

<xsl:template match="reference">
\item <xsl:apply-templates/>
</xsl:template>

<xsl:template match="italic">
\textit{<xsl:apply-templates/>}
</xsl:template>

<xsl:template match="bold">
\textbf{<xsl:apply-templates/>}
</xsl:template>

</xsl:stylesheet>
<!--$Id: qtfmfunction2latex.xsl 1004 2017-11-15 17:14:09Z sangwine $-->
