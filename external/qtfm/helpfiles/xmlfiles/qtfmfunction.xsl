<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE function SYSTEM "qtfmfunction.dtd">
<!-- This XML stylesheet converts an XML description of a QTFM function
     into an HTML documentation page for display in the Matlab help
     browser. See the file readme.txt for further details.

     (c) Stephen J. Sangwine and Nicolas Le Bihan, 2008-2019. -->
     
<xsl:stylesheet version="1.0"
                xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
                
<xsl:output method="html" 
            version="4.01"
            encoding="UTF-8" 
            omit-xml-declaration="yes" 
            indent="no"/>
                
<xsl:template match="function">
<xsl:comment> Generated automatically from an XML file of the same name.
     Copyright: Stephen J. Sangwine and Nicolas Le Bihan, 2008-2019.
</xsl:comment>
<html>
<head>
<title>
<xsl:value-of select="@name"/> :: Functions (Quaternion Toolbox Function Reference)
</title>
<link rel="stylesheet" href="qtfmstyle.css" type="text/css"></link>
</head>
<body>
<h1>Quaternion Function Reference</h1>
<h2><xsl:value-of select="@name"/></h2>
<xsl:apply-templates/>
<h4>(c) 2008-2019 Stephen J. Sangwine and Nicolas Le Bihan</h4>
<p><a href="license.html">License terms.</a></p>
</body>
</html>
</xsl:template>

<xsl:template match="short">
<p>
<xsl:value-of select="text()"/>
<xsl:if test = "/function/@overload = 'true' ">
<br/>(Quaternion overloading of standard &matlab; function)
</xsl:if>
<xsl:if test = "/function/@overload = 'both' ">
<br/>(Quaternion and octonion overloadings of standard &matlab; function)
</xsl:if>
</p>
</xsl:template>

<xsl:template match="syntax">
<h2>Syntax</h2>
<p><tt><xsl:value-of select="text()"/></tt></p>
</xsl:template>

<xsl:template match="long">
<h2>Description</h2>
<xsl:apply-templates/>
</xsl:template>

<xsl:template match="para">
<p><xsl:apply-templates/></p>
</xsl:template>

<xsl:template match="code">
<tt><xsl:value-of select="text()"/></tt>
</xsl:template>

<xsl:template match="precode">
<pre><xsl:value-of select="text()"/></pre>
</xsl:template>

<xsl:template match="figure">
<xsl:element name="img">
<xsl:attribute name="src"><xsl:value-of select="text()"/></xsl:attribute>
</xsl:element>
</xsl:template>

<xsl:template match="examples">
<h2>Examples</h2>
<xsl:apply-templates/>
</xsl:template>

<xsl:template match="seealso">
<h2>See Also</h2>
<!-- Here we output a list of comma separated hyperlinks, one for QTFM functions and one for Matlab
     functions. What we output differs according to whether there are 0, 1 or more functions to be
     listed. -->
<xsl:if test="count(qtfmfun)=1">QTFM function: <xsl:apply-templates select="qtfmfun"/><br/></xsl:if>
<xsl:if test="count(matlabfun)=1">&matlab; function: <xsl:apply-templates select="matlabfun"/><br/></xsl:if>
<xsl:if test="count(qtfmfun)>1">QTFM functions: <xsl:apply-templates select="qtfmfun"/><br/></xsl:if>
<xsl:if test="count(matlabfun)>1">&matlab; functions: <xsl:apply-templates select="matlabfun"/><br/></xsl:if>
</xsl:template>

<xsl:template match="matlabfun">
<xsl:element name="a">
<xsl:attribute name="href">matlab:doc <xsl:value-of select="@name"/></xsl:attribute>
<xsl:value-of select="@name"/>
</xsl:element>
<xsl:if test="not(position()=last())">, </xsl:if>
</xsl:template>

<xsl:template match="qtfmfun">
<xsl:element name="a">
<!-- Here we output different links depending on whether the alias attribute
     was present or not. If absent, the HTML file which is referenced has
     the same name as the XML file. If present, the HTML file has the alias
     name rather than the function name. -->
<xsl:if test="@alias">
  <xsl:attribute name="href"><xsl:value-of select="@alias"/>.html</xsl:attribute>
</xsl:if>
<xsl:if test="not(@alias)">
  <xsl:attribute name="href"><xsl:value-of select="@name"/>.html</xsl:attribute>
</xsl:if>
<xsl:value-of select="@name"/>
</xsl:element>
<xsl:if test="not(position()=last())">, </xsl:if>
</xsl:template>

<xsl:template match="references"><h2>References</h2><ol><xsl:apply-templates select="reference"/></ol></xsl:template>

<xsl:template match="reference"><li><xsl:apply-templates/></li></xsl:template>

<xsl:template match="doi">
DOI: <xsl:element name="a">
<xsl:attribute name="href">http://dx.doi.org/<xsl:value-of select="text()"/></xsl:attribute>
<xsl:value-of select="text()"/>
</xsl:element>
</xsl:template>

<xsl:template match="www">
<xsl:element name="a"><xsl:attribute name="href"><xsl:value-of select="text()"/></xsl:attribute><xsl:value-of select="text()"/></xsl:element>
</xsl:template>

<xsl:template match="italic">
<i><xsl:apply-templates/></i>
</xsl:template>

<xsl:template match="bold">
<b><xsl:apply-templates/></b>
</xsl:template>

</xsl:stylesheet>
<!--$Id: qtfmfunction.xsl 1018 2019-02-21 12:01:28Z sangwine $-->
