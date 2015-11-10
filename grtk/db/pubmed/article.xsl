<xsl:stylesheet version="1.0"
    xmlns:xsl="http://www.w3.org/1999/XSL/Transform">
<xsl:output method="text" />

<xsl:template match="MedlineCitation">
    <xsl:value-of select="PMID" />
    <xsl:text>&#x9;</xsl:text>
    <xsl:value-of select="translate(Article/ArticleTitle, '&#x9;&#xa;', ' ')" />
    <xsl:text>&#x9;</xsl:text>
    <xsl:value-of select="translate(Article/Abstract/AbstractText, '&#x9;&#xa;', ' ')" />
</xsl:template>

</xsl:stylesheet>
