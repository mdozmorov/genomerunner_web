import collections
import gzip
import io
import itertools
import re

import pandas
import numpy

ExpressionSet = collections.namedtuple("ExpressionSet",
                                       "phenotype_data,feature_data,expression")

def as_float(item):
    try:
        return float(item)
    except:
        return numpy.nan
        
def read_table(lines, end_tag):
    buffer = io.StringIO()
    for line in lines:
        if line.startswith(end_tag):
            buffer.seek(0)
            return pandas.io.parsers.read_csv(buffer, sep="\t")
        buffer.write(line)

class Sample(object):
    """
    Represents a GEO GSM.
    """
    def __init__(self, id, expression, characteristics=[]):
        self.id = id
        self.expression = expression
        self.characteristics = characteristics

    def __repr__(self):
        return "<Sample %s with %d probes>" % (self.id, len(self.expression))
        
class SOFTParser(object):
    """
    A class to read from whole-platform SOFT datasets, such
    as those found at ftp://ftp.ncbi.nlm.nih.gov/geo/platforms/ . 

    It assumes the platform data will precede the samples.
    """
    def __init__(self, path):
        self._handle = gzip.open(path)
        self._lines = self._read_lines()
        self._read_until("!platform_table_begin")
        self.feature_data = self._read_table("!platform_table_end")

        expression = {}
        properties = {}

        for sample in self:
            print(sample.id)
            # Try to parse basic sample attributes, like age, gender, etc.
            p = {}
            for ch in sample.characteristics:
                age_match = re.match("^age: ([0-9]+\.*[0-9]*)", ch)
                if age_match:
                    p["age"] = float(age_match.group(1))
                if ch.startswith("gender:") or ch.startswith("sex:"):
                    if "M" in ch or "male" in ch.lower():
                        p["gender"] = "M"
                    elif "F" in ch or "female" in ch.lower():
                        p["gender"] = "F"
            if p:
                expression[sample.id] = sample.expression
                properties[sample.id] = pandas.Series(p)

        self.expression = pandas.DataFrame(expression).T
        self.phenotype_data = pandas.DataFrame(properties).T
        self.phenotype_data["age"] = self.phenotype_data["age"].astype("float64")

    def _read_lines(self):
        for line in self._handle:
            try:
                yield line.decode("UTF-8")
            except UnicodeDecodeError:
                continue
                
    def _read_table(self, end_tag):
        """Read a platform or series table into a DataFrame."""
        buffer = io.StringIO()
        for line in self._lines:
            if line.startswith(end_tag):
                buffer.seek(0)
                return pandas.io.parsers.read_csv(buffer, sep="\t")
            buffer.write(line)

    def _read_until(self, tag):
        """Advance the internal lines iterator to the given tag."""
        for line in self._lines:
            if line.startswith(tag):
                return

    def _read_value(self, line):
        return "".join(line.strip().split(" = ")[1:])

    def __iter__(self):
        """Iterate through the Sample objects in this SOFT file."""
        sample_id = None
        characteristics = []
        for line in self._lines:
            if line.startswith("^SAMPLE"):
                sample_id = line.strip().split()[2]
            elif line.startswith("!Sample_characteristics_ch1"):
                characteristics.append(self._read_value(line))
            elif line.startswith("!sample_table_begin"):
                try:
                    expression = self._read_table("!sample_table_end")
                    #expression_vector = expression["VALUE"]
                    #expression_vector.index = expression["ID_REF"]
                    expression_vector = pandas.Series([as_float(item) for item in expression["VALUE"]],
                                                      index=expression["ID_REF"])
                    yield Sample(sample_id, expression_vector, characteristics=characteristics)
                    sample_id = None
                    characteristics = []
                except Exception as e:
                    print(e)
                    continue

def parseSOFT(path):
   parser = SOFTParser(path) 
   return ExpressionSet(parser.phenotype_data,
                        parser.feature_data, parser.expression)
