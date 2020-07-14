"""This is the implementaton of a SAM format interpreter. Write all utilities
here, and you can use what you need in other programs via import.
"""
import re

class SamLine:
  """This class handles one data line of the SAM format.
  """
  def __init__(self, line):
    self.parse(line)

  def parse(self, line):
    """Parse a SAM data line, and fill attributes with parsed data.
    Attribute names correspond to the SAM specifications.
    """
    f = line.split()
    self.qname, self.flag, self.rname, self.pos, self.mapq, self.cigar = f[:6]
    self.seq, self.qual = f[9:11]
    self.flag = int(self.flag)
    self.pos = int(self.pos)
    self.mapq = int(self.mapq)

  def isPaired(self):
    """Returns true if the flag indicates this read is paired-end.
    """
    return self.flag & 1

  def isUnMapped(self):
    """Returns true if the flag indicates this read is unmapped.
    """
    return self.flag & 4

  def isReverseComplemented(self):
    """Returns true if the flag indicates this read is reverse complemented.
    """
    return self.flag & 16
  
  def applyCIGAR(self):
    """Apply the CIGAR string to align the read sequence with the reference.
    This function should only insert gap/delete nucleotides and mismatch is not
    considered.

    See the table of CIGAR operations in the SAM specification document for
    insertion/deletion related opeartions.

    Expected output: see the examples on the first page of SAM spec.
    """
    # break the CIGAR string to individual operations, if there are multiple.
    # regular expression is used
    pattern = "[0-9]+[MIDNSHPX=]"
    ops = [(int(r.group(0)[:-1]), r.group(0)[-1]) \
        for r in re.finditer(pattern, self.cigar)]
    # this part above should convert a cigar string to individual tuples,
    # each being (number, operation)
    idx = 0
    output_seq = ""
    for num, op in ops:
      if op == "I":
        idx += num
        # insertion. do nothing (for now)
      elif op == "D":
        # deletion. use "." to represent deleted nucleotides
        # change the pass below to something correct
        output_seq += '.' * num
      elif op == "N":
        # skipped. same as above
        # change the pass below to something correct
        output_seq += '.' * num
      elif op == "S":
        # soft clipping. change clipped nucleotides to lower case
        # change the pass below to something correct
        output_seq += self.seq[idx:idx+num].lower()
        idx += num
      elif op == "P":
        pass
      elif op == "H":
        pass
      else:
        # don't apply changes to mismatches for now
        output_seq += self.seq[idx:idx+num]
        idx += num

    return output_seq
  
def findReverseComplement(seq):
  dna = [char for char in seq]
  comp= {'A':'T', 'G':'C', 'T':'A', 'C':'G', 'a':'t', 'g':'c', 't':'a', 'c':'g', 'N':'N', '.':'.'}
  complist = []
  for val in dna:
    complist.append(comp[val])
  revcomplist = complist[::-1]
  revcomp = ""
  return(revcomp.join(revcomplist))


