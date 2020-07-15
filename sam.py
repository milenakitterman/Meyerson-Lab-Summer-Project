import re

class SamLine:
  def __init__(self, line):
    self.parse(line)

  def parse(self, line):
    f = line.split()
    self.qname, self.flag, self.rname, self.pos, self.mapq, self.cigar = f[:6]
    self.seq, self.qual = f[9:11]
    self.flag = int(self.flag)
    self.pos = int(self.pos)
    self.mapq = int(self.mapq)

  def isPaired(self):
    return self.flag & 1

  def isUnMapped(self):
    return self.flag & 4

  def isReverseComplemented(self):
    return self.flag & 16
  
  def applyCIGAR(self):
    pattern = "[0-9]+[MIDNSHPX=]"
    ops = [(int(r.group(0)[:-1]), r.group(0)[-1]) \
        for r in re.finditer(pattern, self.cigar)]
    idx = 0
    output_seq = ""
    for num, op in ops:
      if op == "I":
        idx += num
      elif op == "D":
        output_seq += '.' * num
      elif op == "N":
        output_seq += '.' * num
      elif op == "S":
        output_seq += self.seq[idx:idx+num].lower()
        idx += num
      elif op == "P":
        pass
      elif op == "H":
        pass
      else:
        output_seq += self.seq[idx:idx+num]
        idx += num

    return output_seq
  
def findReverseComplement(seq):
  comp= {'A':'T', 'G':'C', 'T':'A', 'C':'G', 'a':'t', 'g':'c', 't':'a', 'c':'g', 'N':'N', '.':'.'}
  dna = [comp[char] for char in seq[::-1]]
  return "".join(dna)



