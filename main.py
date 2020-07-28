#!/usr/bin/env python3
import sys
import re
from sam import SamLine

def main():
  if len(sys.argv) < 2:
    sys.stdout.write("Usage: %s <input.sam>\n"%sys.argv[0])
    # exit the program if no argument was provided. this is a convention for
    # programs, because users usually want to try to run the program without
    # any arguments and expect to see a welcome/help message.
    sys.exit()

  with open(sys.argv[1], 'r') as inf:
    sam_data = []
    for l in inf:
      if not l.startswith('@'):
        # convert a data line into a SamLine object and store it in the list
        sam_data.append(SamLine(l))

  # more stats can be added here
  paired = 0
  for sam_read in sam_data:
    pattern = "[0-9]+[MIDNSHPX=]"
    ops = [(int(r.group(0)[:-1]), r.group(0)[-1]) \
      for r in re.finditer(pattern, sam_read.cigar)]
    a = sam_read.pos
    b = 0
    for num, op in ops:
      if op == "S":
        b = num
    x = sam_read.applyCIGAR()
    if x[1].islower():
      a = a-b
    sys.stdout.write("".join((sam_read.qname, " " * a, x))+"\n")
    # the following line is not working because the applyCIGAR function in
    # incomplete.
    #sys.stdout.write("\t".join((sam_read.qname, sam_read.applyCIGAR))+"\n")
    if sam_read.isPaired:
      paired += 1

  sys.stdout.write("Total reads: %d; total paired: %d\n"\
    %(len(sam_data), paired))

# usually you only need these codes below if you are writing some executable,
# e.g. you expect to "run" this program. Otherwise, for example the sam.py
# which you only need through import, does not need them.
if __name__ == "__main__":
  main()
