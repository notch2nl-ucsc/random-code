import os
import sys
import pysam
from collections import defaultdict, namedtuple
from itertools import tee, izip, groupby

def fastq_iter(fq):
  n = 4
  def grouped(iterator):
    while True:
      vals = tuple(next(iterator, None) for _ in xrange(n))
      if None not in vals:
        yield vals
      else:
        raise StopIteration

  assert os.path.isfile(fq)
  with open(fq) as f:
    for lines in grouped(f):
      txt = ''.join(lines)
      qname, bcode = lines[0].strip().split()
      yield (bcode, qname, txt)

  raise StopIteration

def sortfq(
  inputFq_path,
  outputFq_path,
):

  print 'loading'
  bcode_entries_map = defaultdict(list)
  for qname, reads_iter in groupby(
    fastq_iter(inputFq_path),
    lambda(x): x[1],
  ):
    for bcode, qname, txt in reads_iter:
      bcode_entries_map[bcode].append(txt)

  print 'writing'
  with open(outputFq_path, 'w') as fout:
    for bcode, entries in bcode_entries_map.items():
      for entry in entries:
        fout.write(entry)

  print 'done!'

def main(argv):

  # argument parsing
  if len(argv) < 3:
    print 'sortfq.py <in fq>  <out fq>'
    sys.exit(1)
  inputFq_path = argv[1]
  outputFq_path = argv[2]

  if False in [
    os.path.isfile(inputFq_path),
  ]:
    print 'error: input fqs not files'
    sys.exit(1)

  sortfq(
    inputFq_path,
    outputFq_path
  ) 

if __name__ == '__main__':
  main(sys.argv)

