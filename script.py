#!/usr/bin/python
from __init__ import *
import sys

def main(fname, outdir=""):
  run_stepminer(fname, outdir)

if __name__ == "__main__":
  args = dict([s.split('=') for s in sys.argv[1:]])
  main(**args)
