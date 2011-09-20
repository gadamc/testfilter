#!/usr/bin/env python

from couchdbkit import Server, Database
import os, sys, copy, shutil
import numpy as np
from ROOT import *
import matplotlib.pyplot as plt


def main(*args):
  
  f = KDataReader(args[0])
  


if __name__ == '__main__':
  main(*sys.argv[1:])