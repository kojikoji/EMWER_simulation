#!/usr/env python
#-*- coding: utf-8 -*-
import random
import re
import argparse
import os
import numpy as np
from cproc.cproc import *
import numpy as np

def sampleCategory(p):
    return np.flatnonzero( np.random.multinomial(1,p,1) )[0]


def add_opt_val(cmd, key, val):
    return(' '.join([cmd, key, str(val)]))


def add_opt(cmd, opt):
    return(' '.join([cmd, opt]))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='make sync into dat')
    parser.add_argument('--tflag', '-t',default=False, action = "store_true",help='test flag not qsub')
    parser.add_argument('--num', '-n',default=2,type = int,help='iteration number')
    parser.add_argument('--tsim', '-i',default="tmp/tmp_sim.txt",type = str,help='file name of tmp simulation')
    parser.add_argument('--hap', '-l',default="tmp/hap_test_one.mimhap",type = str,help='file name of tmp haplotype file')
    parser.add_argument('--recomb', '-c',default="2.0",type = str,help='Recombination rate')
    parser.add_argument('--gen', '-g',default="0,10",type = str,help='list of generation')
    parser.add_argument('--rep', '-r',default="3",type = str,help='number of replicate')
    parser.add_argument('--pop', '-p',default="500",type = str,help='number of population')
    parser.add_argument('--bound', '-b',default="0.05,0.95",type = str,help='lower bound of initial allele frequency') # 
    parser.add_argument('--slc', '-s',default="0.1",type = str,help='selection')
    parser.add_argument('--dom', '-d',default="0.5",type = str,help='dominance')
    parser.add_argument('--afs', '-a',default="",type = str,help='file of AF spectrum')
    parser.add_argument('--snp', '-v',default="tmp/test/snp_test.snp",type = str,help='list of generation')
    parser.add_argument('--output', '-o',type=str,help='haplotype file name')
    parser.add_argument('--dist', '-e',default="",type = str,help='list of generation')
    parser.add_argument('--depth',default="50",type = str,help='mean sequencing depth')
    parser.add_argument('--mflag', '-m',default=False, action = "store_true",help='make minus selection')
    parser.add_argument('--rflag', '-u',default=False, action = "store_true",help='make selection random')
    parser.add_argument('--cont', default=False, action = "store_true",help='considering linkage')
    parser.add_argument('--randomdom',default=False, action = "store_true",help='make dominance random')
    parser.add_argument('--snum', '-j',default="1",type = str,help='number of slected SNPs')
    args = parser.parse_args()
    cp = cproc(False,args.tflag)
    # command for recomination file
    cmdrecomb = "cat constant_recombination.txt | awk '{print $1\"\t\"" + args.recomb + "\"\t\"" + args.recomb + "\"\t\"" + args.recomb + "}' > tmp/recomb.txt"
    cp.add(cmdrecomb)
    # command for mimicree simulation
    if args.cont:
        cmdsim = "python script/mimicree_run.py"
    else:
        cmdsim = "python script/mimicree_independent.py"
    cmdsim = cmdsim +  " -i tmp/sim.osync -l tmp/haplotype.mimhap -c tmp/recomb.txt -a afs/afs_unif.afs"
    cmdsim = add_opt_val(cmdsim, "-n", args.num)
    cmdsim = add_opt_val(cmdsim, "-p", args.pop)
    cmdsim = add_opt_val(cmdsim, "-g", args.gen)
    cmdsim = add_opt_val(cmdsim, "-r", args.rep)
    cmdsim = add_opt_val(cmdsim, "-b", args.bound)
    cmdsim = add_opt_val(cmdsim, "-s", args.slc)
    cmdsim = add_opt_val(cmdsim, "-d", args.dom)
    cmdsim = add_opt_val(cmdsim, "-v", args.snp)
    cmdsim = add_opt_val(cmdsim, "-o", "tmp/presim.sync")
    if args.rflag:
        cmdsim = add_opt(cmdsim, "-u")
    if args.mflag:
        cmdsim = add_opt(cmdsim, "-m")
    if args.randomdom:
        cmdsim = add_opt(cmdsim, "--randomdom")
    if args.cont:
        cmdsim = add_opt_val(cmdsim, "-j", args.snum)
    cp.add(cmdsim)
    # command for poisson sampling
    cmdpoisson = "python script/poisson-3fold-sample.py --input tmp/presim.sync"
    cmdpoisson = add_opt_val(cmdpoisson, "--coverage", args.depth)
    cmdpoisson = cmdpoisson + " > "  + args.output
    cp.add(cmdpoisson)
    cp.exe()
