#!/usr/env python
#-*- coding: utf-8 -*-
import random
import re
import argparse
import os
from cproc.cproc import *


def process_cmd(mkhapcmd, mimcrcmd, adrltcmd, num, rflag=False):
    cp = cproc(False, False)
    for i in range(num):
        if rflag:
            slc = random.uniform(-0.3, 0.3)
            mkhapcmd = " ".join([mkhapcmd, "-s", str(slc), "-n", str(i)])
        # Make Haplotype
        cp.add(mkhapcmd)
        # Run simulation
        cp.add(mimcrcmd)
        # Add result to all results file
        cp.add(adrltcmd)
    cp.exe()


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='make sync into dat')
    parser.add_argument('--tflag', '-t', default=False,
                        action="store_true", help='test flag not qsub')
    parser.add_argument('--num', '-n', default=2,
                        type=int, help='iteration number')
    parser.add_argument('--tsim', '-i', default="tmp/tmp_sim.txt",
                        type=str, help='file name of tmp simulation')
    parser.add_argument('--hap', '-l', default="tmp/hap_test_one.mimhap",
                        type=str, help='file name of haplotype file')
    parser.add_argument('--recomb', '-c', default="'constant_recombination.txt'",
                        type=str, help='file name of recombination rate file')
    parser.add_argument('--gen', '-g', default="0,10",
                        type=str, help='list of generation')
    parser.add_argument('--rep', '-r', default="3",
                        type=str, help='number of replicate')
    parser.add_argument('--pop', '-p', default="500",
                        type=str, help='number of population')
    parser.add_argument('--bound', '-b', default="0.1,0.9",
                        type=str, help='lower bound of initial allele frequency')
    parser.add_argument('--snp', '-v', default="tmp/test/snp_test.snp",
                        type=str, help='list of generation')
    parser.add_argument('--slc', '-s', default="0.1",
                        type=str, help='selection')
    parser.add_argument('--dom', '-d', default="0.5",
                        type=str, help='dominance')
    parser.add_argument('--afs', '-a', default="",
                        type=str, help='file of AF spectrum')
    parser.add_argument('--output', '-o', type=str, help='haplotype file name')
    parser.add_argument('--mflag', '-m', default=False,
                        action="store_true", help='make minus selection')
    parser.add_argument('--rflag', '-u', default=False,
                        action="store_true", help='make selection random')
    parser.add_argument('--snum', '-j',default="-1",type = str,help='number of slected SNPs')
    args = parser.parse_args()
    cp = cproc(False, args.tflag)
    # make empty file
    if os.path.exists(args.output) == True:
        cp.add("rm " + args.output)
        cp.exe()
    if os.path.exists(args.snp) == True:
        cp.add("rm " + args.snp)
        cp.exe()
    # make command of haplotype
    mkhapcmd = "python script/make_hap.py"
    mkhapcmd = " ".join([mkhapcmd, "-p", args.pop])
    mkhapcmd = " ".join([mkhapcmd, "-b", args.bound])
    mkhapcmd = " ".join([mkhapcmd, "-o", args.hap])
    mkhapcmd = " ".join([mkhapcmd, "-v", args.snp])
    mkhapcmd = " ".join([mkhapcmd, "-d", args.dom])
    mkhapcmd = " ".join([mkhapcmd, "-a", args.afs])
    mkhapcmd = " ".join([mkhapcmd, "-u", args.snum])
    mkhapcmd = " ".join([mkhapcmd, "-c"])
    mkhapcmd = " ".join([mkhapcmd, "-s", str(args.slc), "-n", str(args.num)])
    # make command of mimicree
    mimcrcmd = "java -Xmx1g -jar  script/MimicrEESummary.jar --output-format sync --threads 1"
    mimcrcmd = " ".join([mimcrcmd, "--recombination-rate", args.recomb])
    mimcrcmd = " ".join([mimcrcmd, "--haplotypes-g0", args.hap])
    mimcrcmd = " ".join([mimcrcmd, "--output-mode", args.gen])
    mimcrcmd = " ".join([mimcrcmd, "--replicate-runs", args.rep])
    mimcrcmd = " ".join([mimcrcmd, "--additive", args.snp])
    mimcrcmd = " ".join([mimcrcmd, "--output-file", args.output])
    cp.add(mkhapcmd)
    cp.add(mimcrcmd)
    cp.exe()
