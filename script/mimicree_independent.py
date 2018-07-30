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

def process_cmd(mkhapcmdp,mimcrcmd,adrltcmd,num,rflag=False,dist="", base = 0, random_dom=False):
    cp = cproc(False,False)
    if dist != "":
        dist_df = np.loadtxt(dist)
    for i in range(base, base+num):
        #Make Haplotype
        mkhapcmd = " ".join([mkhapcmdp,"-n",str(i)])
        if rflag:
            slc = random.uniform(0,0.3)
            mkhapcmd = " ".join([mkhapcmd,"-s",str(slc)])
            if random_dom:
                dom = random.uniform(0,1)
                mkhapcmd = " ".join([mkhapcmd,"-d",str(dom)])
        elif dist != "":
            slc = dist_df[sampleCategory(dist_df[:,1]),0]
            dom = random.uniform(0,1)
            mkhapcmd = " ".join([mkhapcmd,"-s",str(slc),"-d",str(dom)])
        cp.add(mkhapcmd)
        #Run simulation 
        cp.add(mimcrcmd)
        #Add result to all results file
        cp.add(adrltcmd)
    cp.exe()

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='make sync into dat')
    parser.add_argument('--tflag', '-t',default=False, action = "store_true",help='test flag not qsub')
    parser.add_argument('--num', '-n',default=2,type = int,help='iteration number')
    parser.add_argument('--tsim', '-i',default="tmp/tmp_sim.txt",type = str,help='file name of tmp simulation')
    parser.add_argument('--hap', '-l',default="tmp/hap_test_one.mimhap",type = str,help='file name of haplotype file')
    parser.add_argument('--recomb', '-c',default="tmp/data/dmel.rr.txt",type = str,help='file name of recombination rate file')
    parser.add_argument('--gen', '-g',default="0,10",type = str,help='list of generation')
    parser.add_argument('--rep', '-r',default="3",type = str,help='number of replicate')
    parser.add_argument('--pop', '-p',default="500",type = str,help='number of population')
    parser.add_argument('--bound', '-b',default="0.1,0.9",type = str,help='lower bound of initial allele frequency') # 
    parser.add_argument('--snp', '-v',default="tmp/test/snp_test.snp",type = str,help='list of generation')
    parser.add_argument('--slc', '-s',default="0.1",type = str,help='selection')
    parser.add_argument('--dom', '-d',default="0.5",type = str,help='dominance')
    parser.add_argument('--afs', '-a',default="",type = str,help='file of AF spectrum')
    parser.add_argument('--output', '-o',type=str,help='haplotype file name')
    parser.add_argument('--mflag', '-m',default=False, action = "store_true",help='make minus selection')
    parser.add_argument('--rflag', '-u',default=False, action = "store_true",help='make selection random')
    parser.add_argument('--randomdom',default=False, action = "store_true",help='make dominance random')
    parser.add_argument('--dist', '-e',default="",type = str,help='list of generation')
    args = parser.parse_args()
    cp = cproc(False,args.tflag)
    #make empty file
    if os.path.exists(args.output)==True:
        cp.add("rm " + args.output)
        cp.exe()
    if os.path.exists(args.snp)==True:
        cp.add("rm " + args.snp)
        cp.exe()
    #make command of haplotype
    mkhapcmd = "python script/make_hap.py"
    mkhapcmd = " ".join([mkhapcmd,"-p",args.pop])
    mkhapcmd = " ".join([mkhapcmd,"-b",args.bound])
    mkhapcmd = " ".join([mkhapcmd,"-o",args.hap])
    mkhapcmd = " ".join([mkhapcmd,"-v",args.snp + "_tmp"])
    mkhapcmd = " ".join([mkhapcmd,"-a",args.afs])
    mkhapcmd = " ".join([mkhapcmd,"-d",args.dom])
    mkhapcmdp = " ".join([mkhapcmd,"-s",args.slc])
    mkhapcmdn = " ".join([mkhapcmd,"-s","-"+args.slc])
    #make command of mimicree
    mimcrcmd = "java -Xmx1g -jar  script/MimicrEESummary.jar --output-format sync --threads 1"
    mimcrcmd = " ".join([mimcrcmd,"--recombination-rate",args.recomb])
    mimcrcmd = " ".join([mimcrcmd,"--haplotypes-g0",args.hap])    
    mimcrcmd = " ".join([mimcrcmd,"--output-mode",args.gen])    
    mimcrcmd = " ".join([mimcrcmd,"--replicate-runs",args.rep])    
    mimcrcmd = " ".join([mimcrcmd,"--additive",args.snp + "_tmp"])
    mimcrcmd = " ".join([mimcrcmd,"--output-file",args.tsim])
    #make command of summarization of results
    adrltcmd = "cat"
    adrltcmd = " ".join([adrltcmd,args.tsim,'>>',args.output])
    #make command of summarization of results
    adsnpcmd = "cat"
    adsnpcmd = " ".join([adsnpcmd,args.snp + "_tmp",'>>',args.snp])
    adrltcmd = adrltcmd + ";" + adsnpcmd
    #selection is minus
    #make command of haplotype negative selection
    if args.rflag or (args.dist != ""):
        process_cmd(mkhapcmd,mimcrcmd,adrltcmd,args.num,args.rflag,args.dist, args.randomdom)
    elif args.mflag:
        process_cmd(mkhapcmdp,mimcrcmd,adrltcmd,args.num/2)
        process_cmd(mkhapcmdn,mimcrcmd,adrltcmd,args.num - args.num/2, base=args.num/2)
    else:
        process_cmd(mkhapcmdp,mimcrcmd,adrltcmd,args.num)
