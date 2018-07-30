#!/usr/env python
#-*- coding: utf-8 -*-
import random
import re
import numpy.random as rd
import argparse
import numpy as np
import sys
import matplotlib.pyplot as plt
def sum_regulalize(numpy_array,regulalized_sumed_value):
    sumed_value = sum(numpy_array)
    print(sumed_value)
    regulalized_numpy_array = (regulalized_sumed_value/sumed_value)*numpy_array
    return(regulalized_numpy_array)
def cumrative(numarray):
    cumrative_array = np.zeros(len(numarray))
    cumrative_array[0] = numarray[0]
    for i in range(1,len(numarray)):
        cumrative_array[i] = cumrative_array[i-1] + numarray[i]
    return(cumrative_array)
def cumrative_sample(cumrative_array,pointer):
    index = 0
    for cumrative_value in cumrative_array:
        if cumrative_value > pointer:
            break;
        index += 1
    return(index)
class Bin:
    def __init__(self,low,high,scale):
        self.low = float(low)
        self.high = float(high)
        self.scale = float(scale)
    def get(self,value):
        scaled_value = ((self.high - self.low)/self.scale)*value
        if(scaled_value < 0 or scaled_value > (self.high - self.low)):
            print >> sys.stderr,'Error:This value is bigger than upper bound'
        value_in_bin = self.low + scaled_value
        return value_in_bin
class Spectrum_uniform:
    def set_low_high(self,low,high):
        self.low = low
        self.high = high
    def set_spectrum(self,spectrum):
        lower_bound = self.low
        self.bin_width = (self.high - self.low)/float(len(spectrum))
        self.cumrative_spectrum = cumrative(spectrum)
        self.cumrative_spectrum_anterior = self.cumrative_spectrum - spectrum
        print(self.cumrative_spectrum)
        self.bin_list = list()
        for measure in spectrum:
            upper_bound = lower_bound + self.bin_width
            self.bin_list.append(Bin(lower_bound,upper_bound,measure))
            lower_bound = upper_bound
    def sample(self,low,high):
        pointer = random.uniform(0,self.cumrative_spectrum[-1])
        sample_index = cumrative_sample(self.cumrative_spectrum,pointer)
        inner_pointer = pointer - self.cumrative_spectrum_anterior[sample_index]
        sample_value = self.bin_list[sample_index].get(inner_pointer)
        if(sample_value < low or sample_value > high):
            sample_value = self.sample(low,high)
        return(sample_value)

#from util_wfe import *
tag = "2L\t300001\tC\t"
def make_hap(al1,al2,pop,lb,hb,afs_file,num, chrom="2L"):
    su = Spectrum_uniform()
    su.set_low_high(lb,hb)
    if afs_file == "":
        afs_vec = np.ones(50)
    else:
        afs_vec = np.loadtxt(afs_file)
    su.set_spectrum(afs_vec)
    alprob = su.sample(lb,hb)
    alnum = int(alprob*(pop*2))
    als = al1*alnum + al2*(pop*2-alnum)
    als = ''.join(random.sample(als,len(als)))
    als = re.sub(r'(..)',r'\1 ',als)
    ref = "C"
    als = "\t".join([chrom,str(num), ref, "C/A",als]) +"\n"
    return(als)
def make_snp(selection,dom,num, chrom="2L"):
    ref = "C"
    return("\t".join([chrom,str(num), ref,str(selection),str(dom)]) + '\n')
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='make sync into dat')
    parser.add_argument('--population', '-p',type = str,help='population')
    parser.add_argument('--selection', '-s',default = '0.1',type = str,help='selection')
    parser.add_argument('--dominance', '-d',default = '0.5',type = str,help='dominance')
    parser.add_argument('--bound', '-b',type = str,help='boundary of allele frequency: lower,upper')
    parser.add_argument('--output', '-o',type=str,help='haplotype file name')
    parser.add_argument('--snp', '-v',default="log.txt",type = str,help='file name of selected snp')
    parser.add_argument('--afs', '-a',default="",type = str,help='file name of allele frequency spectrum')
    parser.add_argument('--num', '-n',default="3001",type = str,help='number of SNPs')
    parser.add_argument('--snum', '-u',default="-1",type = str,help='number of slected SNPs')
    parser.add_argument('--cont', '-c',default=False, action = "store_true",help='continuous snp file')
    args = parser.parse_args()
            # loci array
    upper_loci_dict = {"2L":23011543, "2R":21146708, "3L":24543557, "3R":27905053}
    num_vec = rd.multinomial(int(args.num),np.array([0.25,0.25,0.25,0.25]))
    num_dict = {"2L": num_vec[0], "2R": num_vec[1],
                "3L": num_vec[2], "3R": num_vec[3]}
    chrom_list = ["2L", "2R", "3L", "3R"]
    loc_vec_list = [np.sort(
        np.random.choice(
            upper_loci_dict[chrom], size=num_dict[chrom], replace=False)) \
                    for chrom in chrom_list]
    loc_vec_all = np.concatenate(loc_vec_list)
    chrom_vec_all = np.concatenate(
        [np.array([chrom for _ in range(num_dict[chrom])])
         for chrom in  chrom_list])
    with open(args.output,'w') as hap_file, open(args.snp,'w') as snp_file:
        lower = float((args.bound).split(',')[0])
        upper = float((args.bound).split(',')[1])
        if args.cont:
            # index of selected index
            snum = int(args.snum) if int(args.snum) > 0 else int(args.num)
            s_index_vec = np.random.choice(int(args.num), size=int(snum), replace=False)
            for i in range(int(args.num)):
                hap_file.write(make_hap(
                    'A','C',int(args.population),
                    lower,upper, args.afs,
                    loc_vec_all[i], chrom=chrom_vec_all[i]))
                slc_val = "0"
                if i in s_index_vec:
                    slc_val = args.selection
                snp_file.write(make_snp(
                    slc_val, args.dominance,
                    loc_vec_all[i], chrom=chrom_vec_all[i]))
        else:
            hap_file.write(make_hap('A','C',int(args.population),lower,upper,args.afs,args.num))
            snp_file.write(make_snp(args.selection,args.dominance,args.num))
