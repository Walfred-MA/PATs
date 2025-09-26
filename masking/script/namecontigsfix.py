#!/usr/bin/env python3

import os
import argparse
import re

def main(args):


        write = open(args.output, mode ='w') 

        with open(args.input, mode = 'r') as f:

                for line in f:

                        if len(line) == 0 or line[0] != ">":
                                write.write(line)
                                continue
                        line = line.replace(args.name,"bugfix") 
                        line = line.split()
                        line[0] = ">"+"#".join([x for x in line[0][1:].split("#") if x != "bugfix"])
                        line[0] = ">"+ args.name +"#"+ line[0][1:]
                        line = "\t".join(line) + "\n"
                        write.write(line)

        write.close()

def run():
        """
                Parse arguments and run
        """
        parser = argparse.ArgumentParser(description="program edit assembly contig names for PATs")
        parser.add_argument("-i", "--input", help="path to input data file",dest="input", type=str, required=True)
        parser.add_argument("-n", "--name", help="path to input data file",dest="name", type=str, required=True)
        parser.add_argument("-o", "--output", help="path to input data file",dest="output", type=str, required=True)

        parser.set_defaults(func=main)
        args = parser.parse_args()
        args.func(args)

if __name__ == "__main__":
        run()
