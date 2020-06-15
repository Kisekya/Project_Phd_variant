# MULLER CORALIE
#Programming Assignment
# 15 -16 /06/2020

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
import argparse
import requests
import json

#######################
##    Functions      ##
#######################
#Function for check the format and print the help if its not a right format
def check_format(file):
    if file.lower().endswith(('.txt')):
        print("txt")
        #f = open('Variants', 'r')
    else:
        print("Sorry, your input is not found. Please refer to the help.")
        parser.print_help() #print the help link to the parser
        





#######################
##      Main         ##
#######################

#Main of the program, allow to enter our parameters
if __name__=='__main__':
    parser = argparse.ArgumentParser(description= "Process the genomic variant program")
    parser.add_argument("input", help="Choose the file in .txt format")
    parser.add_argument("impact", help="Select the level of putative impact: HIGH,MEDIUM,LOW")
    parser.add_argument("value", help="Select the value to determine the Minor Allele Frequency ")
    parser.add_argument("output", help="Give the name of your file, use the .csv format")
    args = parser.parse_args()

    #save our parameters in variables
    file = args.input
    putative_impact=args.impact
    maf_value = args.value
    output= args.output

     check_format(file)

   
    























