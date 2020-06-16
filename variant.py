# MULLER CORALIE
#Programming Assignment
# 15 -16 /06/2020

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
import argparse
import requests
from requests import ReadTimeout, ConnectTimeout, HTTPError, Timeout, ConnectionError
import urllib.request
import json
import pandas as pd
from pandas import DataFrame
import csv


#Test
#chr9:g.140695360G>A
#chr10:g.61802416T>C
#chr16:g.28883241A>G
#chr1:g.2938265C>T
#chr6:g.31378977G>A
#chr3:g.97806944T>C
#chr1:g.3496479T>C
#chr3:g.98073591A>T
#chr1:g.152280670G>C
#chr19:g.43268067C>T
#chr11:g.118592881C>T



############################
##    Functions  Files    ##
############################
#Function for check the format and print the help if its not a right format
def check_format(file,putative_impact,output,maf_value):
    variant=""
    if file.lower().endswith(('.txt')):
        parse_file(file,putative_impact,output,maf_value)
    else:
        print("Sorry, your input is not found. Please refer to the help.")
        parser.print_help() #print the help link to the parser

def parse_file(file,putative_impact,output,maf_value):
    tmp=0
    x=-1
    variant=''
    df = pd.read_csv(file, sep="\t", dtype = {"Position" : "object"},header=0)
    for index, row in df.iterrows():
        print(row)
        for i in row:
            tmp+=1
            x+=1
            # x+=1
            if tmp==4:
                variant+=i
                request_variant(variant,df,putative_impact,output,maf_value)
                tmp=0
                x=0
                variant=''
            elif tmp<=4:
                variant+=i
                if tmp==1:
                    variant+=":g."
                elif tmp==3:
                    variant+=">"
              



#####################################
##    Functions  Variant  effects  ##
#####################################

def request_variant(variant,df,putative_impact,output,maf_value): 
    url='http://myvariant.info/v1/variant/'
    url=url+variant
    print(url)
    tmp=0
    parameters= 'fields=snpeff'
    for i in range(1):
        try:
            r = requests.get(url,params=parameters,timeout=6.0)
            json_data=r.json()
            true='success' in json_data
            if true==True:
                print("bad")
            else:
                tmp=2
                continue
        except requests.HTTPError as e:
            print("Sorry!! Connection Error. \n")
            print(str(e))
            continue  
        except requests.ConnectionError as e:
            print("Sorry!! Connection Error. Make sure you are connected to Internet. Technical Details given below.\n")
            print(str(e))            
            continue
        except requests.Timeout as e:
            print("Sorry!! Timeout Error")
            print(str(e))
            continue
        except requests.RequestException as e:
            print("Sorry!! General Error")
            print(str(e))
            continue
        except KeyboardInterrupt:
            print("Someone closed the program")
    if tmp==2:
        variant_putative_impact(variant, df,json_data,putative_impact,output,maf_value)



def variant_putative_impact(variant, df,json_data,putative_impact,output,maf_value):
    num=0
    try:
        variable=json_data['snpeff']['ann']['putative_impact']
        if variable==putative_impact:
            num=1
            print("work")
            # add_annotation(json_data, num,df,output)
            maf_function(df,variant,output,maf_value)
    except:
        for i in json_data['snpeff']['ann']:
            if i['putative_impact']== putative_impact:
                num=2
            else:
                continue
        if num==2:
            # add_annotation(json_data, num,df)
            maf_function(df,variant,output,maf_value)
            print(2)


def add_annotation(json_data, num,df,output):
    print("annotate")



#######################################
##    Functions  Variant  frequency  ##
#######################################
def variant_maf(variant, df,json_MAF,output,maf_value):
        print(json_MAF)
        try:
            variable=json_MAF['exac']['af']
            print (variable)
            # maf_value=int(maf_value)
            if variable<0.001:
                print("right")
                write_value(json_MAF)    
        except:
            print("Don't have af value for this variant")
        


def maf_function(df,variant,output,maf_value):
    url='http://myvariant.info/v1/variant/'
    url=url+variant
    print(url)
    tmp=0
    parameters= 'fields=exac.alleles,exac.af'
    for i in range(2):
        try: 
            r = requests.get(url,params=parameters,timeout=6.0)   
            json_MAF=r.json()
            true='success' in json_MAF
            if true==True:
                print("bad")
            else:
                tmp=2
                continue
        except requests.HTTPError as e:
            print("Sorry!! Connection Error. \n")
            print(str(e))  
        except requests.ConnectionError as e:
            print("Sorry!! Connection Error. Make sure you are connected to Internet. Technical Details given below.\n")
            print(str(e))            
            continue
        except requests.Timeout as e:
            print("Sorry!! Timeout Error")
            print(str(e))
            continue
        except requests.RequestException as e:
            print("Sorry!! General Error")
            print(str(e))
            continue
        except KeyboardInterrupt:
            print("Someone closed the program")
    if tmp==2:
        variant_maf(variant, df,json_MAF,output,maf_value)
    
    
#######################################
##    Functions write outputfile     ##
#######################################

def write_value(json_MAF):
    print("write")







# def fusion_file():




    

#######################
##      Main         ##
#######################

#Main of the program, allow to enter our parameters
if __name__=='__main__':
    parser = argparse.ArgumentParser(description= "Process the genomic variant program")
    parser.add_argument("input", help="Choose the file in .txt format")
    parser.add_argument("impact", help="Select the level of putative impact: HIGH,MEDIUM,LOW")
    parser.add_argument("value", help="Select the value to determine the Minor Allele Frequency you want ")
    parser.add_argument("output", help="Give the name of your file, use the .csv format")
    args = parser.parse_args()

    #save our parameters in variables
    file = args.input
    putative_impact=args.impact
    maf_value = args.value
    output= args.output
    check_format(file,putative_impact,output,maf_value)
  


   
    























