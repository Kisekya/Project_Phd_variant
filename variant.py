# MULLER CORALIE
#Programming Assignment
# 15 -16 /06/2020

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
import argparse
import requests
from requests import ReadTimeout, ConnectTimeout, HTTPError, Timeout, ConnectionError
import urllib.response
import json
# import fileinput
import re
import pandas as pd
from pandas import DataFrame


#Test
# chr9:g.140695360G>A
#chr16:g.28883241A>G

############################
##    Functions  Files    ##
############################
#Function for check the format and print the help if its not a right format
def check_format(file):
    variant=""
    if file.lower().endswith(('.txt')):
        parse_file(file)
        print("txt")
    else:
        print("Sorry, your input is not found. Please refer to the help.")
        parser.print_help() #print the help link to the parser

def parse_file(file):
    lines=[]
    tmp=0
    x=-1
    variant=''
    df = pd.read_csv(file, sep="\t", dtype = {"Position" : "object"},header=0)
    for index, row in df.head().iterrows():
        for i in row:
            print(i)
            tmp+=1
            x+=1
            # x+=1
            if tmp==4:
                variant+=i
                request_variant(variant)
                tmp=0
                x=0
                variant=''
            elif tmp<=4:
                variant+=i
                if tmp==1:
                    variant+=":g."
                elif tmp==3:
                    variant+=">"
                







    

#     with open(file) as f:
#         file=f.read().split("\n")
#         print(file)
#         i=0
#         for line in file:
#             print(file[i])
#             i+=1 TO FINISH



#####################################
##    Functions  Variant  effects  ##
#####################################

def request_variant(variant): 
    url='http://myvariant.info/v1/variant/'
    url=url+variant
    print(url)
    parameters= 'fields=snpeff'
    for i in range(2):
        try:
            r = requests.get(url,params=parameters,timeout=6.0)
            print(r.url)
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
    json_data=r.json()
    num=0
    try:
        variable=json_data['snpeff']['ann']['putative_impact']
        if variable=="HIGH":
            num=1
            print("work")
        # add_annotation(json_data, num)
    except:
        for i in json_data['snpeff']['ann']:
            if i['putative_impact']== "HIGH":
                num=2
            else:
                continue
        if num==2:
        # add_annotation(json_data, num)
            print(2)

# def add_annotation(json_data, num):
#     if num ==1:
#         result = {}
#         for key in json_data.keys():
#             if not isinstance(json_data[key], dict):
#                 result[key] = json_data[key]
#             else:
#                 result.update(add_annotation(json_data[key]))
#         return result
#     elif num==2:
#         print("nnnn")


# json_datas = json.load()
# json_list = [j[1][0] for j in json_datas.iterrows()]
# parsed_list = [parse_nested_json(j) for j in json_list]
# result = pd.DataFrame(parsed_list)
# print (result)
# # result.to_csv("my_csv_file.csv", index=False)
    

    

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
    # request_variant()
    check_format(file)

   
    























