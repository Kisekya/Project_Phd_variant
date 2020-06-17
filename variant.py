# MULLER CORALIE
#Programming Assignment
# 15 -16 /06/2020

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import sys
import argparse
import requests
from requests import ReadTimeout, ConnectTimeout, HTTPError, Timeout, ConnectionError
import json
import pandas as pd
from pandas import DataFrame
import csv


#Id to test in Postman for security
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

#Function for check the format and print the help for parameters if its not a right format
def check_format(file,putative_impact,output):
    if file.lower().endswith(('.txt')): #verify the format of the input file, must be a txt 
        parse_file(file,putative_impact,output)
    else:
        print("Sorry, your input is not found. Please refer to the help.")
        parser.print_help() #print the help link to the parser


#Function for the creation of the variant id from the table, ex :chr9:g.140695360G>A:
def parse_file(file,putative_impact,output):
    #Initialize
    tmp=0
    x=-1
    variant=''
    df = pd.read_csv(file, sep="\t", dtype = {"Position" : "object"},header=0) 
    #print(df) #print the table
    for index, row in df.iterrows(): # permit to have key and values by rows
        #print(row) #print each row
        for i in row: # permit to have one by one the values to save in a variable
            tmp+=1
            x+=1
            if tmp==4: #for the four columns: to be sure we selecte one complete ID
                variant+=i
                #print(variant) #Print the ID of the variant
                request_variant(variant,df,putative_impact,output)
                tmp=0
                x=0
                variant=''
            elif tmp<=4: #add values and parameters for the creation of the id at specific position 
                variant+=i
                if tmp==1:
                    variant+=":g."
                elif tmp==3:
                    variant+=">"
              



#######################################
##    Functions write outputfile     ##
#######################################

def write_table(json_MAF,json_data,df,output,selection_ann,annotate_af):
    write_csv("csv1.csv", selection_ann)
    write_csv("csv2.csv",json_data)
    write_csv("af.csv",annotate_af)


def write_csv(output,selection):
    data_to_file = open(output, 'w', newline='')
    csv_writer = csv.writer(data_to_file,delimiter=";")
    # Counter variable used for writing  
    # headers to the CSV file 

    #test append
    count = 0
    
    ann_data = selection
    for emp in ann_data: 
        if count == 0: 
        
            # Writing headers of CSV file 
            header = emp.keys() 
            csv_writer.writerow(header) 
            count += 1
    
        # Writing data of CSV file 
        csv_writer.writerow(emp.values()) 
    
    data_file.close() 

#####################################
##    Functions  Variant  effects  ##
#####################################

#Function for requests and verifiy if the connection id good or if the url is correct
def request_variant(variant,df,putative_impact,output): 
    url='http://myvariant.info/v1/variant/'
    url=url+variant #specific variant ID for each request
    # print(url)
    tmp=0
    parameters= 'fields=snpeff' #parameters of the url
    for i in range(1): #Try two time the request
        try:
            r = requests.get(url,params=parameters,timeout=6.0) #request line with parameters and the timeout limits
            json_data=r.json()  #stock the json data from the request
            #print(json_data) #print the json 
            true='success' in json_data #Verify if the json is correct or not, if not contnue with others variants ID 
            if true==True:
                continue
            else:
                tmp=2
                continue
                #### Excepts for Error Gestion ########
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
        variant_putative_impact(variant, df,json_data,putative_impact,output)
        

#Function for select only the data with the specific parameters for the putative_impact, parameter define by the user
def variant_putative_impact(variant, df,json_data,putative_impact,output):
    try: #Permit to verify if we have multiple trand
        variable=json_data['snpeff']['ann']['putative_impact']
        if variable==putative_impact:
            selection_ann=json_data['snpeff']['ann']
            #print(selection_ann) #show the json with the data only with a putative_impact : HIGH
            maf_function(df,json_data,variant,output,selection_ann)
    except:
        for i in json_data['snpeff']['ann']:
            for j in json_data['snpeff']['ann']:
                #print(j)
                if j['putative_impact']== putative_impact:
                    selection_ann2=j
                    maf_function(df,json_data,variant,output,selection_ann2)
                    print(selection_ann2) #show the json with the data only with a putative_impact : HIGH
                else:
                    continue

            






#######################################
##    Functions  Variant  frequency  ##
#######################################
def variant_maf(variant,json_data, df,json_MAF,output,selection_ann):
        # print(json_MAF)
        try:
            annotate_af=json_MAF['exac']['af']
            # print (annotate_af)
            if annotate_af<0.001:
                # print("right")
                write_table(json_MAF,json_data,df,output,selection_ann,annotate_af)    
        except:
            print("Don't have af value for this variant")
        


def maf_function(df,json_data,variant,output,selection_ann):
    url='http://myvariant.info/v1/variant/'
    url=url+variant
    # print(url)
    tmp=0
    parameters= 'fields=exac.alleles,exac.af'
    for i in range(2):
        try: 
            r = requests.get(url,params=parameters,timeout=6.0)   
            json_MAF=r.json()
            true='success' in json_MAF
            if true==True:
                continue
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
        variant_maf(variant,json_data, df,json_MAF,output,selection_ann)
        
    
    


    

#######################
##      Main         ##
#######################


#Main of the program, allow to enter our parameters
if __name__=='__main__':
    parser = argparse.ArgumentParser(description= "Process the genomic variant program")
    parser.add_argument("input", help="Choose the file in .txt format")
    parser.add_argument("impact", help="Select the level of putative impact: HIGH,MEDIUM,LOW")
    parser.add_argument("output", help="Give the name of your file, use the .csv format")
    args = parser.parse_args()

    #save our parameters in variables
    file = args.input
    putative_impact=args.impact
  
    output= args.output
    check_format(file,putative_impact,output)
  


   
    























