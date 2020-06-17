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

def write_table(json_MAF,json_data,df,output,selection_ann,annotate_af, variant):
    write_csv("csv1.csv", selection_ann,output,annotate_af, variant)
    
#Add the MAF values and convert in csv file
def write_csv(output_selection,selection_ann,output,annotate_af,variant):
    df = pd.DataFrame(selection_ann)
    df["MAF"]=annotate_af
    df["ID"]=variant
    print(df) #print a data frame with the transcips and the MAF in the terminal
    df.to_csv(output_selection, mode='a', header=0)
#mode a to not errase last row    
    from more_itertools import unique_everseen # Delete possible doublon
    with open(output_selection,'r') as f, open(output,'w') as out_file:
        out_file.writelines(unique_everseen(f))
        


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
    try: #Permit to verify if we have simple or multiple  transcripts
        variable=json_data['snpeff']['ann']['putative_impact']
        if variable==putative_impact:
            selection_ann=json_data['snpeff']['ann']
            #print(selection_ann) #show the json with the data only with a putative_impact : HIGH
            maf_function(df,json_data,variant,output,selection_ann)
    except:
        for i in json_data['snpeff']['ann']:
            for j in json_data['snpeff']['ann']: # For select each transcripts in the ann of the json 
                if j['putative_impact']== putative_impact:
                    selection_ann2=j
                    maf_function(df,json_data,variant,output,selection_ann2)
                    #print(selection_ann2) #show the json with the data only with a putative_impact : HIGH
                else:
                    continue

            






#######################################
##    Functions  Variant  frequency  ##
#######################################

#Function for select only the data with the specific parameters for maf , 
#could be a parameter define by the user but in this program, define directly in the script
def variant_maf(variant,json_data, df,json_MAF,output,selection_ann):
        # print(json_MAF)  #show the json with the data of MAF
        try:
            # false='exac' not in json_MAF
            annotate_af=json_MAF['exac']['af']
            if annotate_af<0.001: #For select af from the preselect data from the variant_putative_impact function
                #print (annotate_af) #show the json with the data only with a putative_impact : HIGH AND MAF<0.001
                write_table(json_MAF,json_data,df,output,selection_ann,annotate_af, variant)  
            elif 'af' not in json_MAF:#Verify if the json had doesn't have a exac 
                annotate_af="NA"
        
                write_table(json_MAF,json_data,df,output,selection_ann,annotate_af, variant)  
        except:
            pass
            
        

#Function for requests for the MAF and verifiy if the connection is good or if the url is correct
def maf_function(df,json_data,variant,output,selection_ann):
    url='http://myvariant.info/v1/variant/'
    url=url+variant
    tmp=0
    parameters= 'fields=exac.alleles,exac.af'
    for i in range(2):
        try: 
            r = requests.get(url,params=parameters,timeout=6.0)   
            json_MAF=r.json()
            #print(json_MAF) #Show the json of the request
            true='success' in json_MAF #Verify if the json is correct or not, if not contnue with others variants ID 
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
  


   
    























