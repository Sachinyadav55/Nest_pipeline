import os
import sys
import glob


fastqc = open('fastqc_data.txt').read().split('>>END_MODULE')
parsed_file = dict()
for vals in fastqc:
    val = vals.split('#') 
    try: 
        parsed_file[val[0]] = val[1]
    except:
        parsed_file[val[0]] = 'PASS'
for vals in parsed_file:
    print (vals)
    #print (parsed_file[vals])
