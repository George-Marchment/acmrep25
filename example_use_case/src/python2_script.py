#!/usr/bin/env python2.7

import csv
import sys

def kg_to_pounds(kg):
    return kg * 2.20462

def convert_column_to_pounds(input_file):

    column_index = 4
    with open(input_file, 'rb') as csvfile:
        reader = csv.reader(csvfile)
        rows = list(reader)

        for row in rows[1:]:  
            if row: 
                kg_value = float(row[column_index])
                row[column_index] = str(kg_to_pounds(kg_value))

                if(row[9]=="True"):
                    row[9] = 1
                else:
                    row[9] = 0
                
    with open("clean_data.csv", 'wb') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerows(rows)

input_file = sys.argv[1]
convert_column_to_pounds(input_file)

