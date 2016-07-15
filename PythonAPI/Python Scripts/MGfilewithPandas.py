from pandas import Series, DataFrame
import pandas as pd
import os

os.chdir("C:\\Users\\fjanz\\Documents\\PythonAPI")

file1 = pd.read_table("4637812Cdata.tsv", sep='\t')
file2 = pd.read_table("CPER4664901.3MGdata.tsv", sep= '\t')

#check to see if the entire dataset was loaded into the dataframe
file1.tail()
file2.tail()