# -*- coding: utf-8 -*-
"""
Code by: Palash Patole
MoMpy Project: External_files

Date of creation: Tue May 26 08:44:17 2020

Version: ver1 (almost exactly same as Base, but in public domain - credentials for SQL DB are removed)

Description:
        1. Reads the .pckl files for creation of databases [step 0]
        2. Creates a database, if it does not exist [step 1]
        3. Creates 3 tables in the database, if they do not exist [step 2]
        4. Fills up data in the generated tables [step 3]
        5. User should provide numbers for the required steps in "steps" array
        6. Output validated by executing SQL queries on the generated DB in MySQL Workbench 
"""

###########################################
#### Importing the required modules  ######   
###########################################

import pickle as pk
import sys
import mysql.connector as SQLc
import numpy as np
import math

sys.path.append('../') # access to the entire MoMpy directory
  
#####################################
########### User inputs #############    
#####################################

steps = np.array([0,2,3])   # 0 - reading the data files
                            # 1 - creating a database if does not exist
                            # 2 - creating all tables in the existing database
                            # 3 - filling up all existing tables of the exsiting database

# Paths to the three .pckl files
Path1 = '../Output/EarthRepeats_data1.pckl'
Path2 = '../Output/EarthRepeats_data2.pckl'
Path3 = '../Output/EarthRepeats_data3.pckl'

# SQL Database details
HostName = "" # <enter host-name here, such as the name of local machine>
UserName = "" # <enter the user name>
Password = "" # <enter the associated passoword>
DataBaseName = "EarthRepeatOrbits" 

# Table names
tableNames = np.array(['ERO_set1','ERO_set2','ERO_set3'])
              
######################################
######## Computations ################
######################################

for step in steps:
    if step == 0: # Extracting the data
        file1 = open(Path1,'rb')
        Results1 = pk.load(file1)
        file1.close()
        
        file2 = open(Path2,'rb')
        Results2 = pk.load(file2)
        file2.close()
        
        file3 = open(Path3,'rb')
        Results3 = pk.load(file3)
        file3.close
        
        # Creating a mega-storage of all the results extracted
        Storage = np.array([Results1,Results2,Results3])
    
    elif step == 1: # Database creation
        # Connecting to a host
        DB_connection = SQLc.connect(
                    host = HostName,
                    user = UserName,
                    passwd= Password)
                
        # Creating a cursor to execute SQL queries
        Cursor = DB_connection.cursor()
        
        # Creating a database
        Query = "CREATE DATABASE " + DataBaseName
        Cursor.execute(Query)
        print('Database creation successful!')

    elif step == 2: # Creation of tables
        # Error handling
        if 'Storage' in locals():
            print('Creating tables in the database for the imported results.')
        else:
            print('Input data missing to create tables.')
            sys.exit
            
        # Connecting to an existing database
        DB_connection = SQLc.connect(
                    host = HostName,
                    user = UserName,
                    passwd= Password,
                    database = DataBaseName)    
        
        # Creating a cursor to execute SQL queries
        Cursor = DB_connection.cursor()
        
        # Creating a table
        assert(tableNames.size==3),"Names for three tables are required"
        for tN in tableNames:
            Query = "CREATE TABLE " + tN + " (Sol_index int, Orbit_revolutions int, Days int, Eccentricity double, Altitude double, Inclination double, Orbit_type varchar(100),PRIMARY KEY (Sol_index))"
            Cursor.execute(Query)
            
        print('Creation of three tables is successful!')
        
    elif step ==3:
        
        # Connecting to an existing database
        DB_connection = SQLc.connect(
                    host = HostName,
                    user = UserName,
                    passwd= Password,
                    database = DataBaseName)    
        
        # Creating a cursor to execute SQL queries
        Cursor = DB_connection.cursor()       
        
        # Fiiling up the data in generated tables using ResultsX variables
        for tableIndex in range (tableNames.size):
            Result = Storage[tableIndex]
            Dimensions = Result.shape
            Solution_count = 0
            Rows_read = 0
            for counter1 in range(Dimensions[0]):
                for counter2 in range(Dimensions[1]):
                    if (not (math.isnan(Result[counter1,counter2,4])) or not(math.isnan(Result[counter1,counter2,5]))):
                        Solution_count += 1

                        Query_1 = "INSERT INTO "+ tableNames[tableIndex] + " (Sol_index,Orbit_revolutions,Days,Eccentricity, Altitude, Inclination, Orbit_type) VALUES("
                                    
                        # Storing the solution count
                        Query_2 = str(Solution_count)+ "," 
                    
                        # Storing the j and k values
                        Query_3 = str(int(Result[counter1,counter2,0])) + "," # j
                        Query_4 = str(int(Result[counter1,counter2,1])) + "," # k
                        
                        # Storing the eccentricity, altitude, and inclination
                        Query_5 = str(Result[counter1,counter2,2]) + "," # eccentricity
                        if tableIndex == 2:
                            incl_col = 3
                            Query_7 = str(Result[counter1,counter2,3]) + "," # inclination
                            if not (math.isnan(Result[counter1,counter2,4])): # altitude
                                Query_6 = str(Result[counter1,counter2,4]) + ","
                            else:
                                Query_6 = str(Result[counter1,counter2,5]) + ","
                            
                        else:    
                            Query_6 = str(Result[counter1,counter2,3]) + "," # altitude    
                            if not (math.isnan(Result[counter1,counter2,4])): # inclination
                                Query_7 = str(Result[counter1,counter2,4]) + ","
                                incl_col = 4
                            else:
                                Query_7 = str(Result[counter1,counter2,5]) + ","
                                incl_col = 5
                        
                        # Deciding the orbit type and storing it
                        if (Result[counter1,counter2,2] == 0.0):
                            Query_8 = "'Circular')"
                            if (Result[counter1,counter2,incl_col] >= 89.0) & (Result[counter1,counter2,incl_col] <= 91.0):
                                Query_8 = "'Circular_near_polar')"
                            elif (Result[counter1,counter2,incl_col] >= -1.0) & (Result[counter1,counter2,incl_col] <= 1.0):
                                Query_8 = "'Circular_near_equatorial')"
                        else:
                            Query_8 = "'Elliptical')"
                            if (Result[counter1,counter2,incl_col] >= 89.0) & (Result[counter1,counter2,incl_col] <= 91.0):
                                Query_8 = "'Elliptical_near_polar')"
                            elif (Result[counter1,counter2,incl_col] >= -1.0) & (Result[counter1,counter2,incl_col] <= 1.0):
                                Query_8 = "'Elliptical_near_equatorial')"
                    
                        Cursor.execute(Query_1 + Query_2+ Query_3 + Query_4 + Query_5 + Query_6 + Query_7 + Query_8)

                        DB_connection.commit() # Committing the queries
                        
                        if (not (math.isnan(Result[counter1,counter2,4]))) & (not(math.isnan(Result[counter1,counter2,5]))):
                            print('Two solutions at the same altitude value, skipped one of them in the table.')
                    
                    Rows_read += 1
                        
            print('Solutions found: ', Solution_count, ' from rows read: ',Rows_read, ' and stored in the table viz. ',tableNames[tableIndex],'.')
        