



#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug  4 10:42:16 2021

@author: kutzolive
"""


###Import needed packages
#External packages
import logging 

#Own scripts
from .DbManager import DBManager
from .Preparation import DataPreparation



"""
KWARGS: argFileSum, argDatabase, argBedfile, argSubsample, 
        argReduce, dtime, strVersion, lsCompatibleVersions
"""
def CompToolDelete(tool, **kwargs):
    
    ###Prepare logging
    compylog = logging.getLogger("ComPy")
    
    
    ###Clean up database
    if tool=="clean":
        """
        Workflow
        1) Get database path
        2) Clean database
        """
        
        """
        1) Get database path
        """
        classDataBase = DBManager(kwargs["argDatabase"])                        #Class defined in script DBManager.py

        
        """
        2) Clean database
        """
        #Function defined in script DBManager.py
        DBManager.DelData(
            dbpath = classDataBase.pathDB, booClean = True, 
            lsCompVers = kwargs["lsCompatibleVersions"]
        )
    



    
    ###Delete files from database
    if tool=="files":
        """
        Workflow:
        1) Get database path
        2) Prepare data
        3) Delete data
        """
        
        """
        1) Get database path
        """
        classDataBase = DBManager(kwargs["argDatabase"])                        #Class defined in script DBManager.py

        
        """
        2) Prepare data
        """
        #Class defined in script Preparation.py
        classPrep = DataPreparation(
            "delete", lsbamid = kwargs["argbamid"], 
            lsvcfid = kwargs["argvcfid"],  booDB = classDataBase.booDB, 
            DBpath= classDataBase.pathDB, nameTable = kwargs["argNameTable"], 
            dtime= kwargs["dtime"]
        )
   
        """
        3) Delete data
        """
        #Funktion defined in script DBManager.py
        DBManager.DelData(
            dbpath = classDataBase.pathDB, 
            bamID = classPrep.dicIDs["bam"], 
            vcfID = classPrep.dicIDs["vcf"]
        )