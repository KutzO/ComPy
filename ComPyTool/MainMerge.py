#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  9 10:46:59 2021

@author: kutzolive
"""
#default scripts
import sys
from glob import glob
import logging
import pandas as pd

#own scripts
from .DbManager import DBManager
from .Preparation import DataPreparation 


###Tool to merge databases
class CompToolMerge():
    def __init__(self, argnewdb = False, argolddb = False, argxlsx = False, 
                 argbed = False):
        
        ###Define global variables
        self.newdb = argnewdb
        self.xlsx = argxlsx
        self.bed = argbed
        
        ###Initiale logging tool
        self.compylog = logging.getLogger("ComPy")
        
        ###Start merging
        self.compylog.info("Start merge")
        self.DoMerge(argolddb)

        
        
    ###Function which manages database merging
    def DoMerge(self, argolddb):
        #Define database to be expanded
        classDataBase = DBManager(argolddb)
        self.olddb = classDataBase.pathDB
        self.compylog.info(f"Database to be expanded {classDataBase.pathDB}")
        
        #Define if data derived from other database or .xlsx
        if self.newdb:    
            self.compylog.info(f"New database provided {self.newdb}")
            self.GetData()
            
        elif self.xlsx:
            self.compylog.info(f"Excel files provided {self.xlsx}, {self.bed}")
            lsPathExcel, lsPathBed = self.ControlExcelTables()
            self.GetData(lsPathExcel, lsPathBed)
        
        #Make sure to not include files already present in database
        self.compylog.info("Check bedfiles")
        dicBedID, self.dfBF = self.CheckBed()
        self.compylog.info("Check bamfiles")
        dicBamID = self.CheckFileID("BamInfo", self.dfBamI, "ID")
        self.compylog.info("Check vcf files")
        dicVcfID = self.CheckFileID("VCFInfo", self.dfVI, "ID")
        
        #Adjust IDs of new data to make sure every ID is unique
        self.compylog.info("Adjust dataframes")
        self.AdjustDf(dicBedID, dicBamID, dicVcfID)
        
        #Save adjusted data to database
        self.compylog.info("Add adjusted data to database")
        self.UpdateDB(self.dfBedI.values, "BedInfo")
        self.compylog.info("BedInfo done!")
        self.UpdateDB(self.dfBF.values, "Bedfiles")
        self.compylog.info("Bedfiles done!")
        self.UpdateDB(self.dfBamI.values, "BamInfo")
        self.compylog.info("BamInfo done!")
        self.UpdateDB(self.dfRM.values, "ReadMapping")
        self.compylog.info("ReadMapping done!")
        self.UpdateDB(self.dfRS.values, "ReadStatistiks")
        self.compylog.info("ReadStatistiks done!")
        self.UpdateDB(self.dfQC.values, "QCmetrics")
        self.compylog.info("QCmetrics done!")
        self.UpdateDB(self.dfVI.values, "VCFInfo")
        self.compylog.info("VCFInfo done!")
        self.UpdateDB(self.dfE.values, "Extracted_Variants")
        self.compylog.info("Extracted_Variants done!")
        self.compylog.info("Merge completed!")

      
            
            
    ###Function to check complete input data
    def ControlExcelTables(self):
        lsAllPathExcel = [glob(x)[0] for x  in self.xlsx]
        pathBed = [glob(x)[0] for x in self.bed]
        lsGivenNames = [x.split("/")[-1] for x in lsAllPathExcel]
        lsGivenNames = [x[:x.find(".xlsx")] for x in lsGivenNames]
        lsGivenNames = lsGivenNames
        lsAllTables = ["ReadMapping", "ReadStatistiks", "QCmetrics", 
                        "Extracted_Variants",
                        "BamInfo", "BedInfo", "VCFInfo"]
        lsLost = []
        if len(lsGivenNames) > len(lsAllTables):
            self.compylog.error("To many .xlsx provided!")
            print(
                "Please make sure to insert only following files: "
                +f"{lsAllTables}"
            )
            sys.exit()
        if len(pathBed) == 0:
            self.compylog.error("No bed file .xlsx provided!")
            print(
                "No .bed file was provided. Please provide a fitting .bed file!"
            )
            sys.exit()
        for table in lsAllTables:
            if table not in lsGivenNames:
                lsLost.append(table)
                self.compylog.error(f"Table {table} was not found")
        if len(lsLost) > 0:
            print(f"Following tables are missing: {lsLost}. Please provide!")
            sys.exit()
        else:
            self.compylog.info("All tables were found. Start merging....")
            return lsAllPathExcel, pathBed      





    ###Function to get new data
    def GetData(self, datafiles = False, bedfiles = False):
        if self.newdb:
            self.dfBedI = DBManager.ExtractData(
                "BedInfo", self.newdb, ALL = True
            )
            self.dfBF = DBManager.ExtractData(
                "Bedfiles", self.newdb, ALL = True
            )
            self.dfBamI = DBManager.ExtractData(
                "BamInfo", self.newdb, ALL = True
            )
            self.dfRM = DBManager.ExtractData(
                "ReadMapping", self.newdb, ALL = True
            )
            self.dfRS = DBManager.ExtractData(
                "ReadStatistiks", self.newdb, ALL = True
            )
            self.dfQC = DBManager.ExtractData(
                "QCmetrics", self.newdb, ALL = True
            )
            self.dfVI = DBManager.ExtractData(
                "VCFInfo", self.newdb, ALL = True
            )
            self.dfE = DBManager.ExtractData(
                "Extracted_Variants", self.newdb, ALL = True
            )
        
        else:
            pathBedI = [x for x in datafiles if "BedInfo" in x][0]
            self.dfBedI = pd.read_excel(pathBedI)
            
            self.dfBF = pd.read_excel(bedfiles[0])
            if len(bedfiles) > 1:
                for files in bedfiles[1:]:
                    self.dfBF = self.dfBF.append(
                        pd.read_excel(files), ignore_index = True, sort = False
                    )
            
            pathBamI = [x for x in datafiles if "BamInfo" in x][0]
            self.dfBamI = pd.read_excel(pathBamI)
            
            pathRM = [x for x in datafiles if "ReadMapping" in x][0]
            self.dfRM = pd.read_excel(pathRM)
            
            pathRS = [x for x in datafiles if "ReadStatistiks" in x][0]
            self.dfRS = pd.read_excel(pathRS)
            
            pathQC = [x for x in datafiles if "QCmetrics" in x][0]
            self.dfQC = pd.read_excel(pathQC)
            
            pathVI = [x for x in datafiles if "VCFInfo" in x][0]
            self.dfVI = pd.read_excel(pathVI)
            
            pathE = [x for x in datafiles if "Extracted_Variants" in x][0]
            self.dfE = pd.read_excel(pathE)
        
        #Make sure info BLOB is now type string
        self.dfE["Info"] = self.dfE["Info"].astype(str)




    ###Function to check bed file IDs to avoid duplicate IDs
    def CheckBed(self):
        #Define dic with IDs to be added
        dicBedID = {}
        
        #Define dic with IDs already present in database
        dicBedPresent = {}
        
        #Get information from both bedfile tables
        dfOldBed = DBManager.ExtractData("Bedfiles", self.olddb, ALL = True)
        dfAdd = pd.DataFrame(columns = ["BedID", "Chrom", "Start", "End"])
        lsBedIDNew = list(set(self.dfBF["BedID"]))
        lsBedIDOld = list(set(dfOldBed["BedID"]))
        
        #Exclude bedfiles already present in database
        for newID in lsBedIDNew:
            dfNewSlice = self.dfBF.loc[self.dfBF["BedID"] == newID]
            lsBooAdd = []
            for intID in lsBedIDOld:
                dfOldSlice = dfOldBed.loc[dfOldBed["BedID"] == intID]
                booAdd = self.CheckHelper(dfNewSlice, dfOldSlice)
                lsBooAdd.append(booAdd)
                if not booAdd:
                    dicBedPresent[newID] = intID
            if not False in lsBooAdd:
                dfAdd = dfAdd.append(
                    dfNewSlice, ignore_index = True, sort = False
                )
            else:
                self.compylog.info(f"Bed file {newID} is already present in database!")
        
        #Make sure to insert new unique IDs
        if len(dfAdd.values) > 0:
            lsBedIDNewAdd = list(set(dfAdd["BedID"]))
            for intID in lsBedIDNewAdd:
                booCheck = False
                if intID in lsBedIDOld:
                    for newid in range(1, len(lsBedIDOld) + 1):
                        if newid not in lsBedIDOld \
                         and newid not in dicBedID.keys():
                            dicBedID[intID] = newid
                            booCheck = True
                            break
                else:
                    dicBedID[intID] = intID
                    booCheck = True
                if not booCheck:
                    dicBedID[intID] = len(lsBedIDOld) + 1 
                    lsBedIDOld.append(len(lsBedIDOld) + 1)
                self.compylog.info(f"Bed file ID {intID} now changed to ID {dicBedID[intID]}")
            for intID in lsBedIDNew:
                if intID not in lsBedIDNewAdd:
                    dicBedID[intID] = dicBedPresent[intID]
                    self.compylog.info(f"Bed file ID {intID} now changed to ID {dicBedID[intID]}")
        else:
            #If no new bedfiles are provided create empty dataframe
            for intID in lsBedIDNew:
                dicBedID[intID] = dicBedPresent[intID]
                self.compylog.info(f"Bed file ID {intID} now changed to ID {dicBedID[intID]}")
            self.dfBedI = pd.DataFrame(columns=["BedID", "Date"])
        return dicBedID, dfAdd




    ###Helper function to compare bedfiles
    def CheckHelper(self, dfNew, dfOld):
        dfNew = dfNew[["Chrom", "Start", "End"]]
        dfOld = dfOld[["Chrom", "Start", "End"]]
        if len(dfNew.values) != len(dfOld.values):
            return True
        else:
            for val in zip(dfNew.values, dfOld.values):
                if str(val[0][0]) != str(val[1][0]) \
                    or str(val[0][1]) != str(val[1][1]) \
                    or str(val[0][2]) != str(val[1][2]):
                    return True
        return False 
        
    
    
            
    ###Function to make sure file IDs are always unique
    def CheckFileID(self, table, newDf, col):
        #Extract dataframe with data present in database
        dfOld = DBManager.ExtractData(table, self.olddb, ALL = True)
        
        #Check if data is duplicated
        newDf = self.SortOutPresentSamples(table, newDf, dfOld)
        
        #Identify unique IDs for every new data
        dicID = {}
        lsNewID = list(set(newDf[col]))
        lsOldID = list(set(dfOld[col]))
        for intID in lsNewID:
            booFound = False
            if intID in lsOldID:
                for newID in range(1, len(lsOldID) + 1):
                    if newID not in lsOldID:
                        dicID[intID] = newID
                        self.compylog.info(f"Bed file ID {intID} now changed to ID {newID}")
                        booFound = True
                        lsOldID.append(newID)
                        break
            else:
                dicID[intID] = intID
                self.compylog.info(f"Bed file ID {intID} now changed to ID {intID}")
                lsOldID.append(intID)
                booFound = True
            if not booFound:
                newID = len(lsOldID) + 1
                self.compylog.info(f"Bed file ID {intID} now changed to ID {newID}")
                lsOldID.append(newID)
                dicID[intID] = newID
        return dicID
        
    
    
    
    
    
    ###Helper function to identify duplicate data
    def SortOutPresentSamples(self, table, dfNew, dfOld):
        #Get data from every table
        lsNewMd5 = list(dfNew["Checksum"])
        lsOldMd5 = list(dfOld["Checksum"])
        lsNew = lsNewMd5[:]
        lsOld = lsOldMd5[:]
        if table == "BamInfo":
            lsNewRed = list(dfNew["Reduced"])
            lsOldRed = list(dfOld["Reduced"])
            lsNewSub = list(dfNew["Subsamples"])
            lsOldSub = list(dfOld["Subsamples"])
            lsNewFlag = list(dfNew["Flag"])
            lsOldFlag = list(dfOld["Flag"])
            lsNew = zip(lsNewMd5, lsNewRed, lsNewSub, lsNewFlag)
            lsOld = zip(lsOldMd5, lsOldRed, lsOldSub, lsOldFlag)
        
        #Identify duplicates
        counter = 0
        for i in lsNew:
            if i in lsOld:
                dfNew = dfNew.drop([counter])
                if table == "BamInfo":
                    self.compylog.info(f"Excluded file {self.dfBamI.values[counter]} because duplicate")
                    self.dfBamI = self.dfBamI.drop([counter])
                    
                elif table == "VCFInfo":
                    self.compylog.info(f"Excluded file {self.dfVI.values[counter]} because duplicate")
                    self.dfVI = self.dfVI.drop([counter])
            counter += 1
        
        #Delete dublicates from every table
        self.DeleteDuplicates(table, dfNew)
        return dfNew
        
        
        
        
    ###Helper tool to delete identified duplicates from every table
    def DeleteDuplicates(self, table, dfNew):
        lsReducedIds = list(dfNew["ID"])
        if table == "BamInfo":
            COUNTER = 0
            for intID in self.dfRM["ID"]:
                if intID not in lsReducedIds:
                    self.dfRM = self.dfRM.drop([COUNTER])
                COUNTER += 1
            COUNTER = 0
            for intID in self.dfQC["ID"]:
                if intID not in lsReducedIds:
                    self.dfQC = self.dfQC.drop([COUNTER])
                COUNTER += 1
            COUNTER = 0
            for intID in self.dfRS["ID"]:
                if intID not in lsReducedIds:
                    self.dfRS = self.dfRS.drop([COUNTER])
                COUNTER += 1
        elif table == "VCFInfo":
            COUNTER = 0
            for intID in self.dfE["ID"]:
                if intID not in lsReducedIds:
                    self.dfE = self.dfE.drop([COUNTER])
                COUNTER += 1
        
        
        
        
    ###Function to adjust bed IDs and file IDs in every table
    def AdjustDf(self, dicBedID, dicBamID, dicVcfID):
        for oldID in dicBedID.keys():
            col = "BedID"
            self.dfBedI[col] = self.dfBedI[col].where(
                self.dfBedI[col] != oldID, dicBedID[oldID]
            )
            self.dfBF[col] = self.dfBF[col].where(
                self.dfBF[col] != oldID, dicBedID[oldID]
            )
            self.dfBamI[col] = self.dfBamI[col].where(
                self.dfBamI[col] != oldID, dicBedID[oldID]
            )
            self.dfVI[col] = self.dfVI[col].where(
                self.dfVI[col] != oldID, dicBedID[oldID]
            )
        for oldID in dicBamID.keys():
            col = "ID"
            self.dfBamI[col] = self.dfBamI[col].where(
                self.dfBamI[col] != oldID, dicBamID[oldID]
            )
            self.dfRM[col] = self.dfRM[col].where(
                self.dfRM[col] != oldID, dicBamID[oldID]
            )
            self.dfRS[col] = self.dfRS[col].where(
                self.dfRS[col] != oldID, dicBamID[oldID]
            )
            self.dfQC[col] = self.dfQC[col].where(
                self.dfQC[col] != oldID, dicBamID[oldID]
            )
        for oldID in dicVcfID.keys():
            col = "ID"
            self.dfVI[col] = self.dfVI[col].where(
                self.dfVI[col] != oldID, dicVcfID[oldID]
            )
            self.dfE[col] = self.dfE[col].where(
                self.dfE[col] != oldID, dicVcfID[oldID]
            )
    
    
    
    ###Function to merge adjusted data
    def UpdateDB(self, data, table):
        if len(data) == 0:
            pass
        else:
            if table in ["BamInfo", "VCFInfo", "BedInfo"]:
                for value in data:
                    DBManager.InsertData(value, table, self.olddb)
            else:
                DBManager.InsertData(data, table, self.olddb)

            
    
    
    
  
        
        
        
