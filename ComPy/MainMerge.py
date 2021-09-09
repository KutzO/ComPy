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

class CompToolMerge():
    def __init__(self, argnewdb = False, argolddb = False, argxlsx = False, 
                 argbed = False):
        self.newdb = argnewdb
        self.xlsx = argxlsx
        self.bed = argbed
        
        self.DoMerge(argolddb)
        
        
        
        
        
        
        
    def DoMerge(self, argolddb):
        classDataBase = DBManager(argolddb)
        self.olddb = classDataBase.pathDB
        if self.newdb:    
            self.GetData()
            
        elif self.xlsx:
            lsPathExcel, lsPathBed = self.ControlExcelTables()
            self.GetData(lsPathExcel, lsPathBed)
        
        dicBedID, self.dfBF = self.CheckBed()
        dicBamID = self.CheckFileID("BamInfo", self.dfBamI, "ID")
        dicVcfID = self.CheckFileID("VCFInfo", self.dfVI, "ID")
        
        self.AdjustDf(dicBedID, dicBamID, dicVcfID)
        self.UpdateDB(self.dfBedI.values, "BedInfo")
        self.UpdateDB(self.dfBF.values, "Bedfiles")
        self.UpdateDB(self.dfBamI.values, "BamInfo")
        self.UpdateDB(self.dfRM.values, "ReadMapping")
        self.UpdateDB(self.dfRS.values, "ReadStatistiks")
        self.UpdateDB(self.dfQC.values, "QCmetrics")
        self.UpdateDB(self.dfVI.values, "VCFInfo")
        self.UpdateDB(self.dfE.values, "Extracted_Variants")
        # self.UpdataDB(self.dfSO, "Sorted_out")
      
            
            
            
    def ControlExcelTables(self):
        compylog = logging.getLogger("ComPy")
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
            compylog.error("To many .xlsx provided!")
            print(
                "Please make sure to insert only following files: "
                +f"{lsAllTables}"
            )
        if len(pathBed) == 0:
            print(
                "No .bed file was provided. Please provide a fitting .bed file!"
            )
        for table in lsAllTables:
            if table not in lsGivenNames:
                lsLost.append(table)
                compylog.error(f"Table {table} was not found")
        if len(lsLost) > 0:
            print(f"Following tables are missing: {lsLost}. Please provide!")
            sys.exit()
        else:
            compylog.info("All tables were found. Start merging....")
            return lsAllPathExcel, pathBed      



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
            # self.dfSO = DBManager.ExtractData("Sorted_out", self.olddb, ALL = True)
        
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
            
            # pathSO = [x for x in datafiles if "Sorted_out" in x][0]
            # self.dfSO = pd.read_excel(pathSO)
            
        self.dfE["Info"] = self.dfE["Info"].astype(str)




        
    def CheckBed(self):
        dicBedID = {}
        dfOldBed = DBManager.ExtractData("Bedfiles", self.olddb, ALL = True)
        dfAdd = pd.DataFrame(columns = ["BedID", "Chrom", "Start", "End"])
        lsBedIDNew = list(set(self.dfBF["BedID"]))
        lsBedIDOld = list(set(dfOldBed["BedID"]))
        for intID in lsBedIDOld:
            dfOldSlice = dfOldBed.loc[dfOldBed["BedID"] == intID]
            for newID in lsBedIDNew:
                dfNewSlice = self.dfBF.loc[self.dfBF["BedID"] == newID]
                booAdd = self.CheckHelper(dfNewSlice, dfOldSlice)
                # print(booAdd)
                # sys.exit()
                if booAdd:
                    dfAdd = dfAdd.append(
                        dfNewSlice, ignore_index = True, sort = False
                    )
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
            for intID in lsBedIDNew:
                if intID not in lsBedIDNewAdd:
                    dicBedID[intID] = intID
        else:
            for intID in lsBedIDNew:
                dicBedID[intID] = intID
            self.dfBedI = pd.DataFrame(columns=["BedID", "Date"])
        return dicBedID, dfAdd



        
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
        
    
            
    
    def CheckFileID(self, table, newDf, col):
        dfOld = DBManager.ExtractData(table, self.olddb, ALL = True)
        newDf = self.SortOutPresentSamples(table, newDf, dfOld)
        dicID = {}
        lsNewID = list(set(newDf[col]))
        lsOldID = list(set(dfOld[col]))
        for intID in lsNewID:
            booFound = False
            if intID in lsOldID:
                for newID in range(1, len(lsOldID) + 1):
                    if newID not in lsOldID:
                        dicID[intID] = newID
                        booFound = True
                        lsOldID.append(newID)
                        break
            else:
                dicID[intID] = intID
                lsOldID.append(intID)
                booFound = True
            if not booFound:
                newID = len(lsOldID) + 1
                lsOldID.append(newID)
                dicID[intID] = newID
        return dicID
        
    
    
    
    
    
    
    def SortOutPresentSamples(self, table, dfNew, dfOld):
        #lsIDNew = list(dfNew["ID"])
        lsNewMd5 = list(dfNew["Checksum"])
        lsOldMd5 = list(dfOld["Checksum"])
        lsNew = lsNewMd5[:]
        lsOld = lsOldMd5[:]
        # print(dfNew)
        # print(dfOld)
        if table == "BamInfo":
            lsNewRed = list(dfNew["Reduced"])
            lsOldRed = list(dfOld["Reduced"])
            lsNewSub = list(dfNew["Subsamples"])
            lsOldSub = list(dfOld["Subsamples"])
            lsNew = zip(lsNewMd5, lsNewRed, lsNewSub)
            lsOld = zip(lsOldMd5, lsOldRed, lsOldSub)
        
        counter = 0
        # print(table)
        # print(lsNew)
        # print(lsOld)
        for i in lsNew:
            if i in lsOld:
                dfNew = dfNew.drop([counter])
                if table == "BamInfo":
                    self.dfBamI = self.dfBamI.drop([counter])
                elif table == "VCFInfo":
                    self.dfVI = self.dfVI.drop([counter])
            counter += 1
        self.DeleteDuplicates(table, dfNew)
        return dfNew
        
        
    def DeleteDuplicates(self, table, dfNew):
        # print(dfNew)
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
        
        
        
        
        
    def AdjustDf(self, dicBedID, dicBamID, dicVcfID):
        for oldID in dicBedID.keys():
            col = "BedID"
            self.dfBedI[col] = self.dfBedI[col].where(
                self.dfBedI[col] != oldID, dicBedID[oldID]
            )
            # sys.exit()
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
            # self.dfSO[col] = self.dfSO.where(self.dfSO[col] != oldID, dicVcfID[oldID])
    
    
    
 
    def UpdateDB(self, data, table):
        if len(data) == 0:
            pass
        else:
            if table in ["BamInfo", "VCFInfo", "BedInfo"]:
                for value in data:
                    DBManager.InsertData(value, table, self.olddb)
            else:
                DBManager.InsertData(data, table, self.olddb)
        # if table == "BamInfo":
            # pass
            
    
    
    
  
        
        
        
