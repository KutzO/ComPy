import random
import sys
import pandas as pd
import os
from .DbManager import DBManager
import logging
import csv
from glob import glob
import getpass


class DataPreparation:
    def __init__(self, tool, bam=False, DBpath=False, booDB=False, bed=False, 
                 subsample=False, dtime=False, outputpath=False, vcf = False, 
                 ReduceBed = False, version = False, checksum = False, 
                 nameTable = False, vcfPlot = False, FileClass = False,
                 intID = False, Datatype = False):
        
        self.dbpath = DBpath
        self.booDB = booDB
        self.subsample = subsample
        self.version = version
        self.dtime = dtime
        
        if FileClass:
            self.FileClass = FileClass
        else:
            self.FileClass = "Default"
        
        #Initiale logging tool
        self.comptoollog = logging.getLogger("ComparisonTool")
        
        
        
        
        if tool == "compare":
            ##Identify outputpath
            self.outputpathTmp, self.outputpath = self.FindOutputPath(
                                                        outputpath
                                                  )
            self.PrepareCompare(
                bed, ReduceBed, bam, vcf, checksum, nameTable
            )
            
        elif tool == "extract":
            ##Identify outputpath
            #self.outputpathTmp, self.outputpath = self.FindOutputPath(
            #                                            outputpath
            #                                      )
            self.PrepareExtraction(
                nameTable, intID
            )
            
        elif tool == "delete":
            self.PrepareDeletion(
                Datatype, intID, nameTable
            )
            
            
            
            
    """
    Prepare deletion of data from database
    """
    def PrepareDeletion(self, Datatype, intID, nameTable):
        if nameTable:
            df = pd.read_csv(nameTable)
            if len(list(df.columns)) == 8:
                self.vcf = True
            elif len(list(df.columns)) == 10:
                self.vcf = False
            self.ID = df["ID"].values
            
        else:
            self.ID = intID
            print(self.ID)
            if Datatype == "bam":
                self.vcf = False
            elif Datatype == "vcf":
                self.vcf = True


    """
    Prepare data extraction from database
    """
    
    def PrepareExtraction(self, nameTable, intID):
        if nameTable:
            dfNameData = pd.read_csv(nameTable)
            self.ID = dfNameData["ID"].values
            
        else:
            self.ID = intID
            # self.CollectInfoForExtraction(
            #     bed, ReduceBed, checksum, FileClass
            # )
            
        
        
    
    # def CollectInfoForExtraction(self, bed, ReduceBed, checksum, FileClass):
    #         try:
    #             int(bed)
    #             self.bedid = (bed,)
    #         except:
    #             self.dictTargets = self.SliceBedFile(bed)
    #             self.bedid = self.CheckBed()
    #             self.bedid = (self.bedid,)
    #         try:
    #             glob(checksum[0])
    #             from .MainCompare import CheckSumGenerator
    #             tmpchecksum = [CheckSumGenerator(file) for file in checksum]
    #             self.checksum = [x[0] for x in tmpchecksum]
    #         except:
    #             self.checksum = checksum
    #         if ReduceBed:
    #             self.Reduce = 1
    #         else:
    #             self.Reduce = 0
            
                
                

    



    """
    Prepare comparison of files
    """
    
    
        
    def PrepareCompare(self, bed, ReduceBed, bam, vcf, checksum, nameTable):
        ##Identify targets present in given BED file
        #Optional: Reduce number of targets
        if bed:
            self.dictTargets = self.SliceBedFile(bed)
            self.bedid = self.CheckBed()
        else:
            self.dictTargets = False
            self.bedid = False
            
        if ReduceBed:
            self.dictTargets = self.ReduceBedTargets()
            self.Reduce = 1
        else:
            self.Reduce = 0        
        
        ##Extract .bam names and pathes from parsed arguments
        self.dicIDs = {}
        if bam:
            self.bamnamestodb, self.allbamnames = self.LookAtArgsCompare(
                bam, "bam", checksum, nameTable
           )
        else:
            self.bamnames = False
            self.allbamnames = False
          
        ##Extract .vcf names and pathes from parsed arguments
        if vcf:
            self.vcfnames, self.allvcfnames = self.LookAtArgsCompare(
                vcf, "vcf", checksum, nameTable
            )
            self.booVCF = True
        else:
            self.booVCF = False
            self.allvcfnames = False



    
    ###Define outputpath and create temp folders
    def FindOutputPath(self, outputpath):
        ##Define path
        UserID = getpass.getuser()
        if outputpath:
            if outputpath[:2] == ".." or outputpath[:2] == "./":
                strOutput = os.getcwd() + "/"+outputpath
                if strOutput[-1] != "/":
                    strOutput = strOutput + "/"
            else:
                strOutput = outputpath
                if strOutput[-1] != "/":
                    strOutput = strOutput + "/"  
            strOutputTmp = strOutput + f"tmp/{self.dtime}/"
            strOutput = strOutput + f"Result_{self.dtime}/"
        else:
            strOutputTmp = str(
                f"/home/{UserID}/ComparisonTool/tmp/{self.dtime}/"
            )
            strOutput = str(
                f"/home/{UserID}/ComparisonTool/Result_{self.dtime}/"
            )
        ##Create folders
        
        self.comptoollog.info(f"Create output folder (tmp): {strOutputTmp}")
        if not os.path.exists(strOutputTmp):
            os.makedirs(strOutputTmp)
        
        self.comptoollog.info(f"Create output folder (results): {strOutput}")
        if not os.path.exists(strOutput):
            os.makedirs(strOutput)
        
        strOutputBigComp = strOutput+"BigCompare/"
        self.comptoollog.info("Create output folder "
                              +f"(DBcomp): {strOutputBigComp}"
        )
        if not os.path.exists(strOutputBigComp):
            os.makedirs(strOutputBigComp)
        
        return strOutputTmp, strOutput





    ###Extract targetdata from BED file: Chromosome, Start, End
    def SliceBedFile(self, bed):
        bedimport = pd.read_csv(bed, sep="\t", header=None, low_memory=False)
        start = 0
        #Check-up if a header is present --> if yes: all info extraction starts from position 1 (Position 1 == row 1)
        if bedimport.values[0][2] != int(bedimport.values[0][2]):
            start = 1

        #Identify chromosomes
        lsChromInt = []
        lsChromStr = []
        for chromos in set(bedimport[0].values[start:]):
            try:
                lsChromInt.append(int(chromos))
            except:
                if "chr" in chromos:
                    chromos = chromos[3:]
                try:
                    lsChromInt.append(int(chromos))
                except Exception as e:
                    self.comptoollog.info(
                        f"A chromosome in .bed file was not integer: {chromos}"
                    )
                    lsChromStr.append(str(chromos))
        lsChromInt.sort()
        lsChromStr.sort()
        lsChrom = [str(x) for x in lsChromInt+lsChromStr]
        
        #Identify targets and assign each to the corresponding chromosome
        dictTarget = {}
        for key in lsChrom:
            lsValues = []
            for value in bedimport.values[start:]:
                if str(value[0]) == key:           
                    lsValues.append([int(value[1]), int(value[2])])
            dictTarget[key]=lsValues
                
        return dictTarget
        
    

    
        
    """
    This functions calculates the BED id!
    It checks:
        1) Is the database containing any BED file information?
            1.1) Are the current BED file targets already saved?
                --> if yes, which ID is associated?
                --> if no, which ID should be used next?
    """
    def CheckBed(self):
        if not self.booDB:      #If a new database has to be created, the starting .bed ID is 1
            bedid = 1
            self.AddBedToDatabase(bedid)
            return bedid
        
        else:
            allbedids = DBManager.ExtractData(
                "Bedfiles", self.dbpath, boobed=True
            )  
            
            lsAllids = list(set([x[0] for x in allbedids]))    #Because db extraction returns (ID,) for each .bed ID in "BamInfo" table
            
            ##Convert current .bed dictionary to list for comparison
            lsNewBed = []   
            for chromosome in self.dictTargets.keys():
                for targets in self.dictTargets[chromosome]:
                    lsNewBed.append(
                        [chromosome,int(targets[0]), int(targets[1])]
                    )
            # print(lsNewBed)
            ##Iterate through all .bed ID found in database
            booIterate = True
            if len(lsAllids) == 0:
                booCheck = False
                lsAllids = [0]
                booIterate = False
            if booIterate:
                for bedid in lsAllids:
                    booCheck = True
                    
                    #Extract all targets from db associated with .bed ID
                    bedinfo = DBManager.ExtractData(
                        "Bedfiles",self.dbpath, bedid=bedid
                    )
                    #Start comparison
                    for comparison in zip(bedinfo, lsNewBed):
                        if list(comparison[0])[1:] != comparison[1]:    
                            booCheck = False
                            # print("FALSE!")
                            
                            break                           #Make sure that new .bed targets and extracted targets are 100% similar
                    
                    if booCheck:                            #booCheck will be false if mismatch was found
                        if len(bedinfo) == len(lsNewBed):   #Make sure that parsed .bed and .bed from database has same length!
                            self.comptoollog.info(
                                "Given BED file was already added to database!"
                            )
                            self.comptoollog.info(
                                f"BED file ID is: {bedid}"
                            )
                            bedid = bedid
                            break                           #Returns found .bed ID to main class instantly    
                        else:
                            booCheck = False
            if not booCheck:       #If no similar .bed files can be found in database
                self.comptoollog.info(
                    "The given BED file was not yet added to the database!"
                )
                
                intID = max(lsAllids) + 1     
                #print(intID)
                bedid = intID
                self.comptoollog.info(f"The BED-ID is: {bedid}")
                self.AddBedToDatabase(bedid)
                #Adding new .bed targets to database 
            return bedid
    
    
    
    def AddBedToDatabase(self, bedid):
        lsADD = []
        for chromosome in self.dictTargets.keys():
            for target in self.dictTargets[chromosome]:
                lsADD.append(
                    [bedid, chromosome, int(target[0]), int(target[1])]
                )
        DBManager.InsertData(lsADD,"Bedfiles", self.dbpath)  
        DBManager.InsertData([bedid, self.dtime], "BedInfo", self.dbpath)
        
        
    ###Function to reduce the .bed file
    def ReduceBedTargets(self):
        newDic ={}
        for chromosomes in self.dictTargets.keys():
            countTarget = len(self.dictTargets[chromosomes])
            countTargetRed = int(round(countTarget * 0.1,0))        #90% reduction of all targets
            lsReducer = [True if x in range(countTargetRed) \
                         else False for x in range(countTarget)
            ]
            
            #Pick targets randomized
            random.shuffle(lsReducer)   
            newDic[chromosomes] = []
            counter = 0
            for boo in lsReducer:
                if boo:
                    newDic[chromosomes].append(
                        self.dictTargets[chromosomes][counter]
                    )
                counter += 1
        return newDic
     

            
    


    ###Checks if Sample is already in database
    def LookAtArgsCompare(self, parsedargs, filetype, dicChecksum, nameTable):
        lsArgnameToDb = []
        lsArgnameComplete = []
        dfNames = self.GetNames(
            nameTable, parsedargs, filetype
        )
        
        for checksum in dicChecksum.keys():
            booAdd, self.ID = self.DbAddAdapter(
                checksum, filetype
            )
            argname = dfNames.loc[
                            dfNames["File"] == \
                            dicChecksum[checksum]
                        ]["Name"].values[0]
            if booAdd:
                lsArgnameToDb.append(argname)
                self.DbInsertNewAdapter(
                    argname, checksum, filetype
                )
                lsArgnameComplete.append(argname)
                self.dicIDs[self.ID] = [
                    argname, dicChecksum[checksum], filetype
                ]
            else:
                self.dicIDs[self.ID] = [
                    argname, dicChecksum[checksum], filetype
                ]
                lsArgnameComplete.append(argname)
            #sys.exit()
        return lsArgnameToDb, lsArgnameComplete
        
        
    def GetNames(self, nameTable, files, filetype):
        if nameTable:
            dfNameTable = pd.read_csv(nameTable)
        else:
            if filetype == "bam":
                values = DBManager.ExtractData(
                    "BamInfo",self.dbpath, bedid = self.bedid, 
                    strReduce = self.Reduce, subsamples = self.subsample, 
                    FileClass = self.FileClass
                )
                intStartCount = len(values)
            elif filetype == "vcf":
                values = DBManager.ExtractData(
                    "VCFInfo", self.dbpath, bedid = self.bedid, 
                    FileClass = self.FileClass
                )
                intStartCount = len(values)
            COUNTER = intStartCount
            lsNames = []
            for file in files:
                if filetype == "bam":
                    lsNames.append([file, f"BAM{COUNTER}"])
                elif filetype == "vcf":
                    lsNames.append([file, f"VCF{COUNTER}"])
                COUNTER += 1
            dfNameTable = pd.DataFrame(lsNames, columns = ["File", "Name"])
            pathCSV = glob(self.outputpath + "NameList.csv")
            try:
                pathCSV[0]
                mode = "a"
            except:
                mode = "w"
            with open(
                    self.outputpath+"NameList.csv", f"{mode}"
                    ) as csvfile:
                writedata = csv.writer(csvfile)
                writedata.writerow(dfNameTable.columns)
                writedata.writerows(dfNameTable.values)
        return dfNameTable
        
        
    def DbAddAdapter(self, checksum, filetype):
        booAdd, intID = DBManager.CheckDB(
            self.bedid, self.dbpath, checksum, filetype, self.Reduce, 
            self.subsample, self.FileClass
        )
        return booAdd, intID


    def DbInsertNewAdapter(self, argname, checksum, filetype):
        if filetype == "bam":
            DBManager.InsertData(
                (
                    argname, self.ID, self.bedid, self.version, self.Reduce, 
                    self.subsample, checksum, self.FileClass, "No",
                    self.dtime
                ), 
                "BamInfo", 
                self.dbpath
            )
        elif filetype == "vcf":
            DBManager.InsertData(
                (
                    argname, self.ID, self.bedid, self.version, checksum, 
                    self.FileClass, "No", self.dtime
                ), 
                "VCFInfo", 
                self.dbpath
            )
                    
                    

       

        