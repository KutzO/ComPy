import random
import sys
import pandas as pd
import os
from .DbManager import DBManager
import logging
import csv
from glob import glob
import getpass
import yaml


###Class to prepare arguments for each main tool
class DataPreparation:
    def __init__(self, tool, bam=False, DBpath=False, booDB=False, bed=False, 
                 subsample=False, dtime=False, outputpath = False, vcf = False, 
                 ReduceBed = False, version = False, checksum = False, 
                 nameTable = False, vcfPlot = False, FileClass = False,
                 intID = False, Datatype = False, lsbamid = False,
                 lsvcfid = False, styleyaml = False, flag = False):
        
        #Define global variables
        self.dbpath = DBpath
        self.booDB = booDB
        self.subsample = subsample
        self.version = version
        self.dtime = dtime
        self.flag = flag
        if FileClass:
            self.FileClass = FileClass
        else:
            self.FileClass = "Default"
        
        #Initiale logging tool
        self.compylog = logging.getLogger("ComPy")
        
        #Prepare arguments depending on which tool is needed
        if tool == "compare":
            self.outputpathTmp, self.outputpath = self.FindOutputPath(
                                                        outputpath
            )
            self.PrepareCompare(
                bed, ReduceBed, bam, vcf, checksum, nameTable, styleyaml
            )
            
        elif tool == "extract" or tool == "delete":
            self.PrepareExtractionDeletion(
                nameTable, lsbamid, lsvcfid
            )
            
        elif tool == "plot":
            self.outputpathTmp, self.outputpath = self.FindOutputPath(
                                                        outputpath
            )
            self.PreparePlot(
                lsbamid, lsvcfid, nameTable, styleyaml
            )
            
    
    
    
    
    ###Helper function to prepare args for tool plot    
    def PreparePlot(self, lsbamid, lsvcfid, nameTable, argyamlpath):
        
        #Prepare style dict for plotting
        self.PrepareStyleYaml(argyamlpath)
        
        ##Prepare id dataframe
        if nameTable:
            dfData = pd.read_csv(nameTable)
        else:
            dfData = pd.DataFrame()
            lsIDs = []
            lsType = []
            lsNames = []
            if lsbamid:
                lsIDs = lsIDs + lsbamid
                lsType = lsType + ["bam" for x in lsbamid]
                lsNames = lsNames + [f"Bam_{x}" for x in lsbamid]
            if lsvcfid:
                lsIDs = lsIDs + lsvcfid
                lsType = lsType + ["vcf" for x in lsvcfid]
                lsNames = lsNames + [f"Vcf_{x}" for x in lsvcfid]
            dfData["ID"] = lsIDs
            dfData["Type"] = lsType
            dfData["Name"] = lsNames
        
        #Save namelist
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
            if mode == "w":
                writedata.writerow(dfData.columns)
            writedata.writerows(dfData.values)
        
        #Get ID dict
        self.dicIDs = {}
        if "bam" in dfData["Type"].values:
            dfSliceData = dfData.loc[dfData["Type"] == "bam"]
            self.FillDictionary(dfSliceData, "bam")
        if "vcf" in dfData["Type"].values:
            dfSliceData = dfData.loc[dfData["Type"] == "vcf"]
            self.FillDictionary(dfSliceData, "vcf")
    
    
    
    ###Helper tool to set up plot dict
    def FillDictionary(self, data, filetype):
        self.dicIDs[filetype] = {}
        for entry in zip(
                data["ID"], 
                data["Name"],
                data["Type"]
                ):
            self.dicIDs[filetype][int(entry[0])] = [
                                    entry[1], entry[2]
                                    ]
        
        
        
        
    ###Helper function to prepare deletion
    def PrepareExtractionDeletion(self, nameTable, bamids, vcfids):
        
        #Define IDs
        if nameTable:
            df = pd.read_csv(nameTable)
            self.dicIDs = {}
            self.dicIDs["bam"] = list(
                df.loc[df["Type"] == "bam"]["ID"]
            )
            self.dicIDs["vcf"] = list(
                df.loc[df["Type"] == "vcf"]["ID"]
            )
            
        else:
            self.dicIDs = {}
            self.dicIDs["bam"] = bamids
            self.dicIDs["vcf"] = vcfids


    
    ###Helper function to prepare comparison 
    def PrepareCompare(self, bed, ReduceBed, bam, vcf, checksum, nameTable,
                       argyamlpath):
        
        #Prepare style dict for plotting
        self.PrepareStyleYaml(argyamlpath)
        
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




    ###Helper tool to mimic style.yaml if not provided
    def PrepareStyleYaml(self, argyamlpath):
        if argyamlpath:
            with open(argyamlpath, "r") as streamfile:
                self.styleyaml = yaml.load(
                    streamfile, Loader = yaml.FullLoader
                )
        else:
            self.styleyaml = {}
            self.styleyaml["plotstyle"] = "seaborn" 
            
            self.styleyaml["compare"] = {}
            self.styleyaml["compare"]["headersize"] = 20
            self.styleyaml["compare"]["xlabel"] = 10
            self.styleyaml["compare"]["ylabel"] = 10
            self.styleyaml["compare"]["axlabel"] = 10
            self.styleyaml["compare"]["title"] = 16
            self.styleyaml["compare"]["font"] = 10
            self.styleyaml["compare"]["titlefont"] = 10
            
            self.styleyaml["database"] = {}
            self.styleyaml["database"]["headersize"] = 20
            self.styleyaml["database"]["xlabel"] = 10
            self.styleyaml["database"]["ylabel"] = 10
            self.styleyaml["database"]["axlabel"] = 10
            self.styleyaml["database"]["title"] = 16
            self.styleyaml["database"]["font"] = 10
            self.styleyaml["database"]["titlefont"] = 10
            
            self.styleyaml["boxplots"] = {}
            self.styleyaml["boxplots"]["colorcode"] = ["#808080"]
            self.styleyaml["boxplots"]["markercolor"] = "orange"
            self.styleyaml["boxplots"]["markerstyle"] = "x"
            self.styleyaml["boxplots"]["coverageyscale"] = "log"
            self.styleyaml["boxplots"]["ymax"] = "max"
            
            self.styleyaml["pieplot"] = {}
            self.styleyaml["pieplot"]["highest"] = 300
            self.styleyaml["pieplot"]["middle"] = 100
            self.styleyaml["pieplot"]["lowest"] = 20
            self.styleyaml["pieplot"]["colorcode"] =  [
                "#55a868", "#4c72b0", "#dd8452", "#c44e52"
            ]
            
            self.styleyaml["qcplot"] = {}
            self.styleyaml["qcplot"]["phredscore"] = "blue"
            self.styleyaml["qcplot"]["lengthdist"] = "orange"
            self.styleyaml["qcplot"]["stderivation"] = "lightblue"
            
            self.styleyaml["barplots"] = {}
            self.styleyaml["barplots"]["vardistyscale"] = "log"
            self.styleyaml["barplots"]["colorcode"] = ["#808080"]
            self.styleyaml["barplots"]["highest"] = 50
            self.styleyaml["barplots"]["middle"] = 30
            self.styleyaml["barplots"]["lower"] = 10
            self.styleyaml["barplots"]["lowest"] = 1
            self.styleyaml["barplots"]["colorthresh"] =  [
                "#55a868", "#4c72b0", "#dd8452", "#c44e52"
            ]
            self.styleyaml["barplots"]["ymax"] = "max"
            
            self.styleyaml["vennplots"] = {}
            self.styleyaml["vennplots"]["color1"] = "green"
            self.styleyaml["vennplots"]["color2"] = "blue"
            self.styleyaml["vennplots"]["color3"] = "red"
            self.styleyaml["vennplots"]["color4"] = "brown"
            self.styleyaml["vennplots"]["color5"] = "purple"
            
            self.styleyaml["tables"] = {}
            self.styleyaml["tables"]["colorcode"] = "#C0C0C0"
    
    
    
    
    
    
    ###Helper function to define outputpath and create temp folders
    def FindOutputPath(self, outputpath):
        
        #Define path
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
            UserID = getpass.getuser()
            strOutputTmp = str(
                f"/home/{UserID}/ComPy/tmp/{self.dtime}/"
            )
            strOutput = str(
                f"/home/{UserID}/ComPy/Results/Result_{self.dtime}/"
            )
        
        #Create folders
        self.compylog.info(f"Create output folder (tmp): {strOutputTmp}")
        if not os.path.exists(strOutputTmp):
            os.makedirs(strOutputTmp)
        
        self.compylog.info(f"Create output folder (results): {strOutput}")
        if not os.path.exists(strOutput):
            os.makedirs(strOutput)
        
        strOutputBigComp = strOutput+"BigCompare/"
        self.compylog.info("Create output folder "
                              +f"(DBcomp): {strOutputBigComp}"
        )
        if not os.path.exists(strOutputBigComp):
            os.makedirs(strOutputBigComp)
        
        return strOutputTmp, strOutput





    ###Helper tool to extract targetdata from BED file
    def SliceBedFile(self, bed):
        
        #Import bedfile
        bedimport = pd.read_csv(bed, sep="\t", header=None, low_memory=False)
        start = 0
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
                    self.compylog.info(
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
        
    

    
        
    ###Helper function to check if bedfile is present
    def CheckBed(self):
        
        #Start at 1 if a new database was created
        if not self.booDB:      
            bedid = 1
            self.AddBedToDatabase(bedid)
            return bedid
        
        else:
            
            #Extract all bed ids
            allbedids = DBManager.ExtractData(
                "Bedfiles", self.dbpath, boobed=True
            )  
            lsAllids = list(set([x[0] for x in allbedids]))
            
            ##Convert current .bed dictionary to list for comparison
            lsNewBed = []   
            for chromosome in self.dictTargets.keys():
                for targets in self.dictTargets[chromosome]:
                    lsNewBed.append(
                        [chromosome,int(targets[0]), int(targets[1])]
                    )
            #Iterate through all .bed ID found in database
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
                            break                           
                    
                    #If bedfile is present
                    if booCheck:                            
                        if len(bedinfo) == len(lsNewBed):   
                            self.compylog.info(
                                "Given BED file was already added to database!"
                            )
                            self.compylog.info(
                                f"BED file ID is: {bedid}"
                            )
                            bedid = bedid
                            break                             
                        else:
                            booCheck = False
            
            #If bedfile is not yet present in database
            if not booCheck:       
                self.compylog.info(
                    "The given BED file was not yet added to the database!"
                )
                intID = max(lsAllids) + 1     
                bedid = intID
                self.compylog.info(f"The BED-ID is: {bedid}")
                self.AddBedToDatabase(bedid)
            return bedid
    
    
    ###Helper tool to add new bedfiles to database
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
            countTargetRed = int(round(countTarget * 0.1,0))        
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
     

            
    


    ###Helper function to check if Sample is already in database
    def LookAtArgsCompare(self, parsedargs, filetype, dicChecksum, nameTable):
        lsArgnameToDb = []
        lsArgnameComplete = []
        
        #Get data names
        dfNames = self.GetNames(
            nameTable, parsedargs, filetype, dicChecksum
        )
        
        #Prepare dic containing file IDs
        self.dicIDs[filetype] = {}
        for checksum in dicChecksum.keys():
            if dicChecksum[checksum][-1] == filetype:
                
                #Check if file is present in database
                booAdd, self.ID = self.CheckDB(
                    self.bedid, self.dbpath, checksum, filetype, self.Reduce, 
                    self.subsample, self.FileClass
                )
                
                #Define file name
                if "ID" in dfNames.columns:
                    argname = dfNames.loc[
                            dfNames["ID"] == self.ID
                        ]["Name"].values[0]
                else:
                    try:
                        argname = dfNames.loc[
                                dfNames["File"] == \
                                    dicChecksum[checksum][0]
                            ]["Name"].values[0]
                    except:
                        raise SystemError("The NameList you provided causes trouble! Please check pathes / names!")
                self.compylog.info(
                    f"File: {argname}, checksum: {checksum}, ID: {self.ID}, "
                    + f"Add: {booAdd}"
                )
                
                #Add info to database if file is not yet present
                if booAdd:
                    lsArgnameToDb.append(argname)
                    self.DbInsertNewAdapter(
                        argname, checksum, filetype
                    )
                
                self.dicIDs[filetype][self.ID] = [
                    argname, dicChecksum[checksum][0], filetype
                ]
                lsArgnameComplete.append(argname)
        return lsArgnameToDb, lsArgnameComplete
        
    
    
    
    ###Helper tool to check if file is present in database
    def CheckDB(self, bedid, dbpath, checksum, filetype, strReduce, 
                intSubsample, fileclass):
        
        #Initiale logging tool
        compylog = logging.getLogger("ComPy")

        #Check bam files
        booAdd = True
        if filetype == "bam":
            data = DBManager.ExtractData(
                "BamInfo", dbpath, ALL = True
            )
            sliceddata = data.loc[
                (data["BedID"] == bedid) 
                & (data["Reduced"] == strReduce)
                & (data["Subsamples"] == intSubsample)
                & (data["Checksum"] == checksum)
                & (data["FileClass"] == fileclass)
                & (data["Flag"] == self.flag)
                & (data["Finished"] == "Yes")
            ]
            if len(sliceddata.values) > 0:
                name = sliceddata["BAM"].values[0]
                version = sliceddata["Version"]
                intID = sliceddata["ID"].values[0]
                compylog.info(
                    f"The bamfile: '{name}' was already added to "
                    +f"the database! Reduced: {strReduce}, "
                    +f"Version: {version}, Subsamples: {intSubsample}, "
                    +f"FileClass: {fileclass}, ID: {intID}, "
                    +f"MD5_Sum: {checksum}, Flag: {self.flag}"
                )
                booAdd = False

        #Check vcf files
        else:
            data = DBManager.ExtractData(
                "VCFInfo", dbpath, ALL = True
            )
            sliceddata = data.loc[
                (data["BedID"] == bedid)
                & (data["Checksum"] == checksum)
                & (data["FileClass"] == fileclass)
                & (data["Finished"] == "Yes")
            ]
            if len(sliceddata.values) > 0:
                name = sliceddata["VCF_name"].values[0]
                intID = sliceddata["ID"].values[0]
                compylog.info(
                    f"The VCFfile: '{name}' was already added "
                    + f"to the database! ID: {intID}, BedID: {bedid},"
                    + f"MD5_Sum: {checksum}, FileClass: {fileclass}"
                )
                booAdd = False
        
        #Define new sample ID
        if booAdd:
            try:
                allIds = [x for x in data["ID"].values]
                allIds.sort()
                booIDbetween = False
                for number in range(1, len(allIds)):
                    if number != allIds[number - 1]:
                        intID = number
                        booIDbetween = True
                        break        
                if not booIDbetween:
                    intID = len(allIds) + 1
            except Exception as e:
                compylog.warning(
                    f"No file yet added to database. Error : {e}"
                )
                intID = 1
        return booAdd, intID
    
    
    
    
    
    
    
    
    
 
    ###Helper tool to get file names
    def GetNames(self, nameTable, files, filetype, dicCheckSum):
        
        #If nametable is provided by user
        if nameTable:
            dfNameTable = pd.read_csv(nameTable)
            if len(dfNameTable.columns) < 3:
                dfNameTable = pd.read_csv(nameTable, sep=";")   
            if filetype in dfNameTable["Type"].values:
                dfNameTable = dfNameTable.loc[dfNameTable["Type"] == filetype]
                if "ID" in dfNameTable.columns:
                    if len(dfNameTable["ID"].values) < len(files):
                        Errormsg = str(
                            f"The number of file IDs in {nameTable} is "
                            + "not equal to the number of files provided by "
                            + "the user!"
                        )
                        raise NameError(Errormsg)
                        self.compylog.error(Errormsg)
                        sys.exit()
                elif "File" in dfNameTable.columns:
                    if len(dfNameTable["File"].values) < len(files):
                        Errormsg = str(
                            f"The number of file names in {nameTable} is "
                            + "not equal to the number of files provided by "
                            + "the user!"
                        )
                        raise NameError(Errormsg)
                        self.compylog.error(Errormsg)
                        sys.exit()
            else:
                nameTable = False
        
        #Create new nametable if nothing was provided
        if not nameTable:
            self.compylog.info("A new nametable will be created")
            dfNameTable = self.CreateNewNameTable(files, filetype, dicCheckSum)
        
        #Save used nametable
        self.SaveNameTable(dfNameTable)
        return dfNameTable
    
            

        
    
    ###Helper tool to create new nametable
    def CreateNewNameTable(self, files, filetype, dicCheckSum):
        if filetype == "bam":
            values = DBManager.ExtractData(
                "BamInfo",self.dbpath, bedid = self.bedid, 
                strReduce = self.Reduce, subsamples = self.subsample, 
                FileClass = self.FileClass, flag = self.flag
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
            booFound = False
            if filetype == "bam":
                for savedfile in values:
                    if savedfile[6] in dicCheckSum.keys():
                        if dicCheckSum[savedfile[6]] == file:
                            name = savedfile[0]
                            booFound = True
                            break   
            
            elif filetype == "vcf":
                for savedfile in values:
                    if savedfile[4] in dicCheckSum.keys():
                        if dicCheckSum[savedfile[4]] == file:
                            name = savedfile[0]
                            booFound = True
                            break
            
            if not booFound:
               name = f"{filetype}{COUNTER}"
               COUNTER += 1
            
            lsNames.append([file, name])
            dfNameTable = pd.DataFrame(
                lsNames, columns = ["File", "Name"]
            )
        return dfNameTable
    
    
    
    
    ###Helper tool to save nametable
    def SaveNameTable(self, dfNameTable):
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
            if mode == "w":
                writedata.writerow(dfNameTable.columns)
            writedata.writerows(dfNameTable.values)
        return dfNameTable
        




    ###Helper tool to insert new file info to info table
    def DbInsertNewAdapter(self, argname, checksum, filetype):
        if filetype == "bam":
            DBManager.InsertData(
                (
                    argname, self.ID, self.bedid, self.version, self.Reduce, 
                    self.subsample, checksum, self.FileClass, self.flag, "No",
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
                    
                    

       

        