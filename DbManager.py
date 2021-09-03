import sqlite3
import os
import sys
import pandas as pd
import getpass
from glob import glob
import logging


"""
Class for manage every step associated with creation/insertion/extraction/modification of the database
The database consists of following tables:
    - ReadStatistiks
        Includes various information about read per target
            e.g. read on target, mean coverage, g-c content
        You may have a look at the function ExtractQCdata() in the class ExtractInfoData() in the script ExtractFromBam.py for more information
        
    - ReadMapping   
        Includes information about mapped/unmapped reads per chromosome
       
    - QCmetrics
        Includes mean PHRED score per base, the corresponding standard derivation, length distribution and read orientation
            see class ExtractInfoData() in script ExtractFromBam.py
            
    - BamInfo 
        Includes information about bamname, bedid, version of db for each added bamfile
            is needed as checkup if bamfile was already added and which bedfile was associated
            
    - Extracted_Variants
        Includes information about each variant in the passed .vcf files
            e.g. position, reference, variant, allelefrequency
            
    - Sorted_out
        Includes Variants from .vcf files causing problems at allelefrequency calculation
            - variants included in this table are NOT plotted or analysed by the program!
        See function HelpExtract() of class ExtractFromVCF() defined in script ExtractFromVCF.py for more information
    
    - VCFinfo
        Same as BamInfo but for passed .vcf files
    
    - Bedfiles
        Important to calculate bed file ID number. Is used at Preparation (see function CheckBed() in class DataPreparaion() defined in script Preparation.py  
"""
class DBManager():
    def __init__(self, DBpath=False):
        ###Check path to given database and if a database is available
        
        #Initiale logging tool
        self.compylog = logging.getLogger("ComPy")
        
        #Checks database
        self.pathDB, self.booDB = self.FindDatabase(DBpath)
        if not self.booDB:
            self.CreateDB()                         #If no database can be detected by FindDatabase() a new one will be created
        if self.booDB:
            self.DeleteNotFinished()
    
        
    
    
        
    ###Delete all samples flagged as unfinished (Table BamInfo Column Finished == No)
    """
    BAM STRING, ID INT, BedID INTEGER, Version STRING, Reduced INTEGER, 
    Subsamples INTEGER, MD5_CheckSum, FileClass, Finished STRING, 
    Date STRING
    """
    def DeleteNotFinished(self):
        conn = sqlite3.connect(self.pathDB)
        cur = conn.cursor()
        
        ##Delete all .bam files not yet finished
        DataNotFinishedBam = cur.execute(
            "SELECT * FROM BamInfo WHERE Finished == 'No';"
        )
        dataBam = DataNotFinishedBam.fetchall()
        conn.commit()
        DataNotFinishedVCF = cur.execute(
            "SELECT * FROM VCFInfo WHERE Finished == 'No';"
        )
        dataVCF = DataNotFinishedVCF.fetchall()
        conn.commit()
        

        #conn.commit()
        conn.close()
        
        for files in dataBam:
            self.compylog.info(
                f"Sample {files} was not finished yet and will be deleted "
                +"from database"
            )
            DBManager.DelData(
                dbpath = self.pathDB, bamID = files[1]
            )

        for files in dataVCF:
            self.compylog.info(
                f"Sample {files} was not finished yet and will be deleted "
                +"from database"
            )
            DBManager.DelData(
                dbpath = self.pathDB, vcfID = files[1]
            )
        

        
        
        
        
        
    ###Takes given dbpath argument and searches for database. If no database path is provided it will search at default path (current wd + 01_Result/Extraction.db)
    def FindDatabase(self, dbpath):
        if dbpath:
            try:
                if dbpath[0] == ".":
                    dbpath = glob(os.getcwd() + "/" + dbpath)[0]  #If user passed e.g. ./path or ../../path via argparser
                else:
                    dbpath = glob(dbpath)[0]
                booDB = True
            except:
                raise NameError(
                    "Provided database could not be found at given path: "
                    +f"{dbpath} !"
                )
                self.compylog.critical(
                    "Provided database could not be found at given path: "
                    +f"{dbpath} !"
                )
                self.compylog.info("System interrupts...")
                sys.exit()
        else:
            dbpath = f"/home/{getpass.getuser()}/ComPy/.database/Extraction.db"
            try:
                glob(dbpath)[0]
                booDB = True
            except:
                if not os.path.exists(
                        f"/home/{getpass.getuser()}/ComPy/.database/"
                        ):
                    os.makedirs(
                        f"/home/{getpass.getuser()}/ComPy/.database/"
                    )
                booDB = False
        self.compylog.info(
            f"ComPy will work with database {dbpath} \n"
        )
        return dbpath, booDB
    
    
    
    
    
    ###Create the main database if no database can be found
    def CreateDB(self):
        """
        Database tables:
        
            - BamInfo
            - UmiInfo
            - ReadStatistiks    
            - ReadMapping
            - QCmetrics
            - Extracted_Variants
            - Sorted_out
            - VCFInfo
        """
        
        ##Define table headers
        tplReadStat = str(
            "(BAM, ID INT, Chrom, Start, End, Total INTEGER, "
            +"Mapped INTEGER, MeanC INTEGER, GC)"
        )
        tplReadMap = str(
            "(BAM, ID INT, Chrom, Total INTEGER, Mapped INTEGER, "
            +"Unmapped INTEGER)"
        )
        lsQC = str(
            "(BAM, ID INT, Readposition INTEGER, PHREDmean FLOAT, "
            +"PHREDsd FLOAT, ReadCount INTEGER, Mate)"
        )
        lsBamInfo = str(
            "(BAM STRING, ID INT, BedID INTEGER, Version STRING, Reduced INTEGER, "
            +"Subsamples INTEGER, Checksum, FileClass, Finished STRING, "
            +"Date STRING)"
        )
        lsVar = str(
            "(VCF_name, ID INT, Variant_ID, Chromosome, "
            +"Position, Length, Reference, Alternative, Type, Genotype, "
            +"Total_reads, Supporting_reads, Allelefrequency, "
            +"Sample, Info BLOB)"
        )
        lsVCFinfo = str(
            "(VCF_name, ID INT, BedID INTEGER, Version, Checksum, "
            +"FileClass, Finished, Date)"
        )
        lsBED = str(
            "(BedID INTEGER, Chrom, Start INTEGER, End INTEGER)"
        )
        lsBedInfo = str(
            "(BedID STRING, Date STRING)"
        )
        
        ##Create database and tables
        conn = sqlite3.connect(self.pathDB)
        cur = conn.cursor()
        cur.execute(
            f"CREATE TABLE IF NOT EXISTS ReadMapping {tplReadMap};"
        )
        cur.execute(
            f"CREATE TABLE IF NOT EXISTS ReadStatistiks {tplReadStat};"
        )
        cur.execute(
            f"CREATE TABLE IF NOT EXISTS QCmetrics {lsQC};"
        )
        cur.execute(
            f"CREATE TABLE IF NOT EXISTS BamInfo {lsBamInfo};"
        ) 
        cur.execute(
            f"CREATE TABLE IF NOT EXISTS Extracted_Variants {lsVar};"
        )
        cur.execute(
            f"CREATE TABLE IF NOT EXISTS Sorted_out {lsVar};"
        )
        cur.execute(
            f"CREATE TABLE IF NOT EXISTS VCFInfo {lsVCFinfo};"
        )
        cur.execute(
            f"CREATE TABLE IF NOT EXISTS Bedfiles {lsBED};"
        )   
        cur.execute(
            f"CREATE TABLE IF NOT EXISTS BedInfo {lsBedInfo};"
        )                                
        conn.commit()
        conn.close()
        
        
        
        
        
    ###Function to insert calculated data from other functions/classes to the database
    def InsertData(data, table, dbpath):
        """
        Depending on which kind of data is provided or which table is filled different INSERT command are used
        The function is usually called by other functions from other classes in other scripts
        """ 
        conn = sqlite3.connect(dbpath)
        cur = conn.cursor()
        
        if table in ["BamInfo", "VCFInfo", "BedInfo"]:
            tplVal = ','.join('?' for _ in data)
            strInsert = f"INSERT INTO {table} VALUES ({tplVal});"
            cur.execute(strInsert, data)
            conn.commit()
        
        else:
        
            tplVal = ','.join('?' for _ in data[0])
            strInsert = f"INSERT INTO {table} VALUES ({tplVal});"
            cur.executemany(strInsert, data)
            conn.commit()
        conn.close()





    ###Change status of data in database if all data was extracted from .bam file
    def SecurityChangeData(dbpath, intID, version, tool):
        conn = sqlite3.connect(dbpath)
        cur = conn.cursor()
        if tool == "bam":
            cur.execute("UPDATE BamInfo "
                        +"SET Finished='Yes' "
                        +f"WHERE ID == '{intID}' "
                        +f"AND Version == '{version}';"
            )
        elif tool == "vcf":
            cur.execute("UPDATE VCFInfo "
                        +"SET Finished='Yes' "
                        +f"WHERE ID == '{intID}' "
                        +f"AND Version == '{version}';")
        conn.commit()
        conn.close()
        

    
    """    
    ###Function to EXTRACT data from database
    kwargs: strReduce, bedid, boobed, boocheck, subsamples, ID, FileClass, ALL 
    """
    def ExtractData(table, dbpath, **kwargs):
        """
        Used by scripts which needs data from database (Preparation.py, PlotData.py, WriteCSV.py, PlotAll.py)
        Depending on which data is needed, different SELECT commands are executed
        """
        #Initiale logging tool
        compylog = logging.getLogger("ComPy")

        if "boobed" in kwargs.keys():  
            data = DBManager.ExtractBedFileID(table, dbpath)
            return data
        elif table == "Bedfiles" and "ALL" not in kwargs.keys() and "df" not in kwargs.keys():
            data = DBManager.ExtractBedFileTargets(
                table, dbpath, kwargs["bedid"]
            )
            return data
        elif table == "Bedfiles" and "df" in kwargs.keys():
            data = DBManager.ExtractBedFileDataframe(
                table, dbpath, kwargs["bedid"]
            )
            return data
        # elif table in ["BamInfo", "VCFInfo", "BedInfo"] \
            # and "ALL" in kwargs.keys():
        elif "ALL" in kwargs.keys():
            # print("YO")
            data = DBManager.ExtractAllInfoData(table, dbpath)
            return data
        elif table in ["BamInfo", "VCFInfo", "BedInfo"] \
            and "ALL" not in kwargs.keys() \
            and "ID" not in kwargs.keys():
            data = DBManager.ExtractFileTypeAdapter(table, dbpath, kwargs)
            return data
        # elif table not in ["BamInfo", "VCFInfo", "BedInfo"] \
            # and "ALL" in kwargs.keys():
                # data = DBManager.ExtractAllData(table, dbpath)
                # return data
        elif "ID" in kwargs.keys():
            dataFrm = DBManager.ExtractDataAdapter(table, dbpath, kwargs)
            return dataFrm
            
    # def ExtractAllData(table, dbpath):
        # conn = sqlite3.connect(dbpath)
        # dfdata = pd.read_sql_query(
            # f"SELECT * FROM {table};",
            # conn
        # )
        # print("YO")
        # cur = conn.cursor()
        # data = cur.execute(f"SELECT BedID FROM {table};")
        # rows = data.fetchall()
        # conn.close()
        # return dfdata
    
    def ExtractBedFileDataframe(table, dbpath, bedid):
        conn = sqlite3.connect(dbpath)
        data = pd.read_sql_query(
            f"SELECT * FROM {table} WHERE BedID == {bedid};", conn
        )
        # print(data)
        conn.close()
        return data
    
    def ExtractBedFileID(table, dbpath):
        conn = sqlite3.connect(dbpath)
        cur = conn.cursor()
        data = cur.execute(f"SELECT BedID FROM {table};")
        rows = data.fetchall()
        conn.close()
        return rows
    
    def ExtractBedFileTargets(table, dbpath, bedid):
        conn = sqlite3.connect(dbpath)
        cur = conn.cursor()
        data = cur.execute(
            f"SELECT * FROM {table} "
            +f"WHERE BedID == {bedid}"
        )
        rows = data.fetchall()
        conn.close()
        return rows
    
    def ExtractAllInfoData(table, dbpath):
        conn = sqlite3.connect(dbpath)
        cur = conn.cursor()
        data = pd.read_sql_query(f"SELECT * FROM {table};", conn)
        #data = cur.execute(f"SELECT * FROM {table};")
        #rows = data.fetchall()
        #names = [x[0] for x in data.description]
        conn.close()
        return data    
    

    def ExtractFileTypeAdapter(table, dbpath, allargs):
        if table == "BamInfo":
            return DBManager.ExtractTargetInfoDataBAM(
                table, dbpath, allargs["bedid"], allargs["FileClass"], 
                allargs["strReduce"], allargs["subsamples"]
            )
        elif table == "VCFInfo":
            return DBManager.ExtractTargetInfoDataVCF(
                table, dbpath, allargs["bedid"], allargs["FileClass"]
            )

    
    def ExtractTargetInfoDataBAM(table, dbpath, bedid, FileClass, strReduce, 
                                 subsamples):
        conn = sqlite3.connect(dbpath)
        cur = conn.cursor()
        data = cur.execute(
            f"SELECT * FROM {table} "
            +f"WHERE BedID == {bedid} "
            +f"AND Reduced == '{strReduce}' "
            +f"AND Subsamples == {subsamples} "
            +f"AND FileClass == '{FileClass}';"
        )
        rows = data.fetchall()
        conn.close()
        return rows
    
    def ExtractTargetInfoDataVCF(table, dbpath, bedid, FileClass):
        conn = sqlite3.connect(dbpath)
        cur = conn.cursor()
        data = cur.execute(
            f"SELECT * FROM {table} "
            +f"WHERE BedID == {bedid} "
            +f"AND FileClass == '{FileClass}';"
        )
        rows = data.fetchall()
        conn.close()
        return rows
        
        
    def ExtractDataAdapter(table, dbpath, allargs):
        try:
            if allargs["ID"] == int(allargs["ID"]):
                allargs["ID"] = [allargs["ID"]]
        except:
            pass
        if len(allargs["ID"]) > 1:
            kword = "in"
            allargs["ID"] = tuple(allargs["ID"])
        else:
            kword = "=="
            allargs["ID"] = allargs["ID"][0]

        dfdata = DBManager.ExtractDataFrame(table, dbpath, kword, allargs)
        return dfdata
        
    def ExtractDataFrame(table, dbpath, kword, allargs):
        conn = sqlite3.connect(dbpath)
        # try: 
        intID = allargs["ID"]
        dfdata = pd.read_sql_query(
            f"SELECT * FROM {table} "
            +f"WHERE ID {kword} {intID} ",
            conn
       )

        
        
        conn.close()
        return dfdata
        
    
    
        

        
        
        
        
        
            
    

    """
    Function to delete data from the database
    Only called if user wants to delete data
        short info: .bed files cannot be removed from database
    
    Two types of data deletion are supported:
        1) User defined deletion   
            Data will be deleted when fit parameters:
                .bam name
                .bed identifier (BED ID)
                Was .bed file reduced?
                How many subsamples?
        2) Database cleaning
            Data will be deleted accourding to version
            Only data with version not appearing in list of compatible versions (see ComPy.py, section parse the arguments) will be deleted
            --CLEAN or -clean parameter has to be given in the command line!
    kwargs:
        bedid, dbpath, checksum, reduceBed, intSub, booClean, lsCompVers, FileClass
    """
    def DelData(**kwargs):
        
        dbpath = kwargs["dbpath"]
        if "bamID" in kwargs.keys():
            if kwargs["bamID"]:
                intID = kwargs["bamID"]
                if type(intID) == int:
                    intID = [intID]
                DBManager.DelDataBam(dbpath, intID)
            
        if "vcfID" in kwargs.keys():
            if kwargs["vcfID"]:
                intID = kwargs["vcfID"]
                if type(intID) == int:
                    intID = [intID]
                DBManager.DelDataVcf(dbpath, intID) 
        
        if "booClean" in kwargs.keys():
            databam = DBManager.GetDataCleanBam(dbpath, kwargs["lsCompVers"])
            for row in databam: 
                DBManager.DelDataBam(
                    dbpath, row[0], row[1], row[2], row[3], row[4]
                )
            datavcf = DBManager.GetDataCleanVCF(dbpath, kwargs["lsCompVers"])
            for row in datavcf: 
                DBManager.CleanDataDelVCF(
                    dbpath, row[0], row[1], row[2]
                )
            
            
            
            
    def DelDataVcf(dbpath, ID): 
        #Initiale logging tool
        compylog = logging.getLogger("ComPy")
        
        conn = sqlite3.connect(dbpath)
        cur = conn.cursor()
        
        lsTableNamesVCF = ["Extracted_Variants","Sorted_out","VCFInfo"]
        for intID in ID:
            for strTable in lsTableNamesVCF:
                cur.execute(
                    f"DELETE FROM {strTable} "
                    +f"WHERE ID == '{intID}' ;"
                )
                conn.commit()
                compylog.info(
                    f"The variants taken from File {intID} was removed "
                    +f"from table {strTable}! \n"
                )
        conn.commit()
        conn.close()
        
  
            
    def DelDataBam(dbpath, ID):
        #Initiale logging tool
        compylog = logging.getLogger("ComPy")
        
        conn = sqlite3.connect(dbpath)
        cur = conn.cursor()

        lsTableNamesBAM = [
            "ReadStatistiks", "QCmetrics", "BamInfo", "ReadMapping"
        ]
        
        if type(ID) == int:
            ID = [ID]
        
        for intID in ID:
            for strTable in lsTableNamesBAM:
                cur.execute(
                    f"DELETE FROM {strTable} "
                    +f"WHERE ID == '{intID}' ;"
                )
                conn.commit()
                compylog.info(
                    f"The bamfile {intID} was removed from "
                    +f"table {strTable}! \n"
                )
        conn.close()
        
        
        
    def GetDataCleanBam(dbpath, lsCompVers):   
        conn = sqlite3.connect(dbpath)
        cur = conn.cursor()
        if len(lsCompVers) > 1:
            data = cur.execute(
                "SELECT ID, BedID, Reduced, Subsamples, "
                +f"FileClass FROM BamInfo WHERE Version NOT IN {lsCompVers}"
            )
        else:
            versions = lsCompVers[0][0]
            data = cur.execute(
                "SELECT ID, BedID, Reduced, Subsamples, "
                +"FileClass FROM BamInfo WHERE Version IS NOT "
                +f"'{versions}'"
            )
        row = data.fetchall()       
        return row
        
    
    
            
        
    def GetDataCleanVCF(dbpath, lsCompVers):
        conn = sqlite3.connect(dbpath)
        cur = conn.cursor()
        if len(lsCompVers) > 1:
                data = cur.execute(
                    "SELECT ID, BedID, FileClass "
                    +"FROM VCFInfo "
                    +f"WHERE Version NOT IN {lsCompVers}"
                )
        else:
            versions = lsCompVers[0][0]
            data = cur.execute(
                "SELECT ID, BedID, FileClass "
                +"FROM VCFInfo "
                +f"WHERE Version IS NOT '{versions}'"
            )
        row = data.fetchall()
        conn.commit()
        conn.close()
        return row
    
    
    
    
    # ###Check if file was already added to database
    # #Is used by script Preparation.py
    # """
    # DB table header BAM:
    # BAM STRING, ID INT, BedID INTEGER, Version STRING, Reduced INTEGER, 
    # Subsamples INTEGER, MD5_CheckSum, FileClass, Finished STRING, 
    # Date STRING
    
    # DB table header VCF
    # VCF_name, ID INT, BedID INTEGER, Version, MD5_CheckSum, 
    # FileClass, Finished, Date
    # """
    # def CheckDB(bedid, dbpath, checksum, filetype, strReduce, intSubsample,
    #             fileclass):
        
    #     #Initiale logging tool
    #     compylog = logging.getLogger("ComPy")
        
    #     conn = sqlite3.connect(dbpath)
    #     cur = conn.cursor()
    #     booAdd = True
    #     if filetype == "bam":
    #         bamnames = cur.execute("SELECT * FROM {}".format("BamInfo"))
    #         rows = bamnames.fetchall()
    #         for bamfiles in rows:
    #             if bamfiles[2] == bedid \
    #                 and bamfiles[4] == strReduce \
    #                     and bamfiles[5] == intSubsample \
    #                         and bamfiles[6] == checksum \
    #                             and bamfiles[7] == fileclass \
    #                                 and bamfiles[8] == "Yes":
    #                 compylog.info(
    #                     f"The bamfile: '{bamfiles[0]}' was already added to "
    #                     +f"the database! Reduced: {strReduce}, "
    #                     +f"Version: {bamfiles[3]}, Subsamples: {intSubsample}, "
    #                     +f"FileClass: {fileclass}, ID: {bamfiles[1]}, "
    #                     +f"MD5_Sum: {checksum}"
    #                 )
    #                 booAdd = False
    #                 intID = bamfiles[1]
    #                 break
    #     else:
    #         vcfnames = cur.execute("SELECT * FROM {}".format("VCFInfo"))
    #         rows = vcfnames.fetchall()
    #         for vcffiles in rows:
    #             #print(vcffiles)
    #             #print(bedid)
    #             #print(checksum)
    #             #print(fileclass)
    #             if vcffiles[2] == bedid \
    #                 and vcffiles[4] == checksum \
    #                     and vcffiles[5] == fileclass \
    #                         and vcffiles[6] == "Yes":
    #                 compylog.info(
    #                     f"The VCFfile: '{vcffiles[0]}' was already added "
    #                     +f"to the database! ID: {vcffiles[1]}, BedID: {bedid},"
    #                     +f"MD5_Sum: {checksum}, FileClass: {fileclass}"
    #                 )
    #                 booAdd = False
    #                 intID = vcffiles[1]
    #                 #print(checksum)
    #                 #print("HIT!")
    #                 break
    #         #print(booAdd)
    #     if booAdd:
    #         try:
    #             allIds = [x[1] for x in rows]
    #             #print(allIds)
    #             #print("bin ich hier?")
    #             #print(checksum)
    #             #print(booAdd)
    #             #sys.exit()
    #             booIDbetween = False
    #             for number in range(1, len(allIds)):
    #                 # print(f"Number: {number}")
    #                 # print(f"Control: {allIds[number-1]}")
    #                 if number != allIds[number-1]:
    #                     intID = number
    #                     booIDbetween = True
    #                     break
    #             if not booIDbetween:
    #                 intID = len(allIds) + 1
    #         except Exception as e:
    #             compylog.info(
    #                 "No file yet added to database"
    #             )
    #             print(e)
    #             intID = 1
        
    #     #print(intID)
    #     #print(booAdd)
    #     #sys.exit()
    #     return booAdd, intID
        
        
        
        
        
    ###Get .bed ID for .bed file recovery
    def GetBedID(bam, pathDB):
        
        #Initiale logging tool
        compylog = logging.getLogger("ComPy")
        
        conn = sqlite3.connect(pathDB)
        cur = conn.cursor()
        if "/" in bam:
            bam = bam.split("/")[-1]
        if "." in bam:
            bam = bam.split(".")[0]
        data = cur.execute(
            f"SELECT BedID FROM BamInfo WHERE BAM == '{bam}';"
        )
        rows = data.fetchall()
        conn.close()
        if len(rows) > 1:
            print("The .bam file {} was associated with more then one .bed file. Please use 'extract info' for more information".format(bam))
            compylog.critical("The .bam file {} was associated with more then one .bed file. Please use 'extract info' for more information".format(bam))
            compylog.info("System interrupts!")
            sys.exit()
        if len(rows) < 1:
            print("No .bam file was found in database with this name ({}). Please use 'extract info' for more information".format(bam))
            compylog.info("No .bam file was found in database with this name ({}). Please use 'extract info' for more information".format(bam))
            compylog.info("System interrupts!")
            sys.exit()
        return rows
        
        
        
        
    