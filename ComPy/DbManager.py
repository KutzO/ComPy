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
        
    - ReadMapping   
        Includes information about mapped/unmapped reads per chromosome
       
    - QCmetrics
        Includes mean PHRED score per base, the corresponding standard derivation, length distribution and read orientation
            
    - BamInfo 
        Includes information about e.g bamname, bedid, version of db for each added bamfile
            
    - Extracted_Variants
        Includes information about each variant in provided .vcf files
            e.g. position, reference, variant, allelefrequency
    
    - VCFinfo
        Same as BamInfo but for .vcf files
    
    - Bedfiles
        Stores bed file targets. 
"""
class DBManager():
    def __init__(self, DBpath=False):
    
        #Initiale logging tool
        self.compylog = logging.getLogger("ComPy")
        
        #Check database path
        self.pathDB, self.booDB = self.FindDatabase(DBpath)
        if not self.booDB:
            self.compylog.info("A new database will be created!")
            self.CreateDB()                         
        if self.booDB:
            self.DeleteNotFinished()
    
        
    
    
        
    ###Delete all samples flagged as unfinished (Table BamInfo Column Finished == No)
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
        

        
        
        
        
        
    ###Function to identify path to database
    def FindDatabase(self, dbpath):
        if dbpath:
            try:
                if dbpath[0] == ".":
                    dbpath = glob(os.getcwd() + "/" + dbpath)[0]  
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
    
    
    
    
    
    ###Function to create the main database if no database can be found
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
            +"Subsamples INTEGER, Checksum, FileClass, Flag, Finished STRING, "
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





    ###Change status of data in database if all data were extracted from .bam file
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
    
        #Initiale logging tool
        compylog = logging.getLogger("ComPy")

        if "boobed" in kwargs.keys():  
            compylog.info("Extract bedids")
            data = DBManager.ExtractBedFileID(table, dbpath)
            return data
            
        elif table == "Bedfiles" and "ALL" not in kwargs.keys() and "df" not in kwargs.keys():
            compylog.info("Extract bed file targets")
            data = DBManager.ExtractBedFileTargets(
                table, dbpath, kwargs["bedid"]
            )
            return data
            
        elif table == "Bedfiles" and "df" in kwargs.keys():
            compylog.info("Extract bedfile targets as dataframe")
            data = DBManager.ExtractBedFileDataframe(
                table, dbpath, kwargs["bedid"]
            )
            return data

        elif "ALL" in kwargs.keys():
            compylog.info("Extract all info data")
            data = DBManager.ExtractAllInfoData(table, dbpath)
            return data
            
        elif table in ["BamInfo", "VCFInfo", "BedInfo"] \
            and "ALL" not in kwargs.keys() \
            and "ID" not in kwargs.keys():
            compylog.info("Extract specific info data")
            data = DBManager.ExtractFileTypeAdapter(table, dbpath, kwargs)
            return data

        elif "ID" in kwargs.keys():
            compylog.info("Extract data showing specific ID")
            dataFrm = DBManager.ExtractDataAdapter(table, dbpath, kwargs)
            return dataFrm
            

    ###Helper function to extract bed targets as dataframe
    def ExtractBedFileDataframe(table, dbpath, bedid):
        conn = sqlite3.connect(dbpath)
        data = pd.read_sql_query(
            f"SELECT * FROM {table} WHERE BedID == {bedid};", conn
        )
        conn.close()
        return data
    
    ###Helper function to extract bed IDs
    def ExtractBedFileID(table, dbpath):
        conn = sqlite3.connect(dbpath)
        cur = conn.cursor()
        data = cur.execute(f"SELECT BedID FROM {table};")
        rows = data.fetchall()
        conn.close()
        return rows
    
    ###Helper function to extract bedfile targets
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
    
    ###Helper function to extract info data as dataframe
    def ExtractAllInfoData(table, dbpath):
        conn = sqlite3.connect(dbpath)
        cur = conn.cursor()
        data = pd.read_sql_query(f"SELECT * FROM {table};", conn)
        conn.close()
        return data    
    
    ###Adapter function if bam or vcf info is needed
    def ExtractFileTypeAdapter(table, dbpath, allargs):
        if table == "BamInfo":
            return DBManager.ExtractTargetInfoDataBAM(
                table, dbpath, allargs["bedid"], allargs["FileClass"], 
                allargs["strReduce"], allargs["subsamples"], allargs["flag"]
            )
        elif table == "VCFInfo":
            return DBManager.ExtractTargetInfoDataVCF(
                table, dbpath, allargs["bedid"], allargs["FileClass"]
            )

    ###Helper function to extract specific bam data
    def ExtractTargetInfoDataBAM(table, dbpath, bedid, FileClass, strReduce, 
                                 subsamples, flag):
        if len(FileClass) > 1 \
         and type(FileClass) == list:
            FileClass = tuple(FileClass)
            key = f"in {FileClass}"
        else:
            if type(FileClass) == list:
                FileClass = FileClass[0]
            key = f"== '{FileClass}'"
        conn = sqlite3.connect(dbpath)
        cur = conn.cursor()
        data = cur.execute(
            f"SELECT * FROM {table} "
            +f"WHERE BedID == {bedid} "
            +f"AND Reduced == '{strReduce}' "
            +f"AND Subsamples == {subsamples} "
            +f"AND Flag == {flag} "
            +f"AND FileClass {key};"
        )
        rows = data.fetchall()
        conn.close()
        return rows
    
    ###Helper function to extract specific vcf data
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
    
    ###Adapter function to determine if one or more data should be extracted
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
    
    ###Helper function to extract data in dataframe
    def ExtractDataFrame(table, dbpath, kword, allargs):
        conn = sqlite3.connect(dbpath) 
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
            for intID in databam: 
                DBManager.DelDataBam(
                    dbpath, intID
                )
            datavcf = DBManager.GetDataCleanVCF(dbpath, kwargs["lsCompVers"])
            for intID in datavcf: 
                DBManager.DelDataVcf(
                    dbpath, intID
                )
            
            
            
    ###Helper function to delete vcf data
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
        
  
    ###Helper function to delete bam data
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
        
        
    ###Helper function to identify outdated bam data
    def GetDataCleanBam(dbpath, lsCompVers):   
        conn = sqlite3.connect(dbpath)
        cur = conn.cursor()
        if len(lsCompVers) > 1:
            data = cur.execute(
                f"SELECT ID FROM BamInfo WHERE Version NOT IN {lsCompVers}"
            )
        else:
            versions = lsCompVers[0][0]
            data = cur.execute(
                "SELECT ID FROM BamInfo WHERE Version IS NOT "
                +f"'{versions}'"
            )
        row = data.fetchall()       
        return row
        
            
    ###Helper function to identify outdated vcf data
    def GetDataCleanVCF(dbpath, lsCompVers):
        conn = sqlite3.connect(dbpath)
        cur = conn.cursor()
        if len(lsCompVers) > 1:
                data = cur.execute(
                    "SELECT ID FROM VCFInfo "
                    +f"WHERE Version NOT IN {lsCompVers}"
                )
        else:
            versions = lsCompVers[0][0]
            data = cur.execute(
                "SELECT ID FROM VCFInfo "
                +f"WHERE Version IS NOT '{versions}'"
            )
        row = data.fetchall()
        conn.commit()
        conn.close()
        return row
    
    