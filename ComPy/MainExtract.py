###Import needed packages
#External packages
import os
import logging 
import shutil

#Own scripts
from .DbManager import DBManager
from .Preparation import DataPreparation
from .WriteCSV import WriteCSV

"""
###Class excepts following kwargs:
argBamfiles, argVcffiles, argBedfile, argDatabase, argOutput, argSubsample, 
argReduce, dtime, strVersion
"""
def CompToolExtract(tool, **kwargs):
    
    ###Prepare logging 
    compylog = logging.getLogger("ComPy")
    

    
    ###If user wants to extract data from database as .xlsx
    #KWARGS: argDatabase, argReduce, argChecksum, argSubsample, argClass, dtime, argOutput, argNameList
    if tool=="data":

                    
        """
        Workflow
            1) Prepare database
            2) Prepare arguments
            3) Extract needed data from database
        """
        
        """
        1) Prepare database
        """
        classDataBase = DBManager(kwargs["argDatabase"])                       #Class defined in script DBManager.py        
        
        
        """
        2) Prepare arguments
        """
        #Class defined in Preparation.py
        classPrep = DataPreparation(
            "extract", booDB = True, DBpath = classDataBase.pathDB, 
            dtime = kwargs["dtime"], version = kwargs["strVersion"], 
            nameTable = kwargs["argNameList"], 
            lsbamid = kwargs["argbamid"], lsvcfid = kwargs["argvcfid"]
        )
        
        
        """
        3) Extract needed data from database
        """
        compylog.info("Converting database to .xlsx")
        
        #Class is defined in WriteCSV.py
        clSaveData = WriteCSV(
            kwargs["argOutput"], classDataBase.pathDB, kwargs["dtime"], 
            bamID = classPrep.dicIDs["bam"], 
            vcfID = classPrep.dicIDs["vcf"]
        )   
        # print(classPrep.dicIDs)
        if classPrep.dicIDs["bam"]:
            # print("HI")
            clSaveData.WriteData("ReadStatistiks")
            clSaveData.WriteData("QCmetrics")
            clSaveData.WriteData("ReadMapping")
        
        else:
            compylog.info("No bam data was asked for")
            
        if classPrep.dicIDs["vcf"]:
            # print("HO")
            clSaveData.WriteData("Extracted_Variants")
            clSaveData.WriteData("Sorted_out")
        else:
            compylog.info("No VCF data was provided")
            
    
    

    
    ###If user wants to recover .bed file from database
    #KWARGS: argBED, argOutput, argDatabase
    if tool=="bedrec":
        """
        Workflow
            1) Prepare database and get .bed id if not provided
            2) Recover .bed file targets from database
        """
    
        """
        1) Prepare database and get .bed id if not provided
        """
        classDataBase = DBManager(kwargs["argDatabase"])              #Class defined in script DBManager.py
        
        
        
        """
        2) Recover .bed file targets from database
        """
        intBedID = kwargs["argBedfile"]
        compylog.info(
            "Recover .bed file targets of BedID : {}".format(intBedID)
        )
        clSaveData = WriteCSV(
            kwargs["argOutput"], classDataBase.pathDB, kwargs["dtime"], 
            bedid = intBedID
        )
        clSaveData.RecoverBED()
        
        
        
       
        
        
        
    ###If user wants to extract info about data in database as .xlsx
    if tool=="info":
        """
        Workflow
            1) Prepare database
            2) Extract info data from database
        """
        
        """
        1) Prepare database
        """
        classDataBase = DBManager(kwargs["argDatabase"])              #Class defined in script DBManager.py          
        
        
        """
        2) Extract info data from database
        """
        compylog.info("Converting database to .xlsx")
        
        #Class is defined in WriteCSV.py
        clSaveData = WriteCSV(
            kwargs["argOutput"], classDataBase.pathDB, kwargs["dtime"]
        )   
        clSaveData.WriteData("BamInfo")
        clSaveData.WriteData("VCFInfo")
        clSaveData.WriteData("BedInfo")
        compylog.info("Saved info files!")
        
        
        
        
        
        
        

    ###If user wants to export database
    if tool == "database":
        compylog.info("The database will be exported!")
        #Find database
        classDataBase = DBManager(kwargs["argDatabase"])          #Class defined in script DBManager.py
        
        #Make sure to use right output
        if kwargs["argOutput"][:2] == ".." or kwargs["argOutput"][:2] == "./":
            strOutput = os.getcwd() + "/" + kwargs["argOutput"]
            if strOutput[-1] != "/":
                strOutput = strOutput + "/"
        else:
            strOutput = kwargs["argOutput"]
            if strOutput[-1] != "/":
                strOutput = strOutput + "/"  
        
        #Copy database
        if not os.path.exists(strOutput):
            os.makedirs(strOutput)
        dtime = kwargs["dtime"]
        pathout = strOutput + f"{dtime}_ExportedDatabase.db"
        # print(pathout)
        # print(classDataBase.pathDB)
        shutil.copyfile(
            classDataBase.pathDB, pathout
        )
        
        compylog.info(
            f"The database was exported from {classDataBase.pathDB} to "
            + f"{pathout}"
        )
        