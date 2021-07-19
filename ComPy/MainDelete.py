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
    comptoollog = logging.getLogger("ComparisonTool")
    
    
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
            "delete", intID = kwargs["argIntID"], 
            Datatype = kwargs["argDatatype"],  booDB = classDataBase.booDB, 
            DBpath= classDataBase.pathDB, nameTable = kwargs["argNameTable"], 
            dtime= kwargs["dtime"]
        )
   
        """
        3) Delete data
        """
        #Funktion defined in script DBManager.py
        DBManager.DelData(
            dbpath = classDataBase.pathDB, vcf = classPrep.vcf, 
            ID = classPrep.ID
        )