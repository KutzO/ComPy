from .DbManager import DBManager
import xlsxwriter
import getpass
import os
import logging



###Extract needed data from database in .xlsx format (Identifies .bam names and .bed ID)
#Only called by ComPy.py if user parses argument -xlsx or --Excel
class WriteCSV():
    def __init__(self, outputpath, pathDB, dtime, bamID = False, 
                 vcfID = False, fileclass = False, bedid = False):
                 
        #Define global variables
        self.pathDB = pathDB                                     
        self.out = self.FindOutputPath(outputpath, dtime)        
        self.bamID = bamID
        self.vcfID = vcfID
        self.fileclass = fileclass
        self.bedid = bedid
        
        #Initiale logging tool
        self.compylog = logging.getLogger("ComPy")


    
    ###Helper function to define result output
    def FindOutputPath(self, givenpath, dtime):
        if givenpath:
            if givenpath[:2] == ".." or givenpath[:2] == "./":
                strResultOut = os.getcwd()+"/"+givenpath
                if strResultOut[-1] != "/":
                    strResultOut = strResultOut + "/"
            else:
                strResultOut = givenpath
                if strResultOut[-1] != "/":
                    strResultOut = strResultOut + "/"  
            strResultOut = strResultOut + f"/Extracted/{dtime}/"
        else:
            strResultOut = str(
                f"/home/{getpass.getuser()}/ComPy/Extracted/{dtime}/"
            )
        
        #Create folders
        if not os.path.exists(strResultOut):
            os.makedirs(strResultOut)
        return strResultOut



    
    ###Main function for converting database to .xlsx
    def WriteData(self, keytable): 
        
        #Collect data and table header
        if keytable in ["BamInfo", "VCFInfo", "BedInfo"]:
            data = DBManager.ExtractData(
                keytable, self.pathDB, ID = False, ALL = True
            )
            names = data.columns
            data = data.values
            
        elif keytable in ["ReadStatistiks", "QCmetrics", "ReadMapping"]:
            data = DBManager.ExtractData(
                keytable, self.pathDB, ID = self.bamID, 
            )
            names = data.columns
            data = data.values
            
        elif keytable in ["Extracted_Variants", "Sorted_out"]:
            data = DBManager.ExtractData(
                keytable, self.pathDB, ID = self.vcfID, 
            )
            names = data.columns
            data = data.values
        
        #If needed data can not be found in database
        if len(data) == 0:
            if keytable == "VCFInfo" or keytable == "BamInfo":
                self.compylog.warning(
                    f"EXCEL ERROR: No data was added to table {keytable} yet!"
                )
                return
            self.compylog.warning(
                "Excel ERROR: The given input files was not found in database "
                +f"table {keytable}!"
            )
            self.compylog.info(
                f"File IDs: BAM: {self.bamID}; VCF: {self.vcfID}"
            )
            return
        
        ##Create .xlsx file
        workbook = xlsxwriter.Workbook(self.out + f"{keytable}.xlsx")
        worksheet = workbook.add_worksheet(keytable)
        
        #Write table header
        for headernum, headerrow in enumerate(names):
            worksheet.write_string(0, headernum, headerrow)
            
        #Write data
        for row_num, row in enumerate(data):
            cCOL = 0
            for col in row:
                worksheet.write_string(row_num + 1, cCOL, str(col))
                cCOL += 1
        workbook.close() 
        self.compylog.info(
            f"Excel sheet {self.out+keytable}.xlsx was saved!"
        )
        
    
    ###Function to write new .bed file (saved as tab separated .bed file)
    def RecoverBED(self):
        
        #Get bedfile targets
        targets = DBManager.ExtractData(
            "Bedfiles", self.pathDB, bedid = self.bedid
        )
        
        #Save bedfile as .bed
        with open(
                self.out + f"RecoveredBedID_{self.bedid}.bed","w"
                  ) as bedfile:
            for target in targets:
                for value in target[1:-1]:
                    bedfile.write(str(value))
                    bedfile.write("\t")
                bedfile.write(str(target[-1]))
                bedfile.write("\n")
        self.compylog.info(
            f"Saved .bed file RecoveredBedID_{self.bedid}.bed to path "
            f"{self.out}"
        )
        
        #Save bedfile as .xlsx
        pathExcel = self.out + f"RecoveredBedID_{self.bedid}.xlsx"
        dfBedFile = DBManager.ExtractData(
        "Bedfiles", self.pathDB, bedid = self.bedid, df = True
        )
        workbook = xlsxwriter.Workbook(pathExcel)
        worksheet = workbook.add_worksheet("Targets")
        
        #Write table header
        names = dfBedFile.columns
        data = dfBedFile.values
        for headernum, headerrow in enumerate(names):
            worksheet.write_string(0, headernum, headerrow)
            
        #Write data
        for row_num, row in enumerate(data):
            cCOL = 0
            for col in row:
                worksheet.write_string(row_num + 1, cCOL, str(col))
                cCOL += 1
        workbook.close() 
        self.compylog.info(
            f"Excel sheet {pathExcel}.xlsx was saved!"
        )
        
        
    