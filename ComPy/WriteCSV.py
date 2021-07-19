from .DbManager import DBManager
import xlsxwriter
import getpass
import os
import logging



###Extract needed data from database in .xlsx format (Identifies .bam names and .bed ID)
#Only called by ComparisonTool.py if user parses argument -xlsx or --Excel
class WriteCSV():
    def __init__(self, outputpath, pathDB, dtime, ID = False, 
                 fileclass = False):
        self.pathDB = pathDB                                     #Path to database
        self.out = self.FindOutputPath(outputpath, dtime)        #Path were .xlsx should be stored
        self.intID = ID
        self.fileclass = fileclass
        
        #Initiale logging tool
        self.comptoollog = logging.getLogger("ComparisonTool")


    
    ###Define result output
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
                f"/home/{getpass.getuser()}/ComparisonTool/Extracted/{dtime}/"
            )
        
        #Create folders
        if not os.path.exists(strResultOut):
            os.makedirs(strResultOut)
        return strResultOut



    
    ###Main function for converting database to .xlsx
    #Is separately called by ComparisonTool.py 
    #Table name has to be provided by calling the function
    def WriteData(self, keytable): 
        
        ##Collect data and table header
        #See script DbManager.py for more information
        if keytable in ["BamInfo", "VCFInfo", "BedInfo"]:
            data, names = DBManager.ExtractData(
                keytable, self.pathDB, ID = self.intID, ALL = True
            )
        else:
            data = DBManager.ExtractData(
                keytable, self.pathDB, ID = self.intID, 
            )
            names = data.columns
            data = data.values

        ##If needed data can not be found in database
        if len(data) == 0:
            if keytable == "VCFInfo" or keytable == "BamInfo":
                self.comptoollog.warning(
                    f"EXCEL ERROR: No data was added to table {keytable} yet!"
                )
                return
            self.comptoollog.warning(
                "Excel ERROR: The given input files was not found in database "
                +f"table {keytable}!"
            )
            self.comptoollog.info(f"File IDs: {self.intID}")
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
        self.comptoollog.info(
            f"Excel sheet {self.out+keytable}.xlsx was saved!"
        )
        
    
    ###Function to write new .bed file (saved as tab separated .bed file)
    def RecoverBED(self):
        targets = DBManager.ExtractData(
            "Bedfiles", self.pathDB, bedid = self.bedid
        )
        with open(
                self.out + f"RecoveredBedID_{self.bedid}.bed","w"
                  ) as bedfile:
            for target in targets:
                for value in target[1:-1]:
                    bedfile.write(str(value))
                    bedfile.write("\t")
                bedfile.write(str(target[-1]))
                bedfile.write("\n")
        self.comptoollog.info(
            f"Saved .bed file RecoveredBedID_{self.bedid}.bed to path "
            f"{self.out}"
        )
        
        
        
    