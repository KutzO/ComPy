import multiprocessing
import gzip
import pysam
from tqdm import tqdm
from .DbManager import DBManager
import logging


#Class to extract the VCF data   
class ExtractFromVCF():
    def __init__(self, threads, bedid, pathDB, version, 
                 lsCompatibleVersions, dicIDs, FileClass, dtime):
        self.threads = threads
        self.bedid = bedid
        self.pathDB = pathDB
        self.version = version
        self.lsVersions = lsCompatibleVersions
        self.dtime = dtime
        
        self.dicIDs = {}
        for intID in dicIDs.keys():
            if dicIDs[intID][2] == "vcf":
                self.dicIDs[intID] = dicIDs[intID]
        if FileClass:
            self.FileClass = FileClass
        else:
            self.FileClass = "Default"
        
        #Initiale logging tool
        self.comptoollog = logging.getLogger("ComparisonTool")
 
    

    
    ###Function to extract variants from parsed .vcf files
    #Uses HelpExtract() to split the extraction and save resources
    def ExtractVariants(self):
        if self.threads > len(self.dicIDs.values()):
            pool = multiprocessing.Pool(processes=len(self.dicIDs.values()))
        else:
            pool = multiprocessing.Pool(processes=self.threads)
        
        #Create subprocesses (one for each .vcf file)
        jobs = []
        for intID in self.dicIDs.keys():
            try:
                jobs.append(
                    pool.apply_async(
                        self.HelpExtract, args=(
                                        self.dicIDs[intID][1],
                                        self.dicIDs[intID][0], 
                                        intID
                                          )
                    )
                )
            except Exception as e:
                self.comptoollog.info("An error occured at extracting variants!")
                self.comptoollog.exception(e)
                pass
        
        #Start the subprocesses and call DBManager.InsertData() defined in script DbManager.py
        self.comptoollog.info("Extracting information from given VCF files!")
        print("Extract data from .vcf files")
        for _ in tqdm(jobs):
            result, intID = _.get()
            #print(result[0])
            DBManager.InsertData(result,"Extracted_Variants", self.pathDB)
            #DBManager.InsertData(sortedout, "Sorted_out", self.pathDB)
            DBManager.SecurityChangeData(
                self.pathDB, intID, self.bedid, self.version, "vcf"
            )
        pool.close()
        
        
        
    def HelpExtract(self, file, name, intID):
        ###Initialize important variables
        lsResult = []
        #lsSortOut = []
        strID = intID
        #Open gzipped VCF and convert INFO fields (Float, Integer) to String
        with pysam.VariantFile(file) as vcffile:
            #Start variant extraction
            
            for rec in vcffile:
                SAMPLE_COUNTER = 0
                for sample in rec.samples:
                    strName = str(name)                     #Sample name
                    strId = str(rec.id)                     #rs number (if given, "." otherwise)
                    if strId == None:
                        strId = "."
                    strChr = str(rec.chrom)                 #Chromosome
                    strPos = str(rec.pos)                   #Position
                    strRef = str(rec.ref)                   #Reference base
                    strLen = str(rec.rlen)
                    strANN = str(rec.info["ANN"])
                    strSample = sample
                                        
                    strAlt = rec.alts                   #Variant base
                    COUNTER_ALT_ALLELES = 1
                    for alt in strAlt:
                        try:
                            strType = str(rec.INFO["TYPE"])     #Variant type (Ins/Del/SNP)
                        except:
                            lenRef = len(rec.ref)
                            lenAlt = len(alt)
                            if lenRef == lenAlt:
                                strType = "snp"
                            elif lenRef > lenAlt:
                                strType = "del"
                            elif lenRef < lenAlt:
                                strType = "ins"
                        strGenot = rec.samples[SAMPLE_COUNTER]["GT"]
                        strGenot = str(
                            f"{strGenot[0]}/{strGenot[COUNTER_ALT_ALLELES]}"
                        )
                        intDP = rec.samples[SAMPLE_COUNTER]["DP"]                #Total mapped reads
                        
                        allAD = rec.samples[SAMPLE_COUNTER]["AD"]
                        if len(allAD) > 1:
                            try:
                                intAD = allAD[COUNTER_ALT_ALLELES]                #Reads carrying the variant
                            except:
                                print(rec.samples[SAMPLE_COUNTER]["AD"])
                                print(strChr)
                                print(strGenot)
                                print(strAlt)
                                print(alt)
                        else:
                            intAD = allAD[0]
                        COUNTER_ALT_ALLELES += 1
                        
                        lsResult.append(
                            [strName, strID, self.bedid, strId, strChr, 
                            strPos, strLen, strRef, alt, strType,
                            strGenot, intDP, intAD, strSample, strANN, 
                            self.FileClass]
                        )
                        
                    SAMPLE_COUNTER += 1
        return lsResult, strID      
        
        
        
    # def HelpExtract(self, file, name, intID):
    #     ###Initialize important variables
    #     lsResult = []
    #     lsSortOut = []
    #     strID = intID
    #     #Open gzipped VCF and convert INFO fields (Float, Integer) to String
    #     with pysam.VariantFile(file) as vcffile:
    #         #Start variant extraction
    #         for rec in vcffile:
    #             strName = str(name)                     #Sample name
    #             strId = str(rec.id)                     #rs number (if given, "." otherwise)
    #             if strId == None:
    #                 strId = "."
    #             strChr = str(rec.chrom)                 #Chromosome
    #             strPos = str(rec.pos)                   #Position
    #             strRef = str(rec.ref)                   #Reference base
    #             strAlt = str(rec.alts)                   #Variant base
    #             strLen = str(rec.rlen)
    #             try:
    #                 strType = str(rec.INFO["TYPE"])     #Variant type (Ins/Del/SNP)
    #             except:
    #                 len1 = len(rec.ref)
    #                 len2 = len(rec.alts)
    #                 if len1 == len2:
    #                     strType = "snp"
    #                 elif len1 > len2:
    #                     strType = "del"
    #                 elif len1 < len2:
    #                     strType = "ins"
                        
    #             strGenot = str(rec.samples[0]["GT"])    #Genotype
                
    #             #Which symbol is used in genotype entry?
    #             GTsplitter = ["/","|"]                   
    #             GTsplitter = [x for x in GTsplitter if x in strGenot][0]
                
    #             #Extracting gene name 
    #             try:
    #                 strGene = str(rec.INFO["GENEINFO"])
    #             except:
    #                 gene = [x for x in rec.INFO["ANN"]]
    #                 for i in gene:
    #                     strGene = i.split("|")[3]
    #                     if len(strGene) < 2:
    #                         strGene = i.split("|")[5]
    #                         break
    #                     else:
    #                         break
                            
    #             #Calculate allele frequency
    #             """
    #             Now it is important to explain some points, since they may cause trouble!
    #                 1) AD has two different types
    #                     a) type Integer
    #                     b) type List
    #                 2) AD is type list and has two values
    #                     a) GT is 1|0
    #                         --> only one AD has to be used for AF calculation
    #                     b) GT is 1|2
    #                         --> both AD values has to be used for AF calculation   
    #                         --> DB entry has to be splitted
    #                 3) AD is type list and has 3 values
    #                     a) GT is 0|1|2
    #                         --> There are two different Alleles
    #                         --> DB entry has to be splitted
    #                     b) GT is 1|2
    #                         --> both AD values has to be used for AF calculation   
    #                         --> DB entry has to be splitted
    #                     c) GT is 1|0
    #                         --> TROUBLE! The entry has to be saved in the "Sorted_out" database
    #             """
    #             tmpSD = rec.samples[0]["DP"]                #Total mapped reads
    #             tmpAD = rec.samples[0]["AD"]                #Reads carrying the variant
                
    #             #Check which case is given
    #             if type(tmpAD) == list:
    #                 if len(tmpAD) > 2:
    #                     tmpAF = tmpAD[1:]                   #Slice the AD values since the first belongs to reference!
    #                     #Case 3c
    #                     if len(rec.samples[0]["GT"].split(GTsplitter)) < 2 \
    #                         and "2" not in rec.samples[0]["GT"].split(GTsplitter):          
    #                             strAF = tmpAD.append(rec.samples[0]["GT"])
    #                             tmpSortOut = [
    #                                 strName, self.bedid, strID, strId, strChr,
    #                                 strPos, strRef, strAlt, strType, strGenot,
    #                                 strGene, strAF, 
    #                                 self.FileClass
    #                             ]
    #                             lsSortOut.append(tmpSortOut)
                                
    #                     #Case 3a                                                      
    #                     elif len(rec.samples[0]["GT"].split(GTsplitter)) < 2 \
    #                         and "2" in rec.samples[0]["GT"].split(GTsplitter):             
    #                             for epoch in range(0,2):        #Split database entries
    #                                 strAF = str(
    #                                     round(int(tmpAF[epoch]) / int(tmpSD), 
    #                                           3
    #                                     )
    #                                 )
    #                                 lsRowAddAll = [
    #                                     strName, self.bedid, strID, strId, strChr, 
    #                                     strPos, strRef, strAlt, strType, 
    #                                     strGenot, strGene, strAF, 
    #                                     self.FileClass
    #                                 ]
    #                                 lsResult.append(lsRowAddAll)
                                    
                        
    #                     #Case 3b)
    #                     elif len(rec.samples[0]["GT"].split(GTsplitter)) > 2:
    #                         for epoch in range(0,2):        #Split database entries
    #                             strAF = str(
    #                                 round(int(tmpAF[epoch])/int(tmpSD), 3)
    #                             )
    #                             lsRowAddAll = [
    #                                 strName, self.bedid, strID, strId, strChr, strPos,
    #                                 strRef, strAlt, strType, strGenot, strGene,
    #                                 strAF, self.FileClass
    #                             ]
    #                             lsResult.append(lsRowAddAll)
                    
                    
    #                 #Case 2a
    #                 elif len(tmpAD) == 2 \
    #                     and "2" in rec.samples[0]["GT"].split(GTsplitter):
    #                         for epoch in range(0,2):        #Split database entries
    #                             strAF = str(
    #                                 round(int(tmpAD[epoch])/int(tmpSD), 3)
    #                             )
    #                             lsRowAddAll = [
    #                                 strName, self.bedid,strID, strId, strChr, strPos,
    #                                 strRef, strAlt, strType, strGenot, strGene,
    #                                 strAF, self.FileClass
    #                             ]
    #                             lsResult.append(lsRowAddAll)
                    
    #                 #Case 2b
    #                 elif len(tmpAD) == 2 and \
    #                     "0" in rec.samples[0]["GT"].split(GTsplitter): 
    #                         strAF = str(round(int(tmpAD[1])/int(tmpSD), 3))
    #                         lsRowAddAll = [
    #                             strName, self.bedid, strID, strId, strChr, strPos, 
    #                             strRef, strAlt, strType, strGenot, strGene, 
    #                             strAF, self.FileClass
    #                         ]
    #                         lsResult.append(lsRowAddAll)                
                    
    #                 #Case 1b
    #                 elif len(tmpAD) == 1:
    #                     strAF = str(round(int(tmpAD[0])/int(tmpSD), 3))
    #                     lsRowAddAll = [
    #                         strName, self.bedid, strID, strId, strChr, strPos, strRef,
    #                         strAlt, strType, strGenot, strGene, strAF, 
    #                         self.FileClass
    #                     ]
    #                     lsResult.append(lsRowAddAll)
                
    #             #Case 1a    
    #             else:
    #                 strAF = str(round(int(tmpAD)/int(tmpSD), 3))
    #                 lsRowAddAll = [
    #                     strName, self.bedid, strID, strId, strChr, strPos, strRef,
    #                     strAlt, strType, strGenot, strGene, strAF,
    #                     self.FileClass
    #                 ]
    #                 lsResult.append(lsRowAddAll)  

                
    #     return lsResult, lsSortOut, strID
        
        
        
        