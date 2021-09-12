import multiprocessing
import gzip
import pysam
from tqdm import tqdm
from .DbManager import DBManager
import logging


###Class to extract the VCF data   
class ExtractFromVCF():
    def __init__(self, threads, bedid, pathDB, version, 
                 lsCompatibleVersions, dicIDs, FileClass, dtime):
        
        ###Define global variables
        self.threads = threads
        self.bedid = bedid
        self.pathDB = pathDB
        self.version = version
        self.lsVersions = lsCompatibleVersions
        self.dtime = dtime
        if FileClass:
            self.FileClass = FileClass
        else:
            self.FileClass = "Default"
            
        ###Identify .vcf IDs
        self.dicIDs = {}
        for intID in dicIDs["vcf"].keys():
            self.dicIDs[intID] = dicIDs["vcf"][intID]

        #Initiale logging tool
        self.compylog = logging.getLogger("ComPy")
 
    

    
    ###Function to extract variants from parsed .vcf files
    def ExtractVariants(self):
    
        self.compylog.info("Start extracting from .vcf")
        
        #Reduce threads if lower number of .vcf files provided
        if self.threads > len(self.dicIDs.values()):
            self.compylog.info(f"Using {len(self.dicIDs.values())} threads")
            pool = multiprocessing.Pool(processes = len(self.dicIDs.values()))
        else:
            self.compylog.info(f"Using {self.threads} threads")
            pool = multiprocessing.Pool(processes = self.threads)
        
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
                self.compylog.info("An error occured at extracting variants!")
                self.compylog.exception(e)
                pass

        print("Extract data from .vcf files")
        for _ in tqdm(jobs):
            result, intID = _.get()
            DBManager.InsertData(result,"Extracted_Variants", self.pathDB)
            DBManager.SecurityChangeData(
                self.pathDB, intID, self.version, "vcf"
            )
        pool.close()
        self.compylog.info("Finished extracting variants!")
        
    ###Helper tool for getting variant information
    def HelpExtract(self, file, name, intID):
        ###Initialize important variables
        lsResult = []
        strID = intID
        
        ##Start extraction
        with pysam.VariantFile(file) as vcffile:
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
                    strLen = str(rec.rlen)                  #Variant length
                    strANN = str(rec.info["ANN"])           #Aminoacid changes
                    strSample = sample
                                        
                    strAlt = rec.alts                       #Variants
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
                        strGenot = rec.samples[SAMPLE_COUNTER]["GT"]            #Genotype
                        strGenot = str(
                            f"{strGenot[0]}/{strGenot[COUNTER_ALT_ALLELES]}"
                        )
                        intDP = rec.samples[SAMPLE_COUNTER]["DP"]                #Total mapped reads
                        
                        allAD = rec.samples[SAMPLE_COUNTER]["AD"]
                        if len(allAD) > 1:
                            try:
                                intAD = allAD[COUNTER_ALT_ALLELES]                #Reads carrying the variant
                            except:
                                self.compylog.error(rec.samples[SAMPLE_COUNTER]["AD"])
                                self.compylog.info(file)
                                self.compylog.info(strChr)
                                self.compylog.info(strGenot)
                                self.compylog.info(strAlt)
                                self.compylog.info(alt)
                        else:
                            intAD = allAD[0]
                        COUNTER_ALT_ALLELES += 1
                        intAF = intAD / intDP                                   #Calculate AF
                        lsResult.append(
                            [strName, strID, strId, strChr, 
                            strPos, strLen, strRef, alt, strType,
                            strGenot, intDP, intAD, intAF, strSample, strANN]
                        )
                        
                    SAMPLE_COUNTER += 1
        return lsResult, strID      

        