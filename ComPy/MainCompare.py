###Import needed packages
#External packages
import multiprocessing
import pysam
from tqdm import tqdm
import logging  
import hashlib
import sys
import gzip

#Own scripts
from .DbManager import DBManager
from .Preparation import DataPreparation
from .ExtractFromBam import ExtractInfoData
from .ExtractFromVCF import ExtractFromVCF
from .PlotData import PlotTheData
from .PlotAll import PlotCompAll




###Helper tool for calculate MD5 - checksum
def TestHash(file):
    a = hashlib.md5()
    with open(file, "rb") as readfile:
        for line in readfile:
            a.update(line)
    return a.hexdigest(), file

"""
##Noticed keyword arguments:
argBamfiles, argVcffiles, argBedfile, argReference, argOutput, argDatabase, 
argReduce, argThreads, argSubsamples, argIndex, dtime, strVersion, 
lsCompatibleVersions, argNameTable, argVcfPlotTable, argClass
"""
def CompToolCompare(tool, **kwargs):
        
        ###Prepare logger
        compylog = logging.getLogger("ComPy")
        
        """
        Prepare the data
            1) Prepare database 
            2) Calculate checksums
            2) Prepare arguments
            3) Create index files
        """
        
        """
        1) Prepare database
        """
            
        ###Prepare database
        compylog.info("Prepare database")
        classDataBase = DBManager(DBpath=kwargs["argDatabase"])                
        
        
        """
        2) Calculate checksums
        """
        ###Calculate checksums        
        dicCheckSums = {}
        compylog.info("Calculate checksums")
        print("Create checksums (md5) of input data")
        
        #Calculate .bam checksums
        if tool in ["bam", "all"]:
            lsCheckSumBam = []
            pool = multiprocessing.Pool(processes=kwargs["argThreads"])
            Jobs = []
            for file in kwargs["argBamfiles"]:
                Jobs.append(pool.apply_async(TestHash, args=(file,)))
            for job in tqdm(Jobs):
                tmpCheckSum = job.get()
                lsCheckSumBam.append(tmpCheckSum)
            for calcsum in lsCheckSumBam:
                dicCheckSums[str(calcsum[0])] = [calcsum[1], "bam"]
                compylog.info(f"md5{calcsum[0]} : name:{calcsum[1]}")
            
            pool.close()

        #Calculate .vcf checksums
        if tool in ["vcf", "all"]:
            lsCheckSumVcf = []
            pool = multiprocessing.Pool(processes=kwargs["argThreads"])
            Jobs = []
            for file in kwargs["argVcffiles"]:
                Jobs.append(pool.apply_async(TestHash, args=(file,)))
            for job in tqdm(Jobs):
                tmpCheckSum = job.get()
                lsCheckSumVcf.append(tmpCheckSum)
            for calcsum in lsCheckSumVcf:
                dicCheckSums[calcsum[0]] = [calcsum[1], "vcf"]
                compylog.info(f"md5{calcsum[0]} : name:{calcsum[1]}")
            
            pool.close()



        """
        3) Prepare the data
        """
        ###Prepare arguments
        if tool == "all":
            compylog.info("Prepare arguments")
            classPrep = DataPreparation(
                "compare", bam=kwargs["argBamfiles"], 
                DBpath=classDataBase.pathDB, booDB=classDataBase.booDB, 
                bed=kwargs["argBedfile"], subsample=kwargs["argSubsamples"], 
                dtime=kwargs["dtime"], outputpath=kwargs["argOutput"], 
                vcf=kwargs["argVcffiles"], ReduceBed=kwargs["argReduce"], 
                version=kwargs["strVersion"], checksum = dicCheckSums, 
                nameTable = kwargs["argNameTable"], 
                FileClass= kwargs["argClass"], styleyaml = kwargs["argyaml"],
                flag = kwargs["argflag"]
            )     
        elif tool == "bam":
            compylog.info("Prepare arguments")
            classPrep = DataPreparation(
                "compare", bam=kwargs["argBamfiles"], 
                DBpath=classDataBase.pathDB, booDB=classDataBase.booDB, 
                bed=kwargs["argBedfile"], subsample=kwargs["argSubsamples"], 
                dtime=kwargs["dtime"], outputpath=kwargs["argOutput"], 
                ReduceBed=kwargs["argReduce"], version=kwargs["strVersion"], 
                checksum = dicCheckSums, nameTable = kwargs["argNameTable"], 
                FileClass= kwargs["argClass"], styleyaml = kwargs["argyaml"],
                flag = kwargs["argflag"]
            ) 
        elif tool == "vcf":
            compylog.info("Prepare arguments")
            classPrep = DataPreparation(
                "compare", DBpath = classDataBase.pathDB, 
                booDB = classDataBase.booDB, dtime = kwargs["dtime"], 
                outputpath = kwargs["argOutput"], vcf = kwargs["argVcffiles"], 
                version = kwargs["strVersion"], checksum = dicCheckSums, 
                nameTable = kwargs["argNameTable"], 
                FileClass = kwargs["argClass"], bed = kwargs["argBedfile"],
                styleyaml = kwargs["argyaml"]
            )     


        """
        4) Generate index files (if --index was parsed)
        """
        ###If index files has to be created
        if tool in ["bam", "all"]: 
            if kwargs["argIndex"]:
                compylog.info(
                    "Index files will be created!"
                )
                print("Creating index files")
                pool = multiprocessing.Pool(processes=kwargs["argThreads"])
                lsJobs = []
                for lsValues in classPrep.dicIDs["bam"].values():
                    if lsValues[0] in classPrep.bamnamestodb:
                        lsJobs.append(
                            pool.apply_async(
                                pysam.index, args = (lsValues[1],)
                            )
                        )
                for job in tqdm(lsJobs):
                    job.get()
                pool.close()
                compylog.info("Finished creating index file")
                    
    
        """
        START OF THE MAIN PROGRAM
        
            1) Extraction of .bam data  (if not done previous)
            2) Extraction of .vcf data  (if provided and not done previous)
            3) Plotting the data        
            4) Collect all plots to generate final output
            5) Compare each sample with whole database
        """
        
        """
        1) Extract from .bam files
        Extracting information from each bam file which is not yet listed in the database
            a) Extract the readnumbers (total, mapped, unmapped) at each chromosome given from the bed file
            b) Extract QC data 
                b.1) Gather read PHRED scores from each chromosome 
                        --> check if fwd or rev read
                        --> Reads are picked using the python pseudorandom generator
                b.2) Gather readlength (All reads!)
                b.3) Count number of mapped reads per target (considering all reads)
                b.4) Calculate mean GC content per target using reference
                b.5) Calculate mean Coverage per target 
                b.6) Calculate read length distribution
                b.7) Calculate mean PHRED score per readbase together with standard derivation
            c) Save all data to database
        """
        ###Initialize the class ExtractFromBam for each parsed .bam file
        if tool in ["all", "bam"]:
            listOfBamClasses = []
            for intID in classPrep.dicIDs["bam"].keys():
                if classPrep.dicIDs["bam"][intID][0] in classPrep.bamnamestodb:
                #Class ExtractInfoData is defined in script ExtractFromBam.py
                    ExtractProg = ExtractInfoData(
                        classPrep.dicIDs["bam"][intID][1], 
                        classPrep.dicIDs["bam"][intID][0], 
                        kwargs["argReference"], kwargs["argThreads"], 
                        classPrep.dictTargets, classPrep.bedid, 
                        classDataBase.pathDB, kwargs["argSubsamples"], 
                        kwargs["strVersion"], kwargs["lsCompatibleVersions"], 
                        kwargs["argReduce"], intID, kwargs["argClass"],
                        kwargs["dtime"], flag = kwargs["argflag"]
                    )
                    listOfBamClasses.append(ExtractProg)
            
            for BamClass in listOfBamClasses:
                compylog.info(
                    f"Start Extraction from BAMfile: {BamClass.bam}"
                )
                print(f"Start Extraction from BAMfile: {BamClass.bam}")
                lsMeanPhredFwd, lsStdDevPhredFwd,lsMeanPhredRev, \
                    lsStdDevPhredRev, dicReadLenFwd, \
                        dicReadLenRev = BamClass.ExtractQCdata()               #Collect QC metrics
                BamClass.SavePhredLen(
                    lsMeanPhredFwd, lsStdDevPhredFwd,lsMeanPhredRev, 
                    lsStdDevPhredRev, dicReadLenFwd, dicReadLenRev
                )                                                              #Save data



        """
        2) Extracting data from given .vcf files
        """
        if tool in ["all", "vcf"]:
            if len(classPrep.vcfnames) != 0: 
                compylog.info("Extracting .vcf files")
                print("Extracting .vcf files")
                VCFextraction = ExtractFromVCF( 
                    kwargs["argThreads"], classPrep.bedid, classDataBase.pathDB, 
                    kwargs["strVersion"], kwargs["lsCompatibleVersions"], 
                    classPrep.dicIDs, kwargs["argClass"], kwargs["dtime"]
                )
                VCFextraction.ExtractVariants()

        
        """
        3) Plotting the data
        """
        compylog.info("Plotting data")
        print("Plotting data")
        
        if tool == "all":
            compylog.info("Plotting all data")
            PlotClass = PlotTheData(
                tool, classPrep.outputpathTmp, classPrep.outputpath, 
                classDataBase.pathDB, threads = kwargs["argThreads"],  
                dicIDs = classPrep.dicIDs, 
                VcfPlotTable = kwargs["argVcfPlotTable"], 
                styleyaml = classPrep.styleyaml
            )

        if tool == "bam":
            compylog.info("Plotting only bam data")
            PlotClass = PlotTheData(
                tool, classPrep.outputpathTmp, classPrep.outputpath, 
                classDataBase.pathDB, threads = kwargs["argThreads"], 
                dicIDs = classPrep.dicIDs, styleyaml = classPrep.styleyaml
            )
       
        if tool == "vcf":
            compylog.info("Plotting only vcf data")
            PlotClass = PlotTheData(
                tool, classPrep.outputpathTmp, classPrep.outputpath, 
                classDataBase.pathDB, threads = kwargs["argThreads"], 
                dicIDs = classPrep.dicIDs, 
                VcfPlotTable = kwargs["argVcfPlotTable"],
                styleyaml = classPrep.styleyaml
            )
        
        """
        4) Merge all plots and generate final output
        """
        print("Create sample comparison")
        compylog.info("Merge plots")
        PlotClass.MergeAllPlots(kwargs["dtime"])

        
        """
        5) Compare each sample with whole database
        """
        if tool in ["all", "bam"]:
            print("Start individual comparison")
            compylog.info("Start individual comparison")
            PlotCompAll(
                classDataBase.pathDB, kwargs["dtime"], 
                classPrep.outputpath, kwargs["argThreads"], 
                classPrep.dicIDs, classPrep.styleyaml
            )