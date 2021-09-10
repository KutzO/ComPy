##!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 22 15:27:16 2020

@author: kutzolive

"""


###Import needed packages
#External packages
import os
import time
import argparse
import logging 
import getpass 
from glob import glob

#Own scripts
from .MainCompare import CompToolCompare
from .MainExtract import CompToolExtract
from .MainDelete import CompToolDelete
from .MainPlot import CompToolPlot
from .MainMerge import CompToolMerge



 
"""
This is the main script, guiding through the workflow processes.
It includes:
    - the argparser to communicate with the user
    - defines all parameters
    - calling necessary scripts
    - time measurement
"""
def main():


######################################################################################################################################
##############                                          Define  argparser                                               ##############
######################################################################################################################################     

    ###Define Argparser arguments
    parser = argparse.ArgumentParser()
    parser.set_defaults(tool = "Main")
    subparser = parser.add_subparsers(help = ("Commands"))
    
    """
    COMPARE
    """
    ##Needed arguments for comparison (main program)
    parse_compare = subparser.add_parser(
        "compare", help = ("Main program"),
    )
    parse_compare_sub = parse_compare.add_subparsers(help = ("Subcommands"))
    
    #Tool to compare both .vcf and .bam files
    parse_compare_all = parse_compare_sub.add_parser(
        "all", help = ("Compares data from both .bam and .vcf files"), 
        add_help = False
    )
    parse_compare_all.set_defaults(tool="all")
    group_all_req = parse_compare_all.add_argument_group("required")
    group_all_req.add_argument(
        "-b", "--bam", nargs = "*", help = ("Path to bam files"), 
        required = True
    )
    group_all_req.add_argument(
        "-e", "--bed", type = str, help = ("Path to bed file"), 
        required = True
    )
    group_all_req.add_argument(
        "-r","--ref", type = str, help = ("Path to reference sequence"), 
        required = True
    )
    group_all_req.add_argument(
        "-v", "--vcf", nargs = "*", help = ("Path to vcf data"), 
        required = True
    )

    
    ###
    group_all_opt = parse_compare_all.add_argument_group("optional")
    group_all_opt.add_argument(
        "-c", "--fileclass", type = str, default = False, 
        help = (
            "Define .bam class (e.g. filtered, non-filtered) to make "
            +"separate comparison possible"
        )
    )
    group_all_opt.add_argument(
        "-o", "--out", type = str, default = False, 
        help = (
            "Output directory (default= current_working_directory)"
        )
    )
    group_all_opt.add_argument(
        "-t", "--threads", type = int, default=10, 
        help = ("Number of threads used (default=10)")
    )
    group_all_opt.add_argument(
        "-d", "--reduce", action = "store_true", default = False, 
        help = (
            "If given, number of targets in BED file will be reduced by 90 %%"
        )
    )
    group_all_opt.add_argument(
        "-i", "--index", action = "store_true", default = False, 
        help = (
            "If given, an indexfile will be produced for each bam file"
        )
    )
    group_all_opt.add_argument(
        "-n", "--number", type = int, default = 2_300_000, 
        help = (
            "Number of reads to be collected from each chromosome for QC "
            +"(default = 2_300_000)"
        )
    )
    group_all_opt.add_argument(
        "-p", "--db", type = str, default = False, 
        help = ("Provide path to previous prepared database")
    )
    group_all_opt.add_argument(
        "-a", "--assign", default = False, 
        help = (
            "A comma separated (.csv) table for assigning names to .bam/.vcf "
            +"files! (default = Bam1, Bam2, Bam3)"
        )
    )
    group_all_opt.add_argument(
        "-y", "--yaml", default = False,
        help = ("Path to the style .yaml.")
    )
    group_all_opt.add_argument(
        "-h", "--help", action = "help", 
        help = ("show this help message and exit")
    )
    group_all_opt.add_argument(
        "-f", "--figure", type = str, 
        help = (
            "Comma separated (.csv) table which .vcf data should be compared"
        )
    )
    
    ###
    #Tool to compare .bam files only
    parse_compare_bam = parse_compare_sub.add_parser(
        "bam", help = ("Compares data from .bam files"), add_help = False
    )
    parse_compare_bam.set_defaults(tool = "bam")
    group_bam_req = parse_compare_bam.add_argument_group("required")
    group_bam_req.add_argument(
        "-b", "--bam", nargs = "*", 
        help = ("Path to bam files"), 
        required = True
    )
    group_bam_req.add_argument(
        "-e", "--bed", type = str, help = ("Path to bed file"), 
        required = True
    )
    group_bam_req.add_argument(
        "-r","--ref", type = str, help = ("Path to reference sequence"), 
        required = True
    )
    
    ###
    group_bam_opt = parse_compare_bam.add_argument_group("optional")
    group_bam_opt.add_argument(
        "-h", "--help", action = "help", 
        help = ("show this help message and exit")
    )
    group_bam_opt.add_argument(
        "-c", "--fileclass", type = str, default = False, 
        help = (
            "Define .bam class (e.g. filtered, non-filtered) to make "
            +"separate comparison possible"
        )
    )
    group_bam_opt.add_argument(
        "-o", "--out", type = str, default = False, 
        help = ("Output directory (default= current_working_directory)")
    )
    group_bam_opt.add_argument(
        "-t", "--threads", type = int, default=10, 
        help = ("Number of threads used (default=10)")
    )
    group_bam_opt.add_argument(
        "-d", "--reduce", action = "store_true", default = False, 
        help = (
            "If given, number of targets in BED file will be reduced by 90 %%"
        )
    )
    group_bam_opt.add_argument(
        "-i", "--index", action = "store_true", default = False, 
        help = ("If given, an indexfile will be produced for each bam file")
    )
    group_bam_opt.add_argument(
        "-n", "--number", type = int, default = 2_300_000, 
        help = (
            "Number of reads to be collected from each chromosome for QC "
            +"(default = 2_300_000)"
        )
    )
    group_bam_opt.add_argument(
        "-p", "--db", type = str, default = False, 
        help = ("Provide path to previous prepared database")
    )
    group_bam_opt.add_argument(
        "-a", "--assign", default = False, 
        help = (
            "A comma separated (.csv) table for assigning names to .bam"
            +" files! (default = Bam1, Bam2, Bam3)"
        )
    )
    group_bam_opt.add_argument(
        "-y", "--yaml", default = False,
        help = ("Path to the style .yaml.")
    )
    
    
    
    ####
    #Tool to compare .vcf data only
    parse_compare_vcf = parse_compare_sub.add_parser(
        "vcf", help = ("Compares data from .vcf files"), add_help = False
    )
    parse_compare_vcf.set_defaults(tool = "vcf")
    group_vcf_req = parse_compare_vcf.add_argument_group("required")
    group_vcf_req.add_argument(
        "-v", "--vcf", nargs = "*", help = ("Path to vcf data"), 
        required = True
    )
    group_vcf_req.add_argument(
        "-e", "--bed", type = str, help = ("Path to bed file"), 
        required = True
    )
    
    
    ###
    group_vcf_opt = parse_compare_vcf.add_argument_group("optional")
    group_vcf_opt.add_argument(
        "-h", "--help", action = "help", 
        help = ("show this help message and exit")
    )
    group_vcf_opt.add_argument(
        "-c", "--fileclass", type = str, default = False, 
        help = (
            "Define .vcf class (e.g. filtered, non-filtered) to make "
            +"separate comparison possible"
        )
    )
    group_vcf_opt.add_argument(
        "-o", "--out", type = str, default = False, 
        help = ("Output directory (default= current_working_directory)")
    )
    group_vcf_opt.add_argument(
        "-t", "--threads", type = int, default = 10, 
        help = ("Number of threads used (default=10)")
    )
    group_vcf_opt.add_argument(
        "-p", "--db", type = str, default = False, 
        help = ("Provide path to previous prepared database")
    )
    group_vcf_opt.add_argument(
        "-a", "--assign", default = False, 
        help = (
            "A comma separated (.csv) table for assigning names to .vcf "
            +"files! (default = VCF1, VCF2, VCF3)"
        )
    )
    group_vcf_opt.add_argument(
        "-y", "--yaml", default = False,
        help = ("Path to the style .yaml.")
    )
    group_vcf_req.add_argument(
        "-f", "--figure", type = str, 
        help = (
            "Comma separated (.csv) table which .vcf data should be compared"
        )
    )
    
    """
    EXTRACT
    """
    ##Needed arguments for extracting data from database (sub program)
    parse_extract = subparser.add_parser(
        "extract", help = ("Convert database data to .xlsx")
    )
    parse_extract_sub = parse_extract.add_subparsers(help = ("Subcommands"))

    #Tool to extract .xlsx
    parse_extract_data = parse_extract_sub.add_parser(
        "data", help = ("Extracts data from database in .xlsx format"), 
        add_help = False
    )
    parse_extract_data.set_defaults(tool = "data")
        ###required
    group_data_req = parse_extract_data.add_argument_group(
        "Required! Either ids or namelist"
    )
    group_data_req.add_argument(
        "-b", "--bam", nargs = "*", type = str, default = False,
        help = ("The .bam file ids to be deleted!")
    )
    group_data_req.add_argument(
        "-v", "--vcf", nargs = "*", type = str, default = False,
        help = ("The .vcf file ids to be deleted!")
    )
    group_data_req.add_argument(
        "-n", "--namelist", type = str, default = False,
        help = (
            "A .csv file containing ids of data to be deleted. A default "
            "list can be found in the ComPy test folder."
        )
    )
    # group_data_req = parse_extract_data.add_argument_group("required")
    # group_data_req.add_argument(
    #     "-i", "--id", default = False, nargs = "*",
    #     help = (
    #         "The unique identifier given to the files of interest. It can be "
    #         +"investigated by executing the EXTRACT INFO command and "
    #         +"searching in the BamInfo or VCFinfo files."
    #     )
    # )
    # group_data_req.add_argument(
    #     "-l", "--list", default = False, 
    #     help = (
    #         "A table containing samples that should be extracted "
    #         +"(has to be SAME as EXTRACT INFO table output)"
    #     )
    # )
    
    
    ###
    group_data_opt = parse_extract_data.add_argument_group("optional")
    group_data_opt.add_argument(
        "-h", "--help", action = "help", 
        help = ("show this help message and exit")
    )
    group_data_opt.add_argument(
        "-o", "--out", default = False, 
        help = (
            "Output path were .xlsx should be stored "
            +"(default= current_working_directory)"
        )
    )
    group_data_opt.add_argument(
        "-p", "--db", default = False, 
        help = (
            "Path to database (default= "
            +"'current_working_directory/00_Database/Extraction.db')"
        )
    )
   
    
    
    
    ###
    #Tool to extract .bed file targets
    parse_extract_bedrec = parse_extract_sub.add_parser(
        "bedrec", help = (
            "Recovers .bed file targets previously saved to the database"
        ), 
        add_help=False
    )
    parse_extract_bedrec.set_defaults(tool = "bedrec")
    group_bedrec_req = parse_extract_bedrec.add_argument_group("required")
    group_bedrec_req.add_argument(
        "-b", "--bed", nargs = "*", 
        help = (
            "The .bed identifier of the .bed file you want to be recovered"
        ), 
        required = True
    )
    
    ###
    group_bedrec_opt = parse_extract_bedrec.add_argument_group("optional")
    group_bedrec_opt.add_argument(
        "-h", "--help", action = "help", 
        help = ("show this help message and exit")
    )
    group_bedrec_opt.add_argument(
        "-o", "--out", default = False, 
        help = (
            "Output path were recovered .bed file should be stored "
            +"(default= current working directory)"
        )
    )
    group_bedrec_opt.add_argument(
        "-p", "--db", default = False, 
        help = (
            "Path to the database (default= "
            +"'current_working_directory/00_Database/Extraction.db')"
        )
    )

    
    
    ###
    #Tool to extract info .xlsx (all .bam and .vcf included in current database)
    parse_extract_info = parse_extract_sub.add_parser(
        "info", 
        help = (
            "Extracts information about all .bam and .vcf files currently "
            +"available in the database"
        )
    )
    parse_extract_info.set_defaults(tool = "info")
    parse_extract_info.add_argument(
        "-o", "--out", default = False, 
        help = (
            "Output path were .xlsx files should be stored "
            +"(default= current working directory)"
        )
    )
    parse_extract_info.add_argument(
        "-p", "--db", default = False, 
        help = (
            "Path to the database (default= "
            +"'current_working_directory/00_Database/Extraction.db')"
        )
    )

    #Tool to export database from hidden folder
    parse_extract_database = parse_extract_sub.add_parser(
        "database", help = (
            "Exports database at current status towards a defined path!"
        )
    )
    parse_extract_database.set_defaults(tool = "database")
    parse_extract_database.add_argument(
        "-o", "--out", default = False, 
        help = ("Output path for exporting database. REQUIRED!"), 
        required = True
    )
    parse_extract_database.add_argument(
        "-p", "--db", default = False, 
        help = (
            "Path to the database (default= "
            +"'current_working_directory/00_Database/Extraction.db')"
        )
    )
    
    """
    DELETE DATA
    """
    
    ##Needed arguments for delete data from database
    parse_delete = subparser.add_parser(
        "remove", help = "Delete data from database"
    )
    parse_delete_sub = parse_delete.add_subparsers(help = ("Subcommands"))
    
    ###
    #Tool to clean data showing not comparable versions
    parse_clean = parse_delete_sub.add_parser(
        "clean", 
        help = (
            "Removes all database entries with a version not including in "
            +"current compatibility list"
        )
    )
    parse_clean.set_defaults(tool = "clean")
    parse_clean.add_argument(
        "-p", "--db", type = str, default = False, 
        help = (
            "Path to the database (default= "
            +"'current_working_directory/00_Database/Extraction.db')"
        )
    )
    
    
    
    ###
    #Tool to delete files
    parse_files = parse_delete_sub.add_parser(
        "files", 
        help = (
            "Delete .bam and .vcf files from database"
        ), 
        add_help = False
    )
    parse_files.set_defaults(tool = "files")
    # group_files_req = parse_files.add_argument_group(
    #     "REQUIRED! Either id & datatype or NameList.csv"
    # )
    ###required
    group_files_req = parse_files.add_argument_group(
        "Required! Either ids or namelist"
    )
    group_files_req.add_argument(
        "-b", "--bam", nargs = "*", type = str, default = False,
        help = ("The .bam file ids to be deleted!")
    )
    group_files_req.add_argument(
        "-v", "--vcf", nargs = "*", type = str, default = False,
        help = ("The .vcf file ids to be deleted!")
    )
    group_files_req.add_argument(
        "-n", "--namelist", type = str, default = False,
        help = (
            "A .csv file containing ids of data to be deleted. A default "
            "list can be found in the ComPy test folder."
        )
    )
    # group_files_req.add_argument(
    #     "-i", "--id", nargs = "*", default = False,
    #     help = ("Unique ID given to saved file (see extract info table)")
    # )
    # group_files_req.add_argument(
    #     "-t", "--type", type = str, default = False,
    #     help = (
    #         "Either bam (to delete bam file) or vcf (to delete vcf file)"
    #     )
    # )
    # group_files_req.add_argument(
    #     "-l", "--list", type = str, default = False,
    #     help=(
    #         "A .csv table exactly the same as produced by extract info! "
    #         +"Containing all data that has to be deleted. If provided, no "
    #         +"more information about the data has to be provided!"
    #     )
    # )
    
    
    ###
    group_files_opt = parse_files.add_argument_group("optional")
    group_files_opt.add_argument(
        "-h", "--help", action = "help", 
        help = ("show this help message and exit")
    )
    group_files_opt.add_argument(
        "-p", "--db", default = False, 
        help = (
            "Path to the database (default= "
            +"'current_working_directory/00_Database/Extraction.db')"
        )
    )


    """
    PLOT data
    """
    group_plot = subparser.add_parser(
        "plot", help = "Plot data from database", add_help = False
    )
    
    group_plot.set_defaults(tool = "plot")
    
    ###required
    group_plot_req = group_plot.add_argument_group(
        "Required! Either ids or namelist"
    )
    group_plot_req.add_argument(
        "-b", "--bam", nargs = "*", type = str, default = False,
        help = ("The .bam file ids to be plotted!")
    )
    group_plot_req.add_argument(
        "-v", "--vcf", nargs = "*", type = str, default = False,
        help = ("The .vcf file ids to be plotted!")
    )
    group_plot_req.add_argument(
        "-n", "--namelist", type = str, default = False,
        help = (
            "A .csv file containing ids of data to be plotted. A default "
            "list can be found in the ComPy test folder."
        )
    )
    
    ###optional
    group_plot_opt = group_plot.add_argument_group(
        "Optional"    
    )
    group_plot_opt.add_argument(
        "-o", "--out", default = False, type = str,
        help = ("Define output path to store files at")
    )
    group_plot_opt.add_argument(
        "-db", "--database", default = False, type = str,
        help = ("Define path to database, if not stored at default!")
    )
    group_plot_opt.add_argument(
        "-t", "--threads", default = 10, type = int, 
        help = ("Define number of threads used for creating plots")
    )
    group_plot_opt.add_argument(
        "-f", "--figure", default = False, type = str,
        help = ("A .csv file defining comparable .vcf files.")
    )
    group_plot_opt.add_argument(
        "-y", "--yaml", default = False,
        help = ("Path to the style .yaml.")
    )
    group_plot_opt.add_argument(
        "-h", "--help", action = "help", 
        help = ("show this help message and exit")
    )
    group_plot_opt.add_argument(
        "-c", "--fileclasses", nargs = "*", default = False, 
        help = ("If different classes should be included in one plot.")
    )
    
    
    """
    MERGE data
    """
    group_merge = subparser.add_parser(
        "merge", help = ("Extend existing database with external data"),
        add_help = False
    )
    group_merge.set_defaults(tool = "merge")
    ###requiered
    group_merge_req = group_merge.add_argument_group(
        "Required! Either another db or files (.xlsx and .bed)"
    )
    group_merge_req.add_argument(
        "-d", "--database", type = str, default = False,
        help = ("Path to a database that should be included")
    )
    group_merge_req.add_argument(
        "-x", "--xlsx", nargs = "*", default = False,
        help = ("The .xlsx files that should be included. It is important "
                + "to use the default .xlsx names provided by the extract "
                + "data and extract info function!")
    )
    group_merge_req.add_argument(
        "-b", "--bed", nargs = "*", default = False,
        help = ("Path to bed file")
    )
    ###optional
    group_merge_opt = group_merge.add_argument_group(
        "Optional"    
    )
    group_merge_opt.add_argument(
        "-db", "--db", type = str, default = False,
        help = ("A path to the basic database, that should uptake data")
    )
    group_merge_opt.add_argument(
        "-h", "--help", action = "help", 
        help = ("show this help message and exit")
    )
    
    
    
    """
    PARSE ALL ARGS
    """
    ##Parse all given args
    args = parser.parse_args() 
        
######################################################################################################################################
##############                                      Define  global variables                                            ##############
######################################################################################################################################        
        
    ###Extract version number and list of compatible versions from __init__.py
    
    from . import __version__,__Compatible__
    strVersion = __version__
    lsCompatibleVersions = __Compatible__

    
    ###Start time measurement
    start = time.time()
    dtime = str(
        time.strftime(
            "%y%b%d%H%M"
        )
    ) 
    
    ###Start logging
    #Create logfolder
    if not os.path.exists(
            f"/home/{getpass.getuser()}/ComPy/logs/"
            ):
        os.makedirs(
            f"/home/{getpass.getuser()}/ComPy/logs/"
        )
    #Delete Folder if >100 logs
    pathlogfiles = glob(f"/home/{getpass.getuser()}/ComPy/logs/*")
    if len(pathlogfiles) > 100:
        for f in pathlogfiles:
            os.remove(f)
        
    
    #Define global log parameters
    logpath = f"/home/{getpass.getuser()}/ComPy/logs/{dtime}_log.txt"
    compylog = logging.getLogger(
        "ComPy"
    )
    logging.basicConfig(
        filename= logpath, format='%(asctime)s : %(levelname)s - %(message)s', 
        datefmt='%y-%m-%d %H:%M:%S'
    )
    compylog.setLevel(
        logging.INFO
    )
    
    
    
######################################################################################################################################
##############                                          Calling  programs                                               ##############
######################################################################################################################################   
     
    ###The main program to extract data from .bam / .vcf files and compare them 
    if args.tool == "all":
        CompToolCompare(
            "all", argBamfiles = args.bam, argVcffiles = args.vcf, 
            argBedfile = args.bed, argReference = args.ref, 
            argOutput = args.out, argDatabase = args.db, 
            argReduce = args.reduce, argThreads = args.threads, 
            argSubsamples = args.number, argIndex = args.index, 
            dtime = dtime, strVersion = strVersion, 
            lsCompatibleVersions = lsCompatibleVersions, 
            argNameTable = args.assign, argVcfPlotTable = args.figure, 
            argClass = args.fileclass, argyaml = args.yaml
        )
    if args.tool == "bam":
        CompToolCompare(
            "bam", argBamfiles = args.bam, argBedfile = args.bed, 
            argReference = args.ref, argOutput = args.out, 
            argDatabase = args.db, argReduce = args.reduce, 
            argThreads = args.threads, argSubsamples = args.number, 
            argIndex = args.index, dtime = dtime, strVersion = strVersion, 
            lsCompatibleVersions = lsCompatibleVersions, 
            argNameTable = args.assign, argClass = args.fileclass,
            argyaml = args.yaml
        )
    if args.tool == "vcf":
        CompToolCompare(
            "vcf", argVcffiles = args.vcf, argOutput = args.out, 
            argDatabase = args.db, argThreads = args.threads, dtime = dtime, 
            strVersion = strVersion, argBedfile = args.bed,
            lsCompatibleVersions = lsCompatibleVersions, 
            argVcfPlotTable = args.figure, argClass = args.fileclass,
            argNameTable = args.assign, argyaml = args.yaml
        )

    ###The program to extract data from an existing database
    
    #To extract data from database:
    if args.tool == "data":
        CompToolExtract(
            args.tool, argDatabase = args.db, dtime = dtime, 
            argOutput = args.out, strVersion = strVersion, 
            argNameList = args.namelist, argbamid = args.bam, argvcfid = args.vcf
        )
    
    #To recover .bed file targets
    if args.tool == "bedrec":
        CompToolExtract(
            args.tool, argDatabase = args.db, argBedfile = args.bed, 
            argOutput = args.out, dtime = dtime
        )
    
    #To extract info files
    if args.tool == "info":
        CompToolExtract(
            args.tool, argDatabase =  args.db, argOutput = args.out, 
            dtime = dtime
        )
    
    #To export database
    if args.tool == "database":
        CompToolExtract(
            args.tool, argOutput = args.out, dtime = dtime, 
            argDatabase = args.db
        )

        

    
    ###The program to delete data from an existing database
    if args.tool == "files":
        CompToolDelete(
            args.tool, argbamid = args.bam, argvcfid = args.vcf,  
            argDatabase = args.db, dtime= dtime, argNameTable = args.namelist
        )
        
    if args.tool == "clean":
        CompToolDelete(
            args.tool, argDatabase = args.db, dtime = dtime, 
            lsCompatibleVersions = lsCompatibleVersions
        )
    
    
    if args.tool == "plot":
        CompToolPlot(
            argbamid = args.bam, argvcfid = args.vcf, 
            argnamelist = args.namelist, argoutput = args.out, 
            argDatabase = args.database, dtime = dtime,
            argthreads = args.threads, argvcfplottable = args.figure, 
            argyaml = args.yaml, argThreads = args.threads,
            argclass = args.fileclasses
        )
    
    if args.tool == "merge":
        CompToolMerge(
            argxlsx = args.xlsx, argnewdb = args.database, argolddb = args.db,
            argbed = args.bed
        )
###############################################################################
###############################################################################
###############################################################################   
    
    
    ##Stop time measurement and print needed time for whole program
    stop = time.time()
        
    intresult = stop - start
    compylog.info(
        "FINISHED run trough"
    )
    compylog.info(
        f"Benötigte Zeit: {int(round(intresult,0))} Sekunden!"
    )
    compylog.info(
        f"Das entspricht ca.: {round(intresult/60,2)} Minuten!"
    )
    compylog.info(
        f"Das wiederum entspricht: {round(intresult/3600,2)} Stunden!"
    )

if __name__ == "__main__":
    main()
