# ComPy
A python based tool to compare .bam and .vcf files using a SQLite3 database.



# 1.1 How to install


ComPy has two possibilities to get the program running. But first of all you have to make sure your operating system is **LINUX** and nothing else, since ComPy was designed for this system only. This is mainly due to package incompatibility. 
If that prerequisite is met, you can either do an installation using Anaconda (**recommended!**). Anaconda is a package organization system using so called virtual environments to make it a lot easier to manage 

	conda create -n ENVNAME compy
	conda activate ENVNAME

and explore the main program functions using

	ComPy --help

Another way to install compy is to download the package directly from pypy.org via PIP installer

	pip install compy matplotlib>=3.3 numpy>=1.19 pandas>=1.2 pip>=20.3 pypdf2>=1.26 pysam>=0.15 python>=3.7,<3.10.0a0 pyyaml>=5.4 seaborn>=0.11 tqdm >=4.55 xlsxwriter >=1.3


___
**Please be aware of the fact, that ComPy will store the database and the logs always at _/home/YOURUSERNAME/ComPy/_ ! This is also the default output folder for any results if you not define a different output path.** A possibility to change that folder is not yet given but will be implemented in the future.
___


# 1.2 Testing

You may noticed that ComPy comes with two more folders: *additional_data* and *test_data*. Those folders are made to help you learning the program. Before you start the testing please make sure to download the **GRCh37 (hg19) reference genome** from any kind of source. You may download it from here [LINK](https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.13/) or here (ftp://ftp.ncbi.nlm.nih.gov/genomes/refseq/vertebrate_mammalian/Homo_sapiens/all_assembly_versions/GCF_000001405.25_GRCh37.p13/GCF_000001405.25_GRCh37.p13_genomic.fna.gz). Afterwards, clone the GIT repository to get access to the test data set. 

	git clone https://g27vmgitlab.med.tu-dresden.de/KUTZOLIVE/compy.git
	cd compy/
	
You may find there three different folders. The ComPy folder just include the raw program scripts and can be ignored or if you installed ComPy via Anaconda/PyPy you can delete this folder. Inside *test_data* you can find different file types, serving as input data for the program. You can find 5 different .bam (*binary alignment map*) files consist of ~100.000 reads. Each of the .bam files is indexed by bwa aligner. Furthermore, 7 different .vcf (*variant call format*) files and a .bed (*browser extendible data*) file is stored in the folder. Since several functions of ComPy ask for a .csv (*comma separated values*) file, you can find an example .csv for each of this functions in the folder *additional_data*. In this manual I will show you how to run the different functions using this test data set.

# 2.1 Compare .bam and .vcf files

ComPy provides the *compare* function which allows you to calculate with both .bam and .vcf files. 

	cd YOURPATHTOCOMPY/compy/
	ComPy compare all -b ./test_data/*.bam -v ./test_data/*.vcf -o ./ComPyOut/ -r YOURPATHTOREFERENCE -e ./test_data/testbed.bed -f ./additional_data/Example_VCFplot.csv -a ./additional_data/Example_NameList_compare.csv 

Here, ComPy will take all .bam and .vcf files from the *test_data* folder and save the results in the database *Extraction.db* at */home/USERNAME/ComPy/.database/*. The optional parameters used here are: .bed file reduction: False, Number of subsamples used for read QC: 2.300.000, file group: default, flag filtered: 772. 
>**Note:** It is the most basic command to get all possible output. Due to the optional parameter -f ComPy also draws vennplots defined by the given .csv file. Furthermore, ComPy takes the names provided by the -a parameter.


Executing this command will create a new folder *ComPyOut* at you current working directory. There you can find another folder *Result_THISDATE* containing a .pdf, a .csv and one more folder *BigCompare*. The .csv contains information about the calculated data and which names are used to plot the results. The *Date_SampleComparison.pdf* includes plots which allows you to compare your samples. Providing both .bam and .vcf files and a VCFplot.csv will result in 6 different Plot types included in this .pdf. Details about how to interpret those can be found in *Detailed_Explainings.pdf*.
Furthermore, you can find a detailed analysis about every sample in *BigCompare*. This detailed analysis will only be performed if **more then 5 samples** sharing the same .bed file, file group, number of subsamples, and .bed file reduction status are provided OR already saved in the database.

## 2.1.1 Compare either .bam or .vcf files
If you want to analyze .bam data only, you have to use a different command.

	ComPy compare bam *args **kwargs

Same if you want to analyze .vcf only. Then you have to use

	ComPy compare vcf *args **kwargs

# 3. Plot saved data again

ComPy provides a function to re-plot already saved data again. We can now use the following command to do so.

	ComPy plot -n ./additional_data/Example_NameList_plot.csv -f ./additional_data/Example_VCFplot.csv -o ./ComPyOut/

Here, ComPy will take all provided sample IDs (See point 4 to learn how to get those sample IDs) and re-calculate the plots. 
>**Note:** As an alternative you may also use the initial compare function to re-plot your data since ComPy will notice if a provided file is already in the database. But then you have to provide all files again!


# 4. Export from database

Since ComPy uses a SQLite3 based database to store results, it is quite complicated to get on the data without at least basic knowledge in SQLite3 or a related language. But no worries, ComPy provides a function to export data from the database to a user-friendly .xlsx table. This can be done by using a quite simple command.

	ComPy extract info -o ./ComPyOut/

Now ComPy will create a new subfolder *./ComPyOut/Extracted/DATE/* including 3 .xslx files. In *BamInfo.xlsx* and *VCFInfo.xlsx* you can find all meta-information (**especially the internal sample ID**) about all samples saved in the database. The third .xslx file, *Bedinfo.xlsx*, is about the different .bed files saved in the database.

If you want to get more details about those files and what ComPy calculated (e.g. for own plots or analysis) there is also a way to export those data.

	ComPy extract data -n ./additional_data/Example_Extraction_Deletion.csv -o ./ComPyOut/

Now, ComPy will create 4 different .xlsx files in *./ComPyOut/Extracted/DATE/*. Of those, 3 are related to saved .bam files. *QCmetrics.xlsx* gives you information about PHRED scores and read length distribution. *ReadMapping.xlsx* contains information about the total number of counted reads per target and *ReadStatistiks.xlsx* is about the coverage, GC content, mapped reads after filtering. The last .xlsx called *Extracted_Variants.xlsx* contains all information about the extracted .vcf files. More details you may find in *Detailed_Explainings.pdf*.

ComPy also allows you to re-build used .bed files.

	ComPy extract bedrec -b 1 -o ./ComPyOut/

Here, ComPy will re-build .bed file with ID number 1 and saves it to *./ComPyOut/Extracted/Date/*. The targets are saved as type .bed and also as type .xlsx.

Finally you can export the database itself.

	ComPy extract database -o ./ComPyOut/

>**Note:** Even after export ComPy will save data at the original database path!

# 5. Delete data from database

Besides export functions ComPy also provides a function to delete data from you database.

	ComPy remove files -n ./additional_data/Example_Extraction_Deletion.csv

Here, ComPy will remove every entry showing the ID defined in the provided .csv data. Of course, this function has to be used carefully since **no backup** of the data is stored elsewhere by default!
>**Note:** If you had to interrupt ComPy while saving data to the database, you DON'T have to check if remaining data fragments are inside the database. ComPy has a routine that deletes incomplete data fragments everytime you invoke any process!


# 6. Merge databases

Maybe you want to compare your data with those of your colleagues or them want to use your data? For this case (and ofc other) ComPy provides a function to merge databases. To do so, you need one of two prerequisites. Either the second database itself (you may export it with the export function) or all of the possible .xlsx files. If you provide the .xlsx files you have to make sure nothing was changed here, since ComPy needs them in the way it exports them. You need in total 7 database .xlsx files (*Bedinfo.xlsx, BamInfo.xlsx, QCmetrics.xlsx, ReadMapping.xlsx, ReadStatistiks.xlsx, VCFInfo.xlsx, Extracted_Variants.xlsx*) with exactly the naming they got by ComPy. Furthermore, you need all associated .bed files as *Bedfiles_Name.xlsx*. Those *Bedfiles_Name.xlsx* files have to be provided to the function separately.

	ComPy merge -x PATHTOBAMVCFXLSX -b PATHTOBEDXLSX

or as alternative

	ComPy merge -d PATHTOOTHERDATABASE


# 7. The style .yaml

You may have noticed, that you can find a *Default.yaml* file in the folder *additional_data*. With the help of this .yaml file you can control different plot parameter to adapt them as you like. Feel free to have a look at this. To use this .yaml you have to provide it in the *compare* or *plot* functions with the -y parameter.

# 8. An ERROR occurs

Even though we have been constantly testing during programming, it might happen that you encounter an error or bug. Upsi :D... Here are some suggestions what to do:

### 1 Stay calm and read the docs

It often happens that an error occurs due to wrong insert files. So please make sure that you provided the right kind of data at each parameter. If you are not sure just re-read the point of interest in this manual.

### 2 Look at the Error message

Some major problems are already caught by ComPy itself. If so, the error message will tell you what to do.

### 3 Open an issue

If all the suggestions have not helped at all you may open an issue. Please open issues concerning ComPy only here on GitHub to make sure I can help you asap. If you open an issue please copy-paste the **whole** ERROR message together with the **complete** command you used. Furthermore, there is a folder at */home/YOURUSERNAME/ComPy/* called *logs/* where every ComPy usage is protocolled (The folder is deleted everytime 100 logs are stored, so you may not run into a storage problem here). Please save the log with a fitting date somewhere that you may provide it to me (it helps me a lot to understand when and why the error occurs).
