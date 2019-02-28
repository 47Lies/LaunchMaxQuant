#
#	use files:
#	-mqparTemplate 
#	-Sample Description file
#
#	other Informations:
#	
#	-path2 Fasta File
#	-Job's Name
#	-MaxQuantWorkingDirectory
#	-Raw Files directory
#

import getopt
import sys
from lxml import etree
import os
import re
os.environ['LANG'] = "en_US.UTF-8"

def CreateOrClean(DirectoryName):
	if os.path.isdir(DirectoryName):
		print("Directory",DirectoryName,"exist\nProceed to cleaning")
		CMD="rm --recursive --verbose --force "+DirectoryName+"/*"
		os.system(CMD)
	else :
		os.makedirs(DirectoryName)

def Usage():
	print("python")

def main():
	try:
		#the option letter -h -o
		# o: means that o take an argument
		opts, args = getopt.getopt(sys.argv[1:], "ho:cm:r:s:p:f:w:vt:", ["help", "output=","clean","mqpar=","rawdirectory=","sample-description=","threads=","fasta=","working-directory=","verbose","temp=",])
	except getopt.GetoptError as err:
		# print help information and exit:
		print(err) # will print something like "option -a not recognized"
		Usage()
		sys.exit(2)
	FastaPath="/data/users/ltaing/DATA_TMP/ltaing/Tumeurs_M5/SP-Human-UP000005640-012018.fasta"
	MaxQuantWorkingDirectory="/data/users/ltaing/DATA_TMP/ltaing/Tumeurs_M5"
	Threads=24
	SampleDescription="/data/users/ltaing/DATA_TMP/ltaing/Tumeurs_M5/SampleDescription.txt"
	RAW_DIR="/data/users/ltaing/DATA_TMP/ltaing/Tumeurs_M5/RAW"
	Template="/data/users/ltaing/DATA_TMP/ltaing/Tumeurs_M5/mqpar.template.1.6.2.6.xml"
	OutputXML=MaxQuantWorkingDirectory+"/mqpar.Silac.AllSamples.1.6.2.6.xml"
	TempFolder=MaxQuantWorkingDirectory+"/TempFolder"
	output = None
	verbose = False
	clean = False
	for o, a in opts:
		if o == "-v":
			verbose = True
		elif o == "-c":
			clean = True
		elif o in ("-h", "--help"):
			Usage()
			sys.exit()
		elif o in ("-f", "--fasta"):
			FastaPath=a
		elif o in ("-w", "--working-directory"):
			MaxQuantWorkingDirectory=a
		elif o in ("-p", "--threads"):
			Threads=a
		elif o in ("-t", "--temp"):
			TempFolder=a
		elif o in ("-s", "--sample-description"):
			SampleDescription=a
		elif o in ("-r", "--raw-directory"):
			RAW_DIR=a
		elif o in ("-m", "--mqpar"):
			Template=a
		elif o in ("-o", "--output"):
			OutputXML = a
		else:
			assert False, "unhandled option"
	# ...
	tree = etree.parse(Template)
	Fasta = tree.xpath("/MaxQuantParams/fastaFiles/FastaFileInfo/fastaFilePath")[0]
	Fasta.text = FastaPath
	
	Session = tree.xpath("/MaxQuantParams/name")[0]
	Session.text="Pimprenelle"
	
	Mods= tree.xpath("/MaxQuantParams/parameterGroups/parameterGroup/labelMods")[0]
	ZeMods=etree.Element("string")
	ZeMods.text="Arg10;Lys8"
	Mods.append(ZeMods)
	MaxLab= tree.xpath("/MaxQuantParams/parameterGroups/parameterGroup/maxLabeledAa")[0]
	MaxLab.text="3"
	#get to the first element of multiplicity and redefine it
	Multiplicity= tree.xpath("/MaxQuantParams/parameterGroups/parameterGroup/multiplicity")[0]
	Multiplicity.text="2"
	#get to the first element of fixed modification and erased anything inside
	fixedModifications=tree.xpath("/MaxQuantParams/parameterGroups/parameterGroup/fixedModifications/string")[0]
	fixedModifications.getparent().remove(fixedModifications)
	#get to the first element of variable modification and add a new string element 
	variableModifications=tree.xpath("/MaxQuantParams/parameterGroups/parameterGroup/variableModifications")[0]
	VarMod=etree.Element("string")
	VarMod.text="Carbamidomethyl (C)"
	variableModifications.append(VarMod)
	
	#Define the fixed search folder and the fixed combined folder 
	JobFSF=MaxQuantWorkingDirectory+"/FixedSearchFolder"
	FSF = tree.xpath("/MaxQuantParams/fixedSearchFolder")[0]
	FSF.text=JobFSF
	if clean:
		CreateOrClean(JobFSF)
	
	
	JobFCF=MaxQuantWorkingDirectory+"/fixedCombinedFolder"
	FCF = tree.xpath("/MaxQuantParams/fixedCombinedFolder")[0]
	FCF.text=JobFCF
	if clean:
		CreateOrClean(JobFCF)
	
	TF=tree.xpath("/MaxQuantParams/tempFolder")[0]
	TF.text=TempFolder
	
	#get to the first element of numThreads and set the number of threads from the conf file
	N =  tree.xpath("/MaxQuantParams/numThreads")[0]
	N.text=str(Threads)
	
	#Empty enything in the file sections
	RawFiles = tree.xpath("/MaxQuantParams/filePaths")[0]
	RawFiles.clear()
	Fraction = tree.xpath("/MaxQuantParams/fractions")[0]
	Fraction.clear()
	Experiments = tree.xpath("/MaxQuantParams/experiments")[0]
	Experiments.clear()
	PTMS = tree.xpath("/MaxQuantParams/ptms")[0]
	PTMS.clear()
	PGI = tree.xpath("/MaxQuantParams/paramGroupIndices")[0]
	PGI.clear()
	
	with open(SampleDescription) as f:
		header_line = next(f)
		for line in f:
			LocalFile,LocalExperiment,LocalFraction,LocalPTMS,LocalParamGroupIndices=line.rstrip().split("\t")
			ZefilePath=etree.Element("string")
			LocalPath=RAW_DIR+"/"+LocalFile
			ZefilePath.text=LocalPath
			RawFiles.append(ZefilePath)
			
			#remove index if exists
			RmIndex=LocalPath.replace("raw","index")
			RmIndex="rm --recursive --verbose --force "+RmIndex
			if clean:
				os.system(RmIndex)
			
			#remove the directory if exists
			RmFolder=LocalPath.replace(".raw","")
			RmFolder="rm --recursive --verbose --force "+RmFolder
			if clean:
				os.system(RmIndex)
			
			#create a new element short with the current file fraction index name and add it to the index 
			Zefraction=etree.Element("short")
			Zefraction.text=LocalFraction
			Fraction.append(Zefraction)
			
			#create a new element string with the current experiment index name and add it to the experiment
			ZeExperiment=etree.Element("string")
			ZeExperiment.text=LocalExperiment
			Experiments.append(ZeExperiment)
			
			#create a new element boolean with the current PTMS and add it to the PTMS always false
			ZePTMS=etree.Element("boolean")
			ZePTMS.text=LocalPTMS
			PTMS.append(ZePTMS)
			
			#create a new element int with the current paramGroupIndices and add it to the paramGroupIndices always 0
			ZePGI=etree.Element("int")
			ZePGI.text=LocalParamGroupIndices
			PGI.append(ZePGI)
	tree.write(OutputXML,xml_declaration=True,encoding='UTF-8',pretty_print=True)


if __name__ == "__main__":
    main()






