#!/usr/bin/env python3
from utils import *;
from Petitions import *;
from Alignments import *;
from Bio import Phylo;
from Bio import AlignIO;
from IO import *;
from time import time;
import sys
import glob, os;
import shutil;


def main():
	# Configuration parameters
	platform = sys.platform;
	foldername = "fastaCOXI";
	fileformat = "fasta";
	consensus_text = "consensus"
	aligning_format = "fasta"
	clustalw2_path = "./aligners/clustalw/" + platform  + "/clustalw2"; # TODO: cambiar la ruta dependiendo del sistema operativo
	output_aling_format = ".aln"
	threshold = 0.7
	output_dir = "consensus";

	
	if not os.path.exists(foldername + "/" + output_dir):
        	os.makedirs(foldername + "/" + output_dir);
	
	numero_cluster = re.compile('.+Cluster-(\d+)');
	
	# Leyendo archivos fasta del directorio	
	for filename in glob.glob(foldername + "/" + "*." + fileformat):
		record = parse(filename, fileformat);
		new_filename = filename.split(".")[0];
		id = new_filename.split("/")[1];
		matcher = numero_cluster.search(new_filename); 
		header = matcher.group(1);
		if len(record) > 1:
			print("Processing " + filename);
			print("Aligning ...");
			clustal_align(filename, clustalw2_path);
			alignment = get_align(new_filename + output_aling_format, "clustal");
			#Obteniendo secuencia consenso
			print("Geting consensus ...");
			consensus_seq_record = SeqRecord(get_consensus(alignment,threshold), id = header, description = "")
			print ("Saving " + new_filename + " consensus ...")		
			writeFile(consensus_seq_record, foldername + "/" + output_dir + "/" + id + "-consensus", aligning_format);
		else:
			record[0].id = header;
			record[0].description = "";
			writeFile(record, foldername + "/" + output_dir + "/" + id + "-consensus",aligning_format);



if __name__ == "__main__":
	main();



