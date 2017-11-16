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
from Bio import pairwise2;
from Bio.pairwise2 import format_alignment;

def main():
	# Configuration parameters
	platform = sys.platform;
	foldername = "SSUrRNA"; #"EF1a";#"COXII"; #"COXI"
	fileformat = "fasta";
	consensus_text = "consensus"
	aligning_format = "fasta"
	clustalw2_path = "./aligners/clustalw/" + platform  + "/clustalw2"; # TODO: cambiar la ruta dependiendo del sistema operativo
	output_align_format = ".aln"
	threshold = 0.7
	output_dir = "consensus";
	clusterNumbers = [6,17,35];

	crearConsensos(platform, foldername, fileformat, consensus_text, aligning_format, clustalw2_path, output_align_format, threshold, output_dir);
	unirClusterConsensus(foldername, output_dir, fileformat, clusterNumbers);
	tablaAlineamiento(foldername, output_dir, fileformat, clusterNumbers, clustalw2_path);

def tablaAlineamiento(foldername, output_dir, fileformat, clusterNumbers, clustalw2_path):
	for k in clusterNumbers:
		print("Calculando tabla para cluster " + str(k));
		pattern = ".+K-" + str(k) + ".+";
		for filename in glob.glob(foldername + "/" + output_dir + "/*500_consensus.fasta"):
			matchObj = re.match(pattern, filename);
			if matchObj:	
				secuencias = parse(filename, fileformat);
				header = "";
				contadorCabecera = 0;
				csv = "";
				for secuencia1 in secuencias:
					linea = secuencia1.id;
					header = ",";
					for secuencia2 in secuencias:
						if contadorCabecera == 0:
							header += secuencia2.id + ",";
						alignment = pairwise2.align.globalxx(secuencia1,secuencia2);
						alignmentResult = format_alignment(*alignment[0]);
						match_score = re.search("Score=(\d+)", alignmentResult);
						score = match_score.group(1);
						linea += "," + score;
					if contadorCabecera == 0:
						print(header);
					print(linea);
					contadorCabecera += 1;

def unirClusterConsensus(foldername, output_dir, fileformat, clusterNumbers):
	for k in clusterNumbers:
		pattern = ".+K-"+ str(k) + ".+";
		records = [];
		new_filename = None;
		for filename in glob.glob(foldername + "/" + output_dir + "/*." + fileformat):
			matchObj = re.match(pattern, filename);
			if matchObj:
				record = parse(filename, fileformat);
				records.append(record[0]);
				new_filename = re.match('.*(Marker.+500_).*',filename);
		writeFile(records, foldername + "/" + output_dir + "/" + new_filename.group(1) + "consensus", fileformat);
		
		
def crearConsensos(plataform, foldername, fileformat, consensus_text, aligning_format, clustalw2_path, output_align_format, threshold, output_dir):	
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
			alignment = get_align(new_filename + output_align_format, "clustal");
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



