import os
import re
from Bio import Entrez
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

"""Begins a NCBI connection using your e-mail address.
	INPUTS: email --> e-mail address"""
def init(email):
	Entrez.email = email

"""Closes an open handle.
	INPUTS: data --> an open handle"""
def close(data):
	data.close();

"""Searches specific data into de NCBI database.
	INPUTS: database  --> is an string that specifies the database where you want to search.		Ex.: "nucleotide"
					term --> it's the name of the term that we want to search.		Ex.: "Human papillomavirus"
					email --> it's your e-mail address for begining the NCBI connection
	OUTPUTS: a list of sequence identifiers"""
def search(database, term, email, retmax=20):
	init(email);
	handle = Entrez.esearch(db=database, term=term, retmax=retmax);
	record = Entrez.read(handle);
	return record["IdList"];

"""Downloads the data based on sequence identifiers
	INPUTS: database  --> it's an string that specifies the database where you want to search.		Ex.: "nucleotide"
					identifier --> it's a string of identifiers separated by commas Ex.: "LC027929,KC470291,FJ237042"
					file_format --> it's a string that specifies the data format.
					email --> it's your e-mail address for begining the NCBI connection
	OUTPUTS: a handle object"""
def download(database, identifier, file_format, email, retmax):
	init(email)
	handle = Entrez.efetch(db=database, id=identifier, rettype=file_format, retmode="text", retMax = 20)
	return handle

"""Parses a handle data.
	INPUTS: data --> it's a handle object. (It usually gotten by the download method)
					type_file --> it's a string that specifies the format of the handle data. 	Ex.:"genbank" 
	OUTPUTS: list of parsed sequnce (SeqIO list)"""
def parse(data,type_file):
	record = list(SeqIO.parse(data, type_file))
	return record

"""Saves the handle data into a file.
	INPUTS: filename --> it's the name which we want to save the file.
			data --> it's a handle data """
def save(filename,data):	
	if not os.path.isfile(filename):
		out_handle = open(filename, "w")
		out_handle.write(data.read())
		close(out_handle)
		close(data)
		print("Successfully saved")

"""Shows some information of a sequences list like: identifier, sequence and size.
	INPUTS: records_list --> a sequences list"""
def show_info(records_list):
	for i in range(len(records_list)):
		node = records_list[i];
		print("Id: " + node.id)
		print("Secuencia: " + repr(node.seq))
		print("Tam: " +str(len(node)) + "\n")

"""Changes a list of identifier into a string of identifiers separated by commas
	INPUTS: record_list -->a list of identifier
	OUTPUTS: string of identifiers separated by commas """
def format(records_list):
	cadena  = ""
	for i in range(len(records_list)):
		cadena = cadena + records_list[i]
		if i+1 != len(records_list):
			cadena = cadena + ","		
	return cadena

"""Gets the Open Reading Frame of a sequence.
	INPUTS: record_list --> it's a list of sequences (they usually are gotten by parse method)
					type_seq --> it's an string that specifies the search data. 		Ex.: "gene"
					ident --> it's an string that specifies the Open Reading Frame.		Ex.: "E6"""
def get_ORF(record_list,type_seq,ident):
	dicc = {type_seq : [ident]}
	seq_list = []
	for record in record_list:
		for feature in record.features:
			if feature.qualifiers == dicc:
				pos = str(feature.location)
				match = re.search('<?(\d+):>?(\d+)', pos)
				lower = int(match.group(1))
				uper = int(match.group(2))
				new = record[lower:uper].seq
				#print("Locus: " + record.id + "\n" + str(feature.qualifiers) + "\tLimite inferior: " + str(lower) + "\tLimite superior: " + str(uper))
				#print(new)
				#print("**********************************************************************")
				#new_record = newSequence(new,record.id,record.name,ident,record.dbxrefs)
				new_record = newSequence(new,record.id,record.name,record.description,record.dbxrefs)
				seq_list.append(new_record)
	return seq_list	
				

"""Crates a new sequence based on the data_sequence input
	INPUTS: data_sequence --> a string data
	OUTPUTS: sequence --> a new sequence created based the data_sequence input"""
def newSequence(data,ident,name,description,references):
	sequence = SeqRecord(data,ident,name,description,references)
	return sequence

"""Saves a SeqRecord list into a file.
	INPUTS: sequence_list --> SeqRecord list
					filename --> a string that specifies the filename
					format --> a string that specifies the format of the file. """
def writeFile(sequence_list,file_name,file_format):
	fasta_filename = file_name + "." + file_format;
	SeqIO.write(sequence_list, fasta_filename, file_format)
	return fasta_filename;

"""Checks if a list is empty.
	INPUTS: record_list --> a list
	OUTPUTS: boolean--> TRUE if is empty and False if it's not."""
def is_empty(record_list):
	if len(record_list) == 0:
		return True
	else:
		return False
	
"""Searches and downloads sequences by the term and database. Returns the file path they were saved."""
def downloadSequences(database, term, file_name, file_format, email, saving_path, retmax = 20):
	print ("Searching ...");
	records = search(database, term, email, retmax);
	print("Search finished");
	print("Downloading " + str(len(records)) + " sequences ...");
	records_str = format(records);
	record_handler = download(database, records_str, file_format, email, retmax);
	file_name_extension = file_name + '.' + file_format;
	if not os.path.exists(saving_path):
		os.makedirs(saving_path);
	if not os.path.exists(saving_path + "/" + file_name):
		os.makedirs(saving_path + "/" + file_name);
	saving_path += "/" + file_name + "/" + file_name_extension;
	save(saving_path, record_handler);
	print("Sequences stored in file " + saving_path);
	return saving_path;



