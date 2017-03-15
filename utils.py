import re;
from Bio.SeqRecord import SeqRecord;

def filterByNCBITagValue(sequence_records, feature_name, feature_tag, tag_value):
	sequences_filtered = [];
	feature_tag_regex = '(' + '|'.join(feature_tag) + ')'; # Se construye la expresion regular para no dejar la responsabilidad al usuario.
	tag_value_regex = '(' + '|'.join(tag_value) + ')'; # Se construye la expresion regular para no dejar la responsabilidad al usuario
	for record in sequence_records:
		tag_value_found = "";
		for feature in record.features:
			n = re.search(feature_name,feature.type);
			if n is not None:
				qualifiers = feature.qualifiers;
				qualifiersStr = ','.join(qualifiers.keys());
				m = re.search(feature_tag_regex, qualifiersStr);
				if m is not None:
					tag_value_found += ','.join(qualifiers[m.group(0)]) + ',';
		values_found = re.findall(tag_value_regex, tag_value_found);
		save = True;
		for value in tag_value:
			save = save and (value in values_found);
		if save :
			sequences_filtered.append(record);
	return sequences_filtered; 

def separate_ORFs_per_sequence(sequence_data, feature_tag, tag_values):
	sequence_grouped_by_ORFs = {}; # Este diccionario tendra dos niveles: identificado de ORF y  Arreglo de ORFs
	for tag in tag_values: # Se inicializa el diccionario para que tenga los tag_values que requiere (E6, E7, ..., L1)
		sequence_grouped_by_ORFs[tag] = [];
	for sequence in sequence_data:
		id_sequence = sequence.id;
		for feature in sequence.features:
			if feature.type == feature_tag:
				orf = feature.qualifiers.get('product'); #pasar product como parametro y ver como preguntamos por gene
				orf = ''.join(orf); # Se cambia de arreglo a string para poder compararlo contra arreglo tag_values
				if orf in tag_values:
					pos = str(feature.location);
					match = re.search('<?(\d+):>?(\d+)', pos);
					lower = int(match.group(1));
					uper = int(match.group(2));
					orf_sequence = sequence[lower:uper].seq;
					orf_seqRecord = SeqRecord(orf_sequence);
					orf_seqRecord.id = sequence.id;
					orf_seqRecord.description = sequence.description;
					orf_seqRecord.name = orf
					sequence_grouped_by_ORFs[orf].append(orf_seqRecord);
	return sequence_grouped_by_ORFs;



