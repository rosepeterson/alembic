import Bio
from Bio.Seq import Seq
from Bio import SeqIO
import Bio.SeqIO as IO
from Bio.SeqRecord import SeqRecord
import re
import csv
import argparse
import os
import subprocess
from random import seed
from random import randint

parser = argparse.ArgumentParser(description="Combine contigs from Alembic pipeline")
parser.add_argument('-f', action='store', dest='fasta_file', help='Input fasta file')
parser.add_argument('-g', action='store', dest='gap', help='Gap CSV')
parser.add_argument('-o', action='store', dest='out', help='Output fasta file name')


result = parser.parse_args()

def reverse_complement(dna):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', '\n':'\n', 'N':'N'}
    return ''.join([complement[base] for base in dna[::-1]])



with open(result.gap) as csv_file:
	csv_reader = csv.reader(csv_file, delimiter = ",")
	
	
	for column in csv_reader:
		
		if column[2] == "forward" and column[4] =="forward":
			if re.match("-", column[5]):
				record_dict= SeqIO.index(result.fasta_file, 'fasta')
				id,code = Seq(record_dict.get_raw(column[1]).decode()).split("\t")
				leftlength= len(code) + int(column[5])
				leftlength1 =len(code) + int(column[5]) + len(id)
				
				if leftlength < 0:
					line = Seq(record_dict.get_raw(column[3]).decode())
					line1= re.sub(r"\t", "\n", str(line))
				else:
					seq_list=Seq(record_dict.get_raw(column[1][0:int(leftlength1)]).decode())+Seq(record_dict.get_raw(column[3]).decode())
					line = re.sub(r"\n>.*\t", "", str(seq_list))
					line1= re.sub(r"\t", "\n", line)
					
				with open(result.out, 'a') as writer:
					writer.writelines(line1)
				writer.close()
				record_dict.close()
					
		if column[2] == "forward" and column[4] =="forward":
			if column[5].isdigit():
				record_dict= SeqIO.index(result.fasta_file, 'fasta')
				seq_list1=Seq(record_dict.get_raw(column[1]).decode()) +Seq(record_dict.get_raw(column[3]).decode())
				numsn = "N" * int(column[5])
				line= re.sub(r"(?<=[A-Z])\n>.*\t", numsn, str(seq_list1))
				line1= re.sub(r"\t", "\n", line)
				with open(result.out, 'a') as writer:
					writer.writelines(line1)
				writer.close()
				record_dict.close()
				

		if column[2] == "forward" and column[4] =="reverse":
			if column[5].isdigit():
				record_dict= SeqIO.index(result.fasta_file, 'fasta')
				rightr = Seq(record_dict.get_raw(column[3]).decode())
				rightside = re.sub(r">.*\t","", str(rightr))
				rightcomp = reverse_complement(rightside)
				seq_list1=Seq(record_dict.get_raw(column[1]).decode()) +rightcomp
				numsn = "N" * int(column[5])
				line = re.sub(r"\n\n", numsn, str(seq_list1))
				line1 = re.sub(r"\t", "\n", line) + "\n"
				
				with open(result.out, 'a') as writer:
					writer.writelines(line1)
				writer.close()
				record_dict.close()
								
				
		if column[2] == "forward" and column[4] =="reverse":
			if re.match("-", column[5]):
				record_dict= SeqIO.index(result.fasta_file, 'fasta')
				id,code = Seq(record_dict.get_raw(column[1]).decode()).split("\t")
				leftlength= len(code) + int(column[5])
				leftlength1 =len(code) + int(column[5]) + len(id)
				
				if leftlength < 0:
					rightr = Seq(record_dict.get_raw(column[3]).decode())
					rightside = re.sub(r">.*\t","", str(rightr))
					rightcomp = reverse_complement(rightside)
					line1 = ">"+str(randint(0,1000000000000))+rightcomp+"\n"
				else:
					rightr = Seq(record_dict.get_raw(column[3]).decode())
					rightside = re.sub(r">.*\t","", str(rightr))
					rightcomp = reverse_complement(rightside)
					seq_list1=Seq(record_dict.get_raw(column[1][0:int(leftlength1)]).decode()) +rightcomp
					line = re.sub(r"\n\n", "", str(seq_list1))
					line1 = re.sub(r"\t", "\n", line) + "\n"
				with open(result.out, 'a') as writer:
					writer.writelines(line1)
				writer.close()
				record_dict.close()
		

		if column[2] == "reverse" and column[4] =="reverse":
			if re.match("-", column[5]):
				record_dict= SeqIO.index(result.fasta_file, 'fasta')
				leftr = Seq(record_dict.get_raw(column[1]).decode())
				leftside = re.sub(r">.*\t","", str(leftr))
				leftlength= len(leftside) + int(column[5])
				
				if leftlength < 0:
					rightr  = Seq(record_dict.get_raw(column[3]).decode())
					rightside = re.sub(r">.*\t","", str(rightr))
					rightcomp = reverse_complement(rightside)
					line1 = ">"+str(randint(0,1000000000000))+rightcomp+"\n"
				else:
					leftr = Seq(record_dict.get_raw(column[1]).decode())
					leftside = re.sub(r">.*\t","", str(leftr))
					leftcomp = reverse_complement(leftside)
					rightr  = Seq(record_dict.get_raw(column[3]).decode())
					rightside = re.sub(r">.*\t","", str(rightr))
					rightcomp = reverse_complement(rightside)
					seq_list1=leftcomp[0:int(leftlength)] +rightcomp
					line1 = ">"+str(randint(0,1000000000000))+re.sub(r"(?<=[A-Z])\n(?=[A-Z])", "", str(seq_list1)+"\n")				
			with open(result.out, 'a') as writer:
				writer.writelines(line1)
			writer.close()
			record_dict.close()
					
	
								
		if column[2] == "reverse" and column[4] =="reverse":
			if column[5].isdigit():
				record_dict= SeqIO.index(result.fasta_file, 'fasta')
				leftr = Seq(record_dict.get_raw(column[1]).decode())
				leftside = re.sub(r">.*\t","", str(leftr))
				leftcomp = reverse_complement(leftside)
				rightr  = Seq(record_dict.get_raw(column[3]).decode())
				rightside = re.sub(r">.*\t","", str(rightr))
				rightcomp = reverse_complement(rightside)			
				seq_list1=leftcomp +rightcomp
				numsn = "N" * int(column[5])
				line1 = ">"+str(randint(0,1000000000000))+re.sub(r"(?<=[A-Z])\n(?=[A-Z])", numsn, str(seq_list1)+'\n')
				with open(result.out, 'a') as writer:
					writer.writelines(line1)
				writer.close()
				record_dict.close()
				

		if column[2] == "reverse" and column[4] =="forward":
			if column[5].isdigit():
				record_dict= SeqIO.index(result.fasta_file, 'fasta')
				leftr = Seq(record_dict.get_raw(column[1]).decode())
				leftside = re.sub(r">.*\t","", str(leftr))
				leftcomp = reverse_complement(leftside)			
				seq_list1=leftcomp + Seq(record_dict.get_raw(column[3]).decode())
				numsn = "N" * int(column[5])
				line1 = ">"+str(randint(0,1000000000000))+re.sub(r">.*\t", numsn, str(seq_list1))
				with open(result.out, 'a') as writer:
					writer.writelines(line1)
				writer.close()
				record_dict.close()
									
		if column[2] == "reverse" and column[4] =="forward":
			if re.match("-", column[5]):
				record_dict= SeqIO.index(result.fasta_file, 'fasta')
				leftr = Seq(record_dict.get_raw(column[1]).decode())
				leftside = re.sub(r">.*\t","", str(leftr))
				leftlength= len(leftside) + int(column[5])
		
				if leftlength < 0:
					line1 = re.sub('\t', '\n', str(Seq(record_dict.get_raw(column[3]).decode())))
					
				else:
					leftr = Seq(record_dict.get_raw(column[1]).decode())
					leftside = re.sub(r">.*\t","", str(leftr))
					leftcomp = reverse_complement(leftside)
					seq_list1=leftcomp[0:int(leftlength)] + Seq(record_dict.get_raw(column[3]).decode())
					line1 = ">"+str(randint(0,1000000000000))+re.sub(r"(?<=[A-Z])>.*\t", "", str(seq_list1))
			
			with open(result.out, 'a') as writer:
				writer.writelines(line1)
			writer.close()
			record_dict.close()


	



