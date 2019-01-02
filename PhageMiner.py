##########General Information##########
# PhageMiner Version 1.0
# 
# PhageMiner is a user-supervised semi-automated bioinformatics 
# tool for prophage identification in bacterial genomes
# 
# Author: Reza Rezaei Javan
# Copyright (C) 2019 University of Oxford
# E-mail: Reza.RezaeiJavan@ndm.ox.ac.uk

##########License Information##########
# This is free software; you can redistribute it under
# the terms of the GNU Lesser General Public License as
# published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
#  This is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY.See the GNU Lesser General
#  Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this software package; if not, write to the
# Free Software Foundation, Inc., 51 Franklin St, Fifth Floor,
# Boston, MA  02110-1301  USA

##########Import Packages##########
#Biopython package; required for extracting genebank files
from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
#Pandas; a data analysis library
import pandas as pd
from pandas import *
#Numpy; a package for scientific computing
import numpy as np
#Sklearn; a package for machine learning. Used here for Mean Shift Clustering
from sklearn.cluster import MeanShift, estimate_bandwidth
#defaultdict; used for storing data
from collections import defaultdict
#needed for creating folders and moving files to different directories
import os
#For removing folders from the directory
import shutil
#Reportlab; python library for generating PDFs and graphics.
from reportlab.lib import colors
from reportlab.lib.units import cm
#Bio.Graphics; needed for creatomg diagrams (images) of coding sequences
from Bio.Graphics import GenomeDiagram
#For command line and batch analyses etc.
import sys
#To check whether the software is running on Mac or Linux
import platform

########## Set parameters - defined by user ##########
#Default Sensetivity - Can be changed during analysis - Set 0.1 for high, 0.15 for medium, and 0.2 for low. 
sensitivity = 0.2 
#Default number of flanking genes requested - Can be changed during analysis 
Num_flanking_genes = 5
#Number of genes needed to for a cluster - Increase the value for higher sensetivity
Minimum_his_needed_for_a_cluster = 3
#Manual investigation (Yes or No)
Manual_investigation = "Yes"
#Manual classification (Yes or No)
Manual_classification = "Yes"
#store clusters with flanking in CSV format (Yes or No)
save_flanking_in_csv = "No"
#remove contigs folder to save space
remove_contigs_folder = "Yes"

########## Set parameters - General ##########
#File address
gbk_file = sys.argv[1]
#List of possible answers to question asked by the software
affirmative_answers = ("Yes","Y","yes","YES","y")
negative_answers = ("No","N","no","NO","n")
#What delimiter to use when creating files
ex_delimiter = str("\t")

########## Keywords for annotation based genome mining - defined by user ##########
#Terms to look for - using "|" format as needed for Pandas package 
phage_gene = ("Phage|phage|bacteriophage|PblB|Bacteriophage|prophage|Prophage|major tail|major head|major capsid|minor capsid|tail tape measure|lytic amidase|Paratox")
#Also, make a simple List in a pythonic way, for downstream analyses 
phage_gene_strings = ("Phage","phage","bacteriophage","PblB","Bacteriophage","prophage","Prophage","major tail","major head","major capsid","minor capsid","tail tape measure","lytic amidase", "Paratox")
#Terms to ignore 
ignore_gene = ("Macrophage", "macrophage", "shock")

########## Creating working folder using os library ##########
def check_if_folder_exist_else_REMOVE_and_RECREATE(name_of_the_folder):
	if os.path.exists(name_of_the_folder):
		shutil.rmtree(name_of_the_folder)
		os.makedirs(name_of_the_folder)
	else:
		os.makedirs(name_of_the_folder)

def OpenPDF(PDF_file):
	if platform.system() == "Linux":
		os.system("xdg-open {} 2>/dev/null".format(PDF_file))
	else:
		os.system("open {}".format(PDF_file))

def extract_contigs():
	global gbk_file, genome_quality
	#Create a working folder to store contigs. If it doesn't exist, create it
	if not os.path.exists("{}/contigs".format(folder_name)):
		os.makedirs("{}/contigs".format(folder_name))

	#Load to the multi record Genbank file, and save each contig into an individual gbk file
	#Keep count of the nuber of contigs
	contig_num = 0
	for seq_record in SeqIO.parse(gbk_file, "genbank"):
		contig_num = contig_num + 1
		try:
			SeqIO.write(seq_record, "{}/contigs/contig_{}.gbk".format(folder_name, contig_num), "genbank")
		except ValueError as e:
			if str(e) != "Locus identifier %r is too long":
				seq_record.name = seq_record.name[:10]
				SeqIO.write(seq_record, "{}/contigs/contig_{}.gbk".format(folder_name,contig_num), "genbank")

	#Combine the contigs together. Start with the first contig, and then add the following if needed
	combined_contigs = SeqIO.read("{}/contigs/contig_1.gbk".format(folder_name), "genbank")
	if contig_num > 1:
		#This keeps record of which contig it is adding. Excluding the first one, it adds the rest one by one
		for contig_index in xrange(1,contig_num):
			next_contig = SeqIO.read("{}/contigs/contig_{}.gbk".format(folder_name, contig_index+1), "genbank")
			combined_contigs = combined_contigs + ("N" * 50) + next_contig

	#Saves all the combined contigs into a single gbk file.
	SeqIO.write(combined_contigs, "{}/{}_edited.gbk".format(folder_name,file_name), "genbank")
	gbk_file = "{}/{}_edited.gbk".format(folder_name,file_name)
	
	if len(os.listdir("{}/contigs".format(folder_name))) > 50:
		if len(os.listdir("{}/contigs".format(folder_name))) > 150:
			print "\nWARNING: Very poor genome quality ({} contigs). PhageMiner is very unlikely to perform well.\n".format(len(os.listdir("{}/contigs".format(folder_name))))
			genome_quality = "VeryPoorQuality"
		else:
			print "\nWARNING: Poor genome quality ({} contigs). PhageMiner may not perform well.\n".format(len(os.listdir("{}/contigs".format(folder_name))))
			genome_quality = "PoorQuality"

	#remove the contigs folder to save space
	if remove_contigs_folder == "Yes":
		shutil.rmtree("{}/contigs".format(folder_name))
	else:
		pass

########## Machine learning ##########
def Cluster_using_Machine_learning():
#Use Machine Learning Library for MeanShift Clustering to identify cluster of phage genes on bacterial chromosome
	global X, labels, n_clusters_
	try:
		X = np.array(zip(index_of_hits,np.zeros(len(df))), dtype=np.int)
		bandwidth = estimate_bandwidth(X, quantile=sensitivity)
		ms = MeanShift(bandwidth=bandwidth, bin_seeding=True)
		ms.fit(X)
		labels = ms.labels_
		cluster_centers = ms.cluster_centers_
		labels_unique = np.unique(labels)
		n_clusters_ = len(labels_unique)

		########## Reporting phage-related gene clusters ##########
		print "Phage-related gene clusters identifed: \n"
		for k in range(n_clusters_):
		    my_members = labels == k
		    #Only report clusters that have at least a certain number of hits
		    if len(X[my_members, 0]) > Minimum_his_needed_for_a_cluster:
		    	print "cluster {0}: {1}".format(k+1, X[my_members, 0])
	except ValueError:
		print "Could not identify any cluster for {}".format(str(file_name))
		raise
	 
########## Preparing csv files ##########
def Prepare_CSV_files():
	global dictionary_hits
	#store clusters with flanking in CSV format (Yes or No)
	if save_flanking_in_csv == "Yes":
		#Create a folder to store data, if it does not already exist
		if not os.path.exists("{}/with_flanking".format(folder_name)):
			os.makedirs("{}/with_flanking".format(folder_name))
		#Open excel files with "W" in order to ensure the content is empty from any previous runs, and add the headlines
		with open("{}/with_flanking/all_hits.txt".format(folder_name), "w") as text_file:
			text_file.write((("Cluster Number" + (ex_delimiter) + "Name" + (ex_delimiter) + "Minimum" + (ex_delimiter) + "Maximum" + (ex_delimiter) + "Sequence" + "\n" )))
	if save_flanking_in_csv == "No":
		pass
	#Create a dictionary with cluster number as Key and the hit indexes as Values 
	dictionary_hits = defaultdict(list)
	for k in range(n_clusters_):
		my_members = labels == k
		index_ids_of_each_cluster = []
		#Only consider clusters that have at least a certain number of hits
		if len(X[my_members, 0]) > Minimum_his_needed_for_a_cluster:
			for key in X[my_members, 0]:
				dictionary_hits[k].append(key)

def Custom_flanking_regions(Custom_name, flanking_num, Min_num, Max_num):
	##Create a ditionary with cluster number as key, and then the flanking regions
	flanking_regions_and_hits_cus = defaultdict(list)
	#defining how many flanking genes to capture
	flanking_min_cord = Min_num-flanking_num
	flanking_max_cord = Max_num+flanking_num
	#If the cluster is inserted at the very begining or at the very end of the genome, then ignore the requested number of flanking genes by user, and capture as much as possible.
	if flanking_min_cord < min(df_list):
		flanking_min_cord = min(df_list)
	if flanking_max_cord > max(df_list):
		flanking_max_cord = max(df_list)
	#add values to the flanking_regions_and_hits dictionary - Cluster name as key, min and max cord. as values 
	flanking_regions_and_hits_cus[Custom_name].append(flanking_min_cord)
	flanking_regions_and_hits_cus[Custom_name].append(flanking_max_cord)
	gbk_min_cord = ((df.ix[flanking_regions_and_hits_cus[Custom_name][0], 'Minimum']))
	gbk_max_cord = ((df.ix[flanking_regions_and_hits_cus[Custom_name][1], 'Maximum']))

	for seq_record in SeqIO.parse(gbk_file, "genbank"):
		records = seq_record[gbk_min_cord:gbk_max_cord]
		records.description = '{}_{}'.format(file_name,Custom_name)

		records.id = 'cluster_custom'

		#Try to save file, but if you get an error that the Locus identifier is too long (could happen if genomes are annotated by Prokka), then shorten the locus identifier to 10 characters
		try:
			SeqIO.write(records, "{}/genebank/{}.gbk".format(folder_name,Custom_name), "genbank")
		except ValueError as e:
			if str(e) != "Locus identifier %r is too long":
				records.name = records.name[:10]
				SeqIO.write(records, "{}/genebank/{}.gbk".format(folder_name,Custom_name), "genbank")

		print "cluster_custom.gbk was created from location {}..{}".format(gbk_min_cord,gbk_max_cord)
		#### Creating figures for visualisation
		#Creating individual figures of clusters of
		record = SeqIO.read("{}/genebank/{}.gbk".format(folder_name,Custom_name), "genbank")
		gd_diagram = GenomeDiagram.Diagram(record.id)
		gd_track_for_features = gd_diagram.new_track(2, name="Name")
		gd_feature_set = gd_track_for_features.new_set()

		for feature in record.features:
			# Colour hits red and the rest blue
			#Try to see if the "product" qualifer exists, if yes use that
		    try:
			    if any(key in feature.qualifiers["product"][0] for key in phage_gene_strings):
			        color = colors.red
			    else:
			        color = colors.blue
			#if it has no product name, and hence an error is encountered, then just go with the default qualifier
		    except KeyError:
		    	if any(key in feature.qualifiers for key in phage_gene_strings):
		        	color = colors.red
		    	else:
		        	color = colors.blue

		    # Define how the CDS should lool like in the figure   
		    gd_feature_set.add_feature(feature, sigil="ARROW",
		                               color=color, label=True,
		                               label_size = 6, label_angle=90,label_position="middle",)
		#Highlight assembly gaps in green 
		for site, name, color in [("NNNNNNNNNNNNNNNNNNNN","assembly gap",colors.green)]:
		    index = 0
		    while True:
		        index  = record.seq.find(site, start=index)
		        if index == -1 : break
		        feature = SeqFeature(FeatureLocation(index, index+len(site)))
		        gd_feature_set.add_feature(feature, color=color, name=name,
		                                   label=True, label_size = 10,
		                                   label_color=color,label_position="middle",)
		        index += len(site)


		gd_diagram.draw(format="linear", pagesize='A4',fragment_size=1, fragments=1,
		                start=0, tracklines=0,track_size=0.2,x=0.05, y=0.35, end=len(record))
		gd_diagram.write("{}/pdf/{}.pdf".format(folder_name,Custom_name), "PDF")
			

		############

########## Capturing the flanking regions of hits ##########
def Capture_flanking_regions():
	##Create a ditionary with cluster number as key, and then the flanking regions
	global flanking_regions_and_hits
	flanking_regions_and_hits = defaultdict(list)
	for c_number in range(n_clusters_):
		my_members = labels == c_number
		if len(X[my_members, 0]) > Minimum_his_needed_for_a_cluster:
			#setting up empty lists for calculations
			hit_cord = [] 
			#getting the minimum and maximum cordination of hits from dictionary_hits
			for key in dictionary_hits[c_number]:
				hit_cord.append(key)
			#defining how many flanking genes to capture
			flanking_min_cord = min(hit_cord)-Num_flanking_genes
			flanking_max_cord = max(hit_cord)+Num_flanking_genes
			#If the cluster is inserted at the very begining or at the very end of the genome, then ignore the requested number of flanking genes by user, and capture as much as possible.
			if flanking_min_cord < min(df_list):
				flanking_min_cord = min(df_list)
			if flanking_max_cord > max(df_list):
				flanking_max_cord = max(df_list)
			#add values to the flanking_regions_and_hits dictionary - Cluster name as key, min and max cord. as values 
			flanking_regions_and_hits[c_number].append(flanking_min_cord)
			flanking_regions_and_hits[c_number].append(flanking_max_cord)
			#Saving files into a tab seperated file
			if save_flanking_in_csv == "Yes":
				for key in df_list[flanking_min_cord:flanking_max_cord]:
					with open("{}/with_flanking/all_hits.txt".format(folder_name), "a") as text_file:
						text_file.write((str(c_number+1)) + ex_delimiter + (df.ix[key, 'Name']) + ex_delimiter + str((df.ix[key, 'Minimum'])) + ex_delimiter + str((df.ix[key, 'Maximum']))+ ex_delimiter + str((df.ix[key, 'Sequence']))+ "\n" )
			if save_flanking_in_csv == "No":
				pass
	print "\n"*5
	print "Looking for {} flanking genes of each cluster".format(int(Num_flanking_genes))
	print "\n"

#For manual investigation: changing the name of identified cluster based on whether it is a full-length prophage (FP), satellite prophage (SP)  or simply phage-related (PR).
def classify_cluster(cluster_name, key_id):
	if Manual_classification == "Yes":
		while True:
			cluster_category = raw_input("\nWhat is this cluster? Full-length prophage (F), satellite prophage (S), not sure (NS)?\n")
			if cluster_category == "F":
				while True:
					Assembly_gap = raw_input("\nDoes the sequence have an assembly gap? \"Y\" or \"N\".\n")
					if Assembly_gap in negative_answers:
						os.rename(("{}/genebank/{}{}.gbk".format(folder_name, cluster_name, key_id+1)),("{}/genebank/{}_FP_{}{}.gbk".format(folder_name, file_name, cluster_name, key_id+1)))
						os.rename(("{}/pdf/{}{}.pdf".format(folder_name, cluster_name, key_id+1)),("{}/pdf/{}_FP_{}{}.pdf".format(folder_name, file_name, cluster_name, key_id+1)))
						break
					if Assembly_gap in affirmative_answers:
						os.rename(("{}/genebank/{}{}.gbk".format(folder_name, cluster_name, key_id+1)),("{}/genebank/{}_FP_AG_{}{}.gbk".format(folder_name, file_name, cluster_name, key_id+1)))
						os.rename(("{}/pdf/{}{}.pdf".format(folder_name, cluster_name, key_id+1)),("{}/pdf/{}_FP_AG_{}{}.pdf".format(folder_name, file_name, cluster_name, key_id+1)))
						break
					else:
						print "\nAnswer not recognised. Select \"Y\" or \"N\"\n"
				break
			if cluster_category == "S":
				while True:
					Assembly_gap = raw_input("\nDoes the sequence have an assembly gap? \"Y\" or \"N\".\n")
					if Assembly_gap in negative_answers:
						os.rename(("{}/genebank/{}{}.gbk".format(folder_name, cluster_name, key_id+1)),("{}/genebank/{}_SP_{}{}.gbk".format(folder_name, file_name, cluster_name, key_id+1)))
						os.rename(("{}/pdf/{}{}.pdf".format(folder_name, cluster_name, key_id+1)),("{}/pdf/{}_SP_{}{}.pdf".format(folder_name, file_name, cluster_name, key_id+1)))
						break
					if Assembly_gap in affirmative_answers:
						os.rename(("{}/genebank/{}{}.gbk".format(folder_name, cluster_name, key_id+1)),("{}/genebank/{}_SP_AG_{}{}.gbk".format(folder_name, file_name, cluster_name, key_id+1)))
						os.rename(("{}/pdf/{}{}.pdf".format(folder_name, cluster_name, key_id+1)),("{}/pdf/{}_SP_AG_{}{}.pdf".format(folder_name, file_name, cluster_name, key_id+1)))
						break
					else:
						print "\nAnswer not recognised. Select \"Y\" or \"N\"\n"
				break
			if cluster_category == "NS":
				while True:
					Assembly_gap = raw_input("\nDoes the sequence have an assembly gap? \"Y\" or \"N\".\n")
					if Assembly_gap in negative_answers:
						os.rename(("{}/genebank/{}{}.gbk".format(folder_name, cluster_name, key_id+1)),("{}/genebank/{}_PR_{}{}.gbk".format(folder_name, file_name, cluster_name, key_id+1)))
						os.rename(("{}/pdf/{}{}.pdf".format(folder_name, cluster_name, key_id+1)),("{}/pdf/{}_PR_{}{}.pdf".format(folder_name, file_name, cluster_name, key_id+1)))
						break
					if Assembly_gap in affirmative_answers:
						os.rename(("{}/genebank/{}{}.gbk".format(folder_name, cluster_name, key_id+1)),("{}/genebank/{}_PR_AG_{}{}.gbk".format(folder_name, file_name, cluster_name, key_id+1)))
						os.rename(("{}/pdf/{}{}.pdf".format(folder_name, cluster_name, key_id+1)),("{}/pdf/{}_PR_AG_{}{}.pdf".format(folder_name, file_name, cluster_name, key_id+1)))
						break
					else:
						print "\nAnswer not recognised. Select \"Y\" or \"N\"\n"
				break
			else:
				print "\nAnswer not recognised. Select \"F\", \"S\" or \"NS\"\n "	

#Create working folders to store data, if they do not already exist
def create_working_folders():
	if not os.path.exists("{}/genebank".format(folder_name)):
		os.makedirs("{}/genebank".format(folder_name))
	if not os.path.exists("{}/pdf/".format(folder_name)):
		os.makedirs("{}/pdf/".format(folder_name))

#Remove working folders to clean up, if they exist and are empty
def remove_empty_folders():
	if os.path.exists("{}/genebank".format(folder_name)):
		if len(os.listdir("{}/genebank".format(folder_name))) is 0:
			shutil.rmtree("{}/genebank".format(folder_name))
	if os.path.exists("{}/pdf".format(folder_name)):
		if len(os.listdir("{}/pdf".format(folder_name))) is 0:
			shutil.rmtree("{}/pdf".format(folder_name))

########## Save hits and flanking regions to genebank files ##########
def save_into_gbk():
	create_working_folders()
	for key in range(len(flanking_regions_and_hits)):
		for seq_record in SeqIO.parse(gbk_file, "genbank"):
			gbk_min_cord = ((df.ix[flanking_regions_and_hits[key][0], 'Minimum']))
			gbk_max_cord = ((df.ix[flanking_regions_and_hits[key][1], 'Maximum']))

			records = seq_record[gbk_min_cord:gbk_max_cord]
			records.description = 'Cluster_{}_{}'.format(key+1,file_name)

			records.id = 'Cluster_{}'.format(key+1)

			#Try to save file, but if you get an error that the Locus identifier is too long (could happen if genomes are annotated by Prokka), then shorten the locus identifier to 10 characters
			try:
				SeqIO.write(records, "{}/genebank/cluster_{}.gbk".format(folder_name, key+1), "genbank")
			except ValueError as e:
				if str(e) != "Locus identifier %r is too long":
					records.name = records.name[:10]
					SeqIO.write(records, "{}/genebank/cluster_{}.gbk".format(folder_name, key+1), "genbank")

			print "cluster_{}.gbk was created from location {}..{}".format(key+1,gbk_min_cord,gbk_max_cord)
			#### Creating figures for visualisation
			#Creating individual figures of clusters of
			record = SeqIO.read("{}/genebank/cluster_{}.gbk".format(folder_name, key+1), "genbank")
			gd_diagram = GenomeDiagram.Diagram(record.id)
			gd_track_for_features = gd_diagram.new_track(2, name="Name")
			gd_feature_set = gd_track_for_features.new_set()

			for feature in record.features:
				# Colour hits red and the rest blue
				#Try to see if the "product" qualifer exists, if yes use that
			    try:
				    if any(key in feature.qualifiers["product"][0] for key in phage_gene_strings):
				        color = colors.red
				    else:
				        color = colors.blue
				#if it has no product name, and hence an error is encountered, then just go with the default qualifier
			    except KeyError:
			    	if any(key in feature.qualifiers for key in phage_gene_strings):
			        	color = colors.red
			    	else:
			        	color = colors.blue

			    # Define how the CDS should lool like in the figure   
			    gd_feature_set.add_feature(feature, sigil="ARROW",
			                               color=color, label=True,
			                               label_size = 6, label_angle=90,label_position="middle",)
			#Highlight assembly gaps in green 
			for site, name, color in [("NNNNNNNNNNNNNNNNNNNN","assembly gap",colors.green)]:
			    index = 0
			    while True:
			        index  = record.seq.find(site, start=index)
			        if index == -1 : break
			        feature = SeqFeature(FeatureLocation(index, index+len(site)))
			        gd_feature_set.add_feature(feature, color=color, name=name,
			                                   label=True, label_size = 10,
			                                   label_color=color,label_position="middle",)
			        index += len(site)


			gd_diagram.draw(format="linear", pagesize='A4',fragment_size=1, fragments=1,
			                start=0, tracklines=0,track_size=0.2,x=0.05, y=0.35, end=len(record))
			gd_diagram.write("{}/pdf/cluster_{}.pdf".format(folder_name, key+1), "PDF")
			
			#A series of nested while/break to select the clusters
			if Manual_investigation == "Yes":
				print "checking cluster {}...".format(key+1)
				OpenPDF("{}/pdf/cluster_{}.pdf".format(folder_name,key+1))
				while True:
					Check_1 = raw_input("\nDoes everything look okay with cluster? (\"Y\" or \"N\")\n")
					if Check_1 in negative_answers:
						while True:
							Check_2 = raw_input("\nWould you like to modify (M) or delete (D) this cluster? \n")
							if Check_2 == "D":
								os.remove("{}/pdf/cluster_{}.pdf".format(folder_name, key+1))
								os.remove("{}/genebank/cluster_{}.gbk".format(folder_name, key+1))
								print "removed cluster_{}".format(key+1)
								break
							if Check_2 == "M":
								os.remove("{}/pdf/cluster_{}.pdf".format(folder_name, key+1))
								os.remove("{}/genebank/cluster_{}.gbk".format(folder_name, key+1))
								while True:
									Check_3 = raw_input("\nWould you like to change the number of flanking genes (F) or set custom region (C)? \n")
									if Check_3 == "F":
										print "\nThe current number of flanking genes is {}\n".format(Num_flanking_genes)
										Check_flanking_num = int(raw_input("How many flanking genes should be used? \n"))
										Custom_flanking_regions("cluster_{}".format(key+1),Check_flanking_num, int(flanking_regions_and_hits[key][0]),int(flanking_regions_and_hits[key][1]))
										OpenPDF("{}/pdf/cluster_{}.pdf".format(folder_name,key+1))
										while True:
											Question2 = raw_input("Do you want to keep this edited file (\"Y\" or \"N\")?\n")
											if Question2 in negative_answers:
												os.remove("{}/genebank/cluster_{}.gbk".format(folder_name, key+1))
												break
											if Question2 in affirmative_answers:
												classify_cluster("cluster_",key)
												break
											else:
												print "\nAnswer not recognised. Use \"Y\" or \"N\".\n"
										break	
									if Check_3 == "C":
										Create_custom_region()
										break
									else:
										print "\nAnswer not recognised. Use \"F\" or \"C\".\n"	
								break	
							else:
								print "\nAnswer not recognised. Use \"D\" for delete or \"M\" for modifying.\n"
								pass
						break
					if Check_1 in affirmative_answers:
						classify_cluster("cluster_",key)	
						break
					else:
						print "\nAnswer not recgonised. Use \"Y\" or \"N\"\n"
						pass				
				break
			############

def make_diagrams():
	record = SeqIO.read(gbk_file, "genbank")
	gd_diagram = GenomeDiagram.Diagram(record.id)
	gdt_features = gd_diagram.new_track(0, scale=1, scale_largeticks = -2 , height = 1, greytrack =1, scale_largetick_interval = 100000)
	gds_features = gdt_features.new_set()


	for key in range(len(flanking_regions_and_hits)):
		for seq_record in SeqIO.parse(gbk_file, "genbank"):
			hit_min_cord = ((df.ix[flanking_regions_and_hits[key][0], 'Minimum']))
			hit_max_cord = ((df.ix[flanking_regions_and_hits[key][1], 'Maximum']))
			hit_size = hit_max_cord - hit_min_cord
			if hit_size > 0:
				feature = SeqFeature(FeatureLocation(hit_min_cord, hit_max_cord))
				gds_features.add_feature(feature, name="Cluster_{}".format(key+1),label_size = 15, label_position="middle", label_color=colors.green, label=True, sigil="BOX")
	#Highlight assembly gaps in green 
	for site, name, color in [("NNNNNNNNNNNNNNNNNNNN"," ",colors.black)]:
	    index = 0
	    while True:
	        index  = record.seq.find(site, start=index)
	        if index == -1 : break
	        feature = SeqFeature(FeatureLocation(index, index+len(site)))
	        gds_features.add_feature(feature, color=colors.red, name=name,
	                                   label=True, label_size = 8,
	                                   label_color=color, label_position="middle")
	        index += len(site)


	gd_diagram.draw(format="circular", circular=True, circle_core = 0.5, pagesize='A4',fragment_size=0, fragments=0,
	                start=0, tracklines=1, track_size=0.5,x=0.05, y=0.05, end=len(record))
	gd_diagram.write("{}/overall.pdf".format(folder_name), "PDF")

#Keeps a log of analyses: File name, number of cluster identified, number of putative hits
def KeepLog():
	with open("Log.csv", "a") as text_file:
		try:
			text_file.write((str(gbk_file) + (ex_delimiter) + str(len(flanking_regions_and_hits)) + (ex_delimiter) + str(len(index_of_hits)) +  "\n" ))
		except:
			text_file.write((str(gbk_file) + (ex_delimiter) + "N/A" + (ex_delimiter) + str(len(index_of_hits)) + "\n" )) 
iteration_num = 0
def Create_custom_region():
	global iteration_num
	create_working_folders()
	while True:
		while True:
			Q_Min = raw_input("What is the minimum gene number? \n")
			Q_Max = raw_input("What is the maximum gene number? \n")
			if (Q_Min.isdigit() == True) and (Q_Max.isdigit() == True):
				break
			else:
				print "Inputs need to be numbers. Try again... \n"
		Custom_flanking_regions("custom_{}".format(iteration_num+1),Num_flanking_genes, int(Q_Min),int(Q_Max))
		OpenPDF("{}/pdf/custom_{}.pdf".format(folder_name,iteration_num+1))
		while True:
			Question2 = raw_input("\nDo you want to keep the custom file (\"Y\" or \"N\")\n")
			if Question2 in negative_answers:
				os.remove("{}/genebank/custom_{}.gbk".format(folder_name, iteration_num+1))
				os.remove("{}/pdf/custom_{}.pdf".format(folder_name, iteration_num+1))
				print "custom_{} file deleted".format(iteration_num+1)
				break
			if Question2 in affirmative_answers:
				classify_cluster("custom_", iteration_num)
				print "custom_{} file saved".format(iteration_num+1)
				iteration_num = iteration_num + 1
				break
			else:
				print "\nAnswer not recognised. Use \"Y\" or \"N\"\n"
		break

def run_analyses():
	Prepare_CSV_files()
	Capture_flanking_regions()
	save_into_gbk()
	make_diagrams()

print "\nNow investigating {} file \n".format(str(gbk_file))



#Set the name of the file and working folder - take the name from the file, ignoring extensions
file_name = gbk_file.rsplit('.', 1)[0]
folder_name = gbk_file.rsplit('.', 1)[0]
#Changed space to underscore to avoid incompatibility problems 
folder_name = folder_name.replace(" ", "_")
file_name = file_name.replace(" ", "_")
#Assess genome quality based on the number of contigs
genome_quality = "default"

########## Creating working folder using os library ##########
check_if_folder_exist_else_REMOVE_and_RECREATE(folder_name)
extract_contigs()
########## Creating the CSV file from genebank file ##########
#Create it with "W" in order to ensure the content is empty from any previous runs, and add the headlines
csv_file = "{}/CDS_list.txt".format(folder_name)
with open(csv_file, "w") as text_file:
	text_file.write(("Minimum" + (ex_delimiter) + "Maximum" + (ex_delimiter) + "Sequence" + (ex_delimiter) + "Name" + "\n" ))
#Parse the data from genebank file and append it to the CSV file
for rec in SeqIO.parse(gbk_file, "genbank"):
    if rec.features:
        for feature in rec.features:
            if feature.type == "CDS":
            	with open(csv_file, "a") as text_file:
            		#Try using product names (ideal if genomes were annotated with Prokka). Otherwise, use default CDS name
            		try:
            			text_file.write((str(feature.location.start) + (ex_delimiter) + str(feature.location.end) + (ex_delimiter) + str(feature.location.extract(rec).seq) + (ex_delimiter) + str(feature.qualifiers["product"][0]) + "\n" ))
            		except KeyError as e:
						if str(e) != "product":
							text_file.write((str(feature.location.start) + (ex_delimiter) + str(feature.location.end) + (ex_delimiter) + str(feature.location.extract(rec).seq) + (ex_delimiter) + str(feature.qualifiers) + "\n" ))

########## Load CSV file in a data structure using Pandas package ##########
df = pd.read_csv(csv_file,sep='\t')
df.head()
df_list = df.index.tolist()

########## Search for phage related genes in CDS list ##########
df_phage_genes_extract = df[df["Name"].str.contains(phage_gene)]
df_phage_genes_extract.to_csv("{}/phage_genes_hits.csv".format(folder_name))
index_of_hits = df_phage_genes_extract.index.tolist()
#Remove the following catches from the hits
for key in index_of_hits:
	if any(x in df.ix[key, 'Name'] for x in ignore_gene):
		index_of_hits.remove(key) 


#For visual checks 
print "Hits for putative phage related genes in the genome: \n"
for key in index_of_hits:
	print key, df.ix[key, 'Name'] 
print "\n"*5

#Try clustering with the default sensitivity. Check if the user is happy with the clustering, otherwise try again with a new sensitivity. 
while True:
	try:
		Cluster_using_Machine_learning()
		TryClusteringAgain = raw_input("\nAre you happy with this clustering? (\"Y\" or \"N\")\n")
		if TryClusteringAgain in negative_answers:
			sensitivity = float(raw_input("\nChoose a new sensitivity (\"0.1\" to \"0.9\")\n"))
		if TryClusteringAgain in affirmative_answers:
			run_analyses()
			break
	except ValueError:
		TryClusteringAgain = raw_input("\nWould you like to try again with a new sensitivity? (\"Y\" or \"N\")\n")
		if TryClusteringAgain in affirmative_answers:
			sensitivity = float(raw_input("\nChoose a new sensitivity (\"0.1\" to \"0.9\")\n"))
		if TryClusteringAgain in negative_answers:
			break
	


if Manual_investigation == "Yes":
	print "\n\nFinal check:"
	while True:
		Question1 = raw_input("\nWould you like to select another custom region? (\"Y\" or \"N\")\n")
		if Question1 in affirmative_answers:
			Create_custom_region()
		if Question1 in negative_answers:
			break
		
remove_empty_folders()

#Change folder name, if the quality of the genome is poor
if genome_quality != "default":
	os.rename(folder_name,"{}_{}".format(folder_name,genome_quality))
else:
	pass

KeepLog()

print "="*30
