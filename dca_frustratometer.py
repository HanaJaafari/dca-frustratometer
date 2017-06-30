import random
import numpy
import os.path
import argparse

scratchdir='/scratch/bs25'
basedir='/home/bs25'

def get_pfamID(pdbID, chain):
    import pandas as pd
    df = pd.read_table('%s/dca-frustratometer/pdb_chain_pfam.lst.txt' % basedir, header=1)
    if sum((df['PDB'] == pdbID.lower()) & (df['CHAIN'] == chain.upper())) != 0:
        pfamID=df.loc[(df['PDB'] == pdbID.lower()) & (df['CHAIN'] == chain.upper())]["PFAM_ID"].values[0]
    else:
        print('cant find pfamID')
        pfamID='null'
    return(pfamID)

def get_uniprotID(pdbID, chain):
    import urllib2
    response = urllib2.urlopen('http://www.bioinf.org.uk/cgi-bin/pdbsws/query.pl?plain=1&qtype=pdb&id=%s&chain=%s' % (pdbID, chain))
    html = response.read()
    info=html.split('\n')
    uniprotID=info[3].split(' ')[1]
    return uniprotID

def get_pfam_map(pdbID, chain):
    import pandas as pd
    df = pd.read_table('%s/dca-frustratometer/pdb_pfam_map.txt' % basedir, header=0)
    if sum((df['PDB_ID'] == pdbID.upper()) & (df['CHAIN_ID'] == chain.upper())) != 0:
        start=df.loc[(df['PDB_ID'] == pdbID.upper()) & (df['CHAIN_ID'] == chain.upper())]["PdbResNumStart"].values[0]
        end=df.loc[(df['PDB_ID'] == pdbID.upper()) & (df['CHAIN_ID'] == chain.upper())]["PdbResNumEnd"].values[0]
    else:
        print('data not found')
        pfamID='null'
    return int(start), int(end)

def download_pfam(pfamID):
    import urllib
    urllib.urlretrieve('http://pfam.xfam.org/family/%s/alignment/full' % pfamID, "%s%s.stockholm" % (directory, pfamID))


def download_pdb(pdbID):
    import urllib
    urllib.urlretrieve('http://www.rcsb.org/pdb/files/%s.pdb' % pdbID, "%s%s.pdb" % (directory, pdbID))

def stockholm2fasta(pfamID):
    from Bio import AlignIO
    # rewrite Stockholm alignment in FASTA format
    input_handle = open("%s%s.stockholm" % (directory, pfamID), "rU")
    output_handle = open("%s%s.fasta" % (directory, pfamID), "w")
    alignments = AlignIO.parse(input_handle, "stockholm")
    AlignIO.write(alignments, output_handle, "fasta")
    output_handle.close()
    input_handle.close()

def filter_fasta(gap_threshold, pfamID, pdbID, chain, seq, resnos):
    from Bio import AlignIO
    import numpy
    import subprocess
    #gap_threshold=0.25
    pfam_start, pfam_end = get_pfam_map(pdbID, chain)
    mapped_seq = seq[resnos.index(pfam_start):resnos.index(pfam_end)+1]

    #print mapped fasta file
    f = open('%s%s%s_pfam_mapped.fasta' % (directory, pdbID, chain), 'w')
    f.write('>%s:%s pdb mapped to pfam\n' % (pdbID, chain))
    f.write(mapped_seq)
    f.close()

    submit=("%s/dca-frustratometer/muscle3.8.31_i86linux64 -profile -in1 %s%s.fasta -in2 %s%s%s_pfam_mapped.fasta -out %s%s%s.fasta" % (basedir, directory,pfamID,directory,pdbID,chain,directory,pdbID,chain))

    #print(submit)

    process = subprocess.Popen(submit.split(), stdout=subprocess.PIPE)
    process.communicate()
    
    # Filter sequences based on gaps in input sequence and gap threshold
    alignment = AlignIO.read(open("%s%s%s.fasta" % (directory, pdbID, chain)), "fasta")
    targetseq=alignment[-1].seq
    targetname=alignment[-1].name
    if targetseq=='':
        print("targetseq not found")

    output_handle = open("%s%s_msa_filtered.fasta" % (directory, pdbID), "w")
    target_array = numpy.array([list(targetseq)], numpy.character)
    bools = target_array != '-'
    sequences_passed_threshold = 0
    for i, record in enumerate(alignment):
            record_array = numpy.array([list(record.seq)], numpy.character)
            aligned_sequence = record_array[bools]
            if float(sum(aligned_sequence=='-'))/len(aligned_sequence) < gap_threshold:
                    output_handle.write(">%s\n" % record.id + "".join(aligned_sequence).upper()+'\n')
                    sequences_passed_threshold += 1
    output_handle.close()

    fastaseq=''.join(target_array[bools]).upper()
    
    stat_output = open(stat_output_file_name, "w")
    stat_output.write("FASTA_alignments " + str(len(alignment)) + "\n")
    stat_output.write("Filtered_alignments " + str(sequences_passed_threshold) + "\n")
    stat_output.close()

    return fastaseq, sequences_passed_threshold

def calc_plm(pdbID):
    import matlab.engine
    eng = matlab.engine.start_matlab()

    #import StringIO
    #out = StringIO.StringIO()

    fastafile=("%s%s_msa_filtered.fasta" % (directory, pdbID))
    outputfile=("%soutputfile.%s" % (directory, pdbID))
    lambda_h=0.01
    lambda_J=0.01
    reweighting_threshold=0.1
    nr_of_cores=1
    outputDistribution=("%soutputDistribution.%s" % (directory, pdbID))
    outputMatrix=("%soutputMatrix.%s" % (directory, pdbID))

    eng.addpath('%s/dca-frustratometer/plm' % basedir, nargout=0)
    eng.addpath('%s/dca-frustratometer/plm/functions' % basedir , nargout=0)
    eng.addpath('%s/dca-frustratometer/plm/3rd_party_code/minFunc' % basedir, nargout=0)
    eng.plmDCA_symmetric_mod7(fastafile,outputfile,lambda_h,lambda_J,reweighting_threshold,nr_of_cores,
                              outputDistribution,outputMatrix, nargout=0)#, stdout=out )


def read_dca(pdbID):
    import scipy.io
    import numpy
    mat = scipy.io.loadmat("%soutputMatrix.%s" % (directory, pdbID))
    h=mat['h']
    J=mat['J']
    q=mat['h'].shape[1]
    N=mat['h'].shape[0]
    fields = numpy.zeros((N, q))
    couplings = numpy.swapaxes(J,1,2)
    fields=h
    return fields, couplings

def vector(p1, p2):
    return [p2[0]-p1[0], p2[1]-p1[1], p2[2]-p1[2]]

def vabs(a):
    from math import sqrt
    return sqrt(pow(a[0],2)+pow(a[1],2)+pow(a[2],2))

def calc_dis(p1, p2):
    v = vector(p1, p2)
    return vabs(v)

def three2one(prot):
    """ translate a protein sequence from 3 to 1 letter code"""

    code = {"GLY" : "G", "ALA" : "A", "LEU" : "L", "ILE" : "I",
            "ARG" : "R", "LYS" : "K", "MET" : "M", "CYS" : "C",
            "TYR" : "Y", "THR" : "T", "PRO" : "P", "SER" : "S",
            "TRP" : "W", "ASP" : "D", "GLU" : "E", "ASN" : "N",
            "GLN" : "Q", "PHE" : "F", "HIS" : "H", "VAL" : "V",
            "M3L" : "K", "MSE" : "M", "CAS" : "C" }

    newprot = ""
    for a in prot:
        newprot += code.get(a, "?")

    return newprot

def calc_distances(pdbID, chainID):
    import numpy
    from Bio.PDB.PDBParser import PDBParser
    import warnings
    from Bio.PDB.PDBExceptions import PDBConstructionWarning
    warnings.simplefilter('ignore', PDBConstructionWarning)

    p = PDBParser(PERMISSIVE=1)
    struct_id = pdbID

    filename = directory + struct_id + ".pdb"

    s = p.get_structure(struct_id, filename)
    chains = s[0].get_list()
    sequence = []
    dis = []
    all_res = []
    resnos = []

    count = 0

    for chain in chains:
        if chain.get_id() == chainID.upper():
            for resid in chain:
                is_regular_res = resid.has_id('CA') and resid.has_id('O')
                res_id = resid.get_id()[0]

                if count==0:
                    pdbstart=resid.get_id()[1]
                count=count + 1

                if (res_id==' ' or res_id=='H_MSE' or res_id=='H_M3L' or res_id=='H_CAS') and is_regular_res:
                    all_res.append(resid)
                    sequence.append(resid.get_resname())
                    resnos.append(resid.get_id()[1])

    number_of_pdb_amino_acids = len(all_res)
    native_distances = numpy.zeros((number_of_pdb_amino_acids, number_of_pdb_amino_acids))

    for i in range(0, len(all_res)):
        dis.append([]);
        ires = all_res[i]
        xyz_CAi = ires['CA'].get_coord()
        for j in range(0, len(all_res)):
            jres = all_res[j]
            xyz_CAj = jres['CA'].get_coord()
            r = calc_dis(xyz_CAi, xyz_CAj)
            dis[i].append(r);
            native_distances[i,j]=r

    sequence=three2one(sequence)
    return native_distances, sequence, resnos, pdbstart

def map_fasta_to_pdb(pdbseq, fastaseq):
    if(fastaseq in pdbseq):
        pdbstart=pdbseq.find(fastaseq)
        fastastart=0
    elif(pdbseq in fastaseq):
        fastastart=fastaseq.find(pdbseq)
        pdbstart=0
        
    else:
        import re
        from Bio import pairwise2
        alignments = pairwise2.align.globalxx(fastaseq, pdbseq)

        fastastart=re.search("[A-Z]", alignments[0][1]).start()
        pdbstart=re.search("[A-Z]", alignments[0][0]).start()

    N=min(len(fastaseq),len(pdbseq))

    pdb_indices=range(pdbstart,N+pdbstart)
    dca_indices=range(fastastart,N+fastastart)
    map_to_dca = dict(zip(pdb_indices, dca_indices))
    map_to_pdb = dict(zip(dca_indices, pdb_indices))

    return dca_indices, pdb_indices, map_to_dca, map_to_pdb, len(pdb_indices)

def evaluate_configurational_frustration(sequence, number_of_decoys, number_of_residues):
    include_fields = True
    native_totals = numpy.zeros((number_of_residues, number_of_residues))
    decoy_means = numpy.zeros((number_of_residues, number_of_residues))
    decoy_stds = numpy.zeros((number_of_residues, number_of_residues))
    configurational_frustration_indices = numpy.zeros((number_of_residues, number_of_residues))
    for i in range(number_of_residues):
        i_dca = dca_indices[i]
        i_pdb = pdb_indices[i]
        amino_acid_i = amino_acid_type_dictionary[sequence[i_dca]]
        for j in range(i+minimum_sequence_separation, number_of_residues):
            j_dca = dca_indices[j]
            j_pdb = pdb_indices[j]
            amino_acid_j = amino_acid_type_dictionary[sequence[j_dca]]
            native_distance = native_distances[i_pdb, j_pdb]
            if native_distance < max_distance_threshold:
                native_coupling = couplings[i_dca, amino_acid_i, j_dca, amino_acid_j]
                native_totals[i, j] = native_coupling
                if include_fields:
                    native_fields = fields[i_dca, amino_acid_i] + fields[j_dca, amino_acid_j]
                    native_totals[i, j] += native_fields
                # Perform decoy calculations
                decoy_fields = numpy.zeros(number_of_decoys)
                decoy_couplings = numpy.zeros(number_of_decoys)
                decoy_totals = numpy.zeros(number_of_decoys)
                for decoy_i in range(number_of_decoys):
                    i_sequence_index = i_dca
                    j_sequence_index = j_dca
                    
                    #randomize indices
                    decoy_distance = max_distance_threshold + 1.0
                    while decoy_distance > max_distance_threshold:
                        # Draw random indices based on number of pdb residues
                        i_sequence_index = random.choice(dca_indices)
                        j_sequence_index = random.choice(dca_indices)
                        # Choose the pdb indices based on those random indices
                        i_sequence_index_pdb = map_to_pdb[i_sequence_index]
                        j_sequence_index_pdb = map_to_pdb[j_sequence_index]
                        # Calculate the distance to make sure it is below the threshold
                        decoy_distance = native_distances[i_sequence_index_pdb, j_sequence_index_pdb]
                        
                    #randomize amino acid types
                    i_amino_acid_type = random.choice(dca_indices)
                    j_amino_acid_type = random.choice(dca_indices)
                    amino_acid_i = amino_acid_type_dictionary[sequence[i_amino_acid_type]]
                    amino_acid_j = amino_acid_type_dictionary[sequence[j_amino_acid_type]]
                    
                    decoy_totals[decoy_i] = couplings[i_sequence_index, amino_acid_i, j_sequence_index, amino_acid_j]
                    if include_fields:
                        decoy_totals[decoy_i] += fields[i_sequence_index, amino_acid_i] + fields[j_sequence_index, amino_acid_j]
                        
                decoy_means[i, j] = numpy.mean(decoy_totals)
                decoy_stds[i, j] = numpy.std(decoy_totals)
                configurational_frustration_indices[i, j] = (native_totals[i, j]-decoy_means[i, j])/decoy_stds[i, j]
                
    return native_totals, decoy_means, decoy_stds, configurational_frustration_indices

def evaluate_mutational_frustration(sequence, number_of_decoys, number_of_residues):
    native_totals = numpy.zeros((number_of_residues, number_of_residues))
    decoy_means = numpy.zeros((number_of_residues, number_of_residues))
    decoy_stds = numpy.zeros((number_of_residues, number_of_residues))
    mutational_frustration_indices = numpy.zeros((number_of_residues, number_of_residues))
    for i in range(number_of_residues):
        i_dca = dca_indices[i]
        i_pdb = pdb_indices[i]
        amino_acid_i = amino_acid_type_dictionary[sequence[i_dca]]
        for j in range(i+minimum_sequence_separation, number_of_residues):
            j_dca = dca_indices[j]
            j_pdb = pdb_indices[j]
            amino_acid_j = amino_acid_type_dictionary[sequence[j_dca]]
            native_distance = native_distances[i_pdb, j_pdb]
            if native_distance < max_distance_threshold:
                native_coupling = couplings[i_dca, amino_acid_i, j_dca, amino_acid_j]
                native_totals[i, j] = native_coupling
                if include_fields:
                    native_fields = fields[i_dca, amino_acid_i] + fields[j_dca, amino_acid_j]
                    native_totals[i, j] += native_fields
                # Perform decoy calculations
                decoy_fields = numpy.zeros(number_of_decoys)
                decoy_couplings = numpy.zeros(number_of_decoys)
                decoy_totals = numpy.zeros(number_of_decoys)
                
                for decoy_i in range(number_of_decoys):
                    i_sequence_index = i_dca
                    j_sequence_index = j_dca
                    
                    #get random amino acids for i and j
                    i_amino_acid_type = random.choice(dca_indices)
                    j_amino_acid_type = random.choice(dca_indices)
                    amino_acid_i = amino_acid_type_dictionary[sequence[i_amino_acid_type]]
                    amino_acid_j = amino_acid_type_dictionary[sequence[j_amino_acid_type]]
                    
                    #compute energy terms for (i,j) pair
                    decoy_totals[decoy_i] += couplings[i_sequence_index, amino_acid_i, j_sequence_index, amino_acid_j]
                    if include_fields:
                        decoy_totals[decoy_i] += fields[i_sequence_index, amino_acid_i] + fields[j_sequence_index, amino_acid_j]
                        
                    #in mutational mode, all (i,k) and (j,k) pairs also contribute to the decoy energy
                    #let's compute those pairs now
                    for k in range(number_of_residues):
                        amino_acid_k = amino_acid_type_dictionary[sequence[k]]
                        #make sure to exclude native interactions that we already counted
                        if (k==i_sequence_index or k==j_sequence_index):
                            continue
                            
                        rik = native_distances[pdb_indices[i], pdb_indices[k]]
                        if (rik < max_distance_threshold):
                            decoy_totals[decoy_i] += couplings[i_sequence_index, amino_acid_i, k, amino_acid_k]
                        rjk = native_distances[pdb_indices[j], pdb_indices[k]]
                        if (rik < max_distance_threshold):
                            decoy_totals[decoy_i] += couplings[j_sequence_index, amino_acid_j, k, amino_acid_k]
                            
                decoy_means[i, j] = numpy.mean(decoy_totals)
                decoy_stds[i, j] = numpy.std(decoy_totals)
                mutational_frustration_indices[i, j] = (native_totals[i, j]-decoy_means[i, j])/decoy_stds[i, j]
                
    return native_totals, decoy_means, decoy_stds, mutational_frustration_indices


def evaluate_single_residue_frustration(sequence, number_of_decoys, number_of_residues):
    native_totals = numpy.zeros(number_of_residues)
    decoy_means = numpy.zeros(number_of_residues)
    decoy_stds = numpy.zeros(number_of_residues)
    single_residue_frustration_indices = numpy.zeros(number_of_residues)
    for i in range(number_of_residues):
        i_dca = dca_indices[i]
        i_pdb = pdb_indices[i]
        amino_acid_i = amino_acid_type_dictionary[sequence[i_dca]]
        
        #compute native energy
        native_totals[i] = evaluate_singleres_native_energy(i_dca, i_pdb, amino_acid_i, sequence)
        
        #compute decoy energies
        decoy_totals = evaluate_singleres_decoy_energy(i_dca, i_pdb, amino_acid_i, sequence, number_of_decoys)
                            
        decoy_means[i] = numpy.mean(decoy_totals)
        decoy_stds[i] = numpy.std(decoy_totals)
        single_residue_frustration_indices[i] = (native_totals[i] - decoy_means[i])/decoy_stds[i]
                
    return native_totals, decoy_means, decoy_stds, single_residue_frustration_indices

def evaluate_singleres_decoy_energy(i_dca, ipdb, amino_acid_i, sequence, number_of_decoys):
    total_energy = 0
    decoy_totals = numpy.zeros(number_of_decoys)
    for decoy_i in range(number_of_decoys):
        #get random amino acids for i
        i_amino_acid_type = random.choice(dca_indices)
        amino_acid_i = amino_acid_type_dictionary[sequence[i_amino_acid_type]]
        
        #compute decoy_energy
        decoy_totals[decoy_i] = evaluate_singleres_native_energy(i_dca, ipdb, amino_acid_i, sequence)

    return decoy_totals

def evaluate_singleres_native_energy(i_dca, i_pdb, amino_acid_i, sequence):
    total_energy = 0
    for j in range(0, number_of_residues):
        j_dca = dca_indices[j]
        j_pdb = pdb_indices[j]
        if i_dca == j_dca:
            continue
        amino_acid_j = amino_acid_type_dictionary[sequence[j_dca]]
        native_distance = native_distances[i_pdb, j_pdb]
        if (native_distance < max_distance_threshold) and (abs(i_dca-j_dca)>minimum_sequence_separation):
            total_energy = couplings[i_dca, amino_acid_i, j_dca, amino_acid_j]
        if include_fields:
            total_energy += fields[i_dca, amino_acid_i]
            
    return total_energy

########MAIN PROGRAM##########

# Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument("pdbID", help="id of pdb")
parser.add_argument("chain", help="chain id of pdb")
args = parser.parse_args()

#Parameters for frustration calculations
randomize_amino_acid_types=True
randomize_indices=True
include_fields=True
compute_configurational_frustration=True
compute_mutational_frustration=True
compute_single_residue_frustration=True
n_decoys=5000
minimum_sequence_separation=2
max_distance_threshold=9.5

pdbID = args.pdbID#'5pti'
chain = args.chain#'a'

directory = ('%s/dca-frustratometer/automated/%s%s/' % (scratchdir, pdbID, chain))
if not os.path.exists(directory):
    os.makedirs(directory)
    
configurational_frustration_output_file_name=("%sconfigurational_frust_output.%s" % (directory, pdbID))
mutational_frustration_output_file_name=("%smutatational_frust_output.%s" % (directory, pdbID))
single_residue_frustration_output_file_name=("%ssing_frust_output.%s" % (directory, pdbID))
stat_output_file_name = ("%s%s.info" % (directory, pdbID))
status_file_name = ("%s%s.status" % (directory, pdbID))

status_output = open(status_file_name, "w")

#get the pfam and uniprot ID from the pdb info
pfamID = get_pfamID(pdbID, chain)
#uniprotID = get_uniprotID(pdbID, chain)

status_output.write("pdbID: " + pdbID + "\n")
print("pdbID: " + pdbID + "\n")
status_output.write("pfamID: " + pfamID + "\n")
print("pfamID: " + pfamID + "\n")
#status_output.write("uniprotID: " + uniprotID + "\n")
#print("uniprotID: " + uniprotID + "\n")

#download pdb file if it doesn't exist
if os.path.isfile(directory + pdbID + '.pdb') == False:
    download_pdb(pdbID)
    status_output.write("PDB downloaded" + "\n")
    print("PDB downloaded" + "\n")


#download pfam msa if it doesn't exist
if os.path.isfile(directory + pfamID + '.stockholm') == False:
    download_pfam(pfamID)
    status_output.write("Pfam MSA downloaded" + "\n")
    print("Pfam MSA downloaded" + "\n")

if os.path.isfile(directory + pfamID + '.fasta') == False:
    stockholm2fasta(pfamID)
    status_output.write("stockholm coverted" + "\n")
    print("stockholm coverted" + "\n")

#obtain information from PDB file (pairwise separation(native_distances) and pdb_sequences)
native_distances, pdb_sequence, resnos, pdb_start = calc_distances(pdbID, chain)

#filter sequences and get the target sequence and number of sequences remaining in the alignment
#set filter threshold
filter_threshold = 0.25
msa_query_sequence, n_sequences = filter_fasta(filter_threshold, pfamID, pdbID, chain, pdb_sequence, resnos)
status_output.write("FASTA filtered" + "\n")
print("FASTA filtered" + "\n")

#run the plm calculation which saves appropriate files
if os.path.isfile("%soutputMatrix.%s" % (directory, pdbID)) == False:
    calc_plm(pdbID)
    status_output.write("PLM computed" + "\n")
    print("PLM computed" + "\n")

#read in plm calculations fields and couplings
fields, couplings = read_dca(pdbID)

#map the fasta sequence on to the sequence of the pdb file
dca_indices, pdb_indices, map_to_dca,map_to_pdb, number_of_residues = map_fasta_to_pdb(pdb_sequence, msa_query_sequence)
status_output.write("parameter initialization complete" + "\n")
print("parameter initialization complete" + "\n")


#print some frustratometer statistics
stat_output = open(stat_output_file_name, "a")
stat_output.write("msa_query_sequence " + msa_query_sequence + "\n")
stat_output.write("msa_query_sequence_length " + str(len(msa_query_sequence)) + "\n")
stat_output.write("pdb_sequence " + pdb_sequence + "\n")
stat_output.write("pdb_sequence_length " + str(len(pdb_sequence)) + "\n")
stat_output.close()

# Define parameters
number_of_amino_acid_types = 21
amino_acid_type_dictionary = {
'-':0,
'A':1,
'C':2,
'D':3,
'E':4,
'F':5,
'G':6,
'H':7,
'I':8,
'K':9,
'L':10,
'M':11,
'N':12,
'P':13,
'Q':14,
'R':15,
'S':16,
'T':17,
'V':18,
'W':19,
'Y':20,
'X':0
}


# Compute and write out configurational frustration indices
target = open(("%s%s_configurational_frustration.tcl" % (directory, pdbID)), 'w')
atomselect=0
if compute_configurational_frustration:
    status_output.write("Computing configurational frustration...\n")
    print("Computing configurational frustration...\n")
    configurational_native_totals, configurational_decoy_means, configurational_decoy_stds, configurational_frustration_indices = evaluate_configurational_frustration(msa_query_sequence, n_decoys, number_of_residues)
    configurational_frustration_output_file = open(configurational_frustration_output_file_name, "w")
    for i in range(number_of_residues):
        for j in range(i+minimum_sequence_separation, number_of_residues):
            if native_distances[pdb_indices[i], pdb_indices[j]] < max_distance_threshold:
                configurational_frustration_output_file.write("%4d %s %4d %s %8.2f %8.2f %8.2f %8.2f \n" % (pdb_indices[i]+pdb_start,
                                                                                                     msa_query_sequence[dca_indices[i]],
                                                                                                     pdb_indices[j]+pdb_start,
                                                                                                     msa_query_sequence[dca_indices[j]],
                                                                                                     configurational_native_totals[i, j],
                                                                                                     configurational_decoy_means[i, j],
                                                                                                     configurational_decoy_stds[i, j],
                                                                                                     configurational_frustration_indices[i, j]))
                
                if (configurational_frustration_indices[i,j] > 0.78 or configurational_frustration_indices[i,j] < -1):
                    pdbi=pdb_indices[i]+pdb_start-1
                    pdbj=pdb_indices[j]+pdb_start-1
                    target.write("set sel%d [atomselect top \"resid %d and name CA\"]\n" % (pdbi, pdbi+1))
                    target.write("set sel%d [atomselect top \"resid %d and name CA\"]\n" % (pdbj, pdbj+1))
                    target.write("lassign [atomselect%d get {x y z}] pos1\n" % atomselect)
                    atomselect = atomselect + 1
                    target.write("lassign [atomselect%d get {x y z}] pos2\n" % atomselect)
                    atomselect = atomselect + 1
                    if configurational_frustration_indices[i,j] > 0.78:
                        target.write("draw color green\n")
                    else:
                        target.write("draw color red\n")
                    target.write("draw line $pos1 $pos2 style solid width 1\n")
    target.write("mol modselect 0 top \"all\"\n")
    target.write("mol modstyle 0 top newcartoon\n")
    target.write("mol modcolor 0 top colorid 15\n")
    configurational_frustration_output_file.close()
target.close()

# Compute and write out mutational frustration indices
target = open(("%s%s_mutational_frustration.tcl" % (directory, pdbID)), 'w')
atomselect=0
if compute_mutational_frustration:
    status_output.write("Computing mutational frustration...\n")
    print("Computing mutational frustration...\n")
    mutational_native_totals, mutational_decoy_means, mutational_decoy_stds, mutational_frustration_indices = evaluate_mutational_frustration(msa_query_sequence, n_decoys, number_of_residues)
    mutational_frustration_output_file = open(mutational_frustration_output_file_name, "w")
    for i in range(number_of_residues):
        for j in range(i+minimum_sequence_separation, number_of_residues):
            if native_distances[pdb_indices[i], pdb_indices[j]] < max_distance_threshold:
                mutational_frustration_output_file.write("%4d %s %4d %s %8.2f %8.2f %8.2f %8.2f \n" % (pdb_indices[i]+pdb_start,
                                                                                                     msa_query_sequence[dca_indices[i]],
                                                                                                     pdb_indices[j]+pdb_start,
                                                                                                     msa_query_sequence[dca_indices[j]],
                                                                                                     mutational_native_totals[i, j],
                                                                                                     mutational_decoy_means[i, j],
                                                                                                     mutational_decoy_stds[i, j],
                                                                                                     mutational_frustration_indices[i, j]))
                if (mutational_frustration_indices[i,j] > 0.78 or mutational_frustration_indices[i,j] < -1):
                    pdbi=pdb_indices[i]+pdb_start-1
                    pdbj=pdb_indices[j]+pdb_start-1
                    target.write("set sel%d [atomselect top \"resid %d and name CA\"]\n" % (pdbi, pdbi+1))
                    target.write("set sel%d [atomselect top \"resid %d and name CA\"]\n" % (pdbj, pdbj+1))
                    target.write("lassign [atomselect%d get {x y z}] pos1\n" % atomselect)
                    atomselect = atomselect + 1
                    target.write("lassign [atomselect%d get {x y z}] pos2\n" % atomselect)
                    atomselect = atomselect + 1
                    if mutational_frustration_indices[i,j] > 0.78:
                        target.write("draw color green\n")
                    else:
                        target.write("draw color red\n")
                    target.write("draw line $pos1 $pos2 style solid width 1\n")
    target.write("mol modselect 0 top \"all\"\n")
    target.write("mol modstyle 0 top newcartoon\n")
    target.write("mol modcolor 0 top colorid 15\n")
    mutational_frustration_output_file.close()
target.close()

# Compute and write out single residue frustration indices
if compute_single_residue_frustration:
    status_output.write("Computing single residue frustration...\n")
    print("Computing single residue frustration...\n")
    target = open(("%s%s_singleres_frustration.tcl" % (directory, pdbID)), 'w')
    atomselect = 0
    single_residue_native_totals, single_residue_decoy_means, single_residue_decoy_stds, single_residue_frustration_indices = evaluate_single_residue_frustration(msa_query_sequence, n_decoys, number_of_residues)
    single_residue_frustration_output_file = open(single_residue_frustration_output_file_name, "w")
    for i in range(number_of_residues):
        single_residue_frustration_output_file.write("%4d %s %8.2f %8.2f %8.2f %8.2f \n" % (pdb_indices[i]+pdb_start,
                                                                                            msa_query_sequence[dca_indices[i]],
                                                                                            single_residue_native_totals[i],
                                                                                            single_residue_decoy_means[i],
                                                                                            single_residue_decoy_stds[i],
                                                                                            single_residue_frustration_indices[i]))
        atomselect = atomselect + 1
        #temp = 0.5*abs(single_residue_frustration_indices[i])
        target.write("mol addrep 0\n")
        target.write("mol modselect %d 0 resid %d\n" % (atomselect, pdb_indices[i]+pdb_start))
        target.write("mol modstyle %d 0 VDW %f 12.000000\n" % (atomselect, 0.5*abs(single_residue_frustration_indices[i])))
        target.write("mol modmaterial %d 0 Transparent\n" % atomselect)
        if single_residue_frustration_indices[i] > 0:
            target.write("mol modcolor %d 0 ColorID 7\n" % atomselect)
        else:
            target.write("mol modcolor %d 0 ColorID 1\n" % atomselect)
    target.write("mol modselect 0 top \"all\"\n")
    target.write("mol modstyle 0 top newcartoon\n")
    target.write("mol modcolor 0 top colorid 15\n")
    target.close()
    single_residue_frustration_output_file.close()

status_output.close()
