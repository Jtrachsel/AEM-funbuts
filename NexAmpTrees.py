from Bio import AlignIO, SeqIO, Seq
from Bio.Align import MultipleSeqAlignment
from Bio.Seq import MutableSeq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils.CheckSum import seguid
from Bio.Seq import Seq
from Bio.Data.IUPACData import ambiguous_dna_values
import pandas as pd
import re
from Bio.Align.Applications import ClustalOmegaCommandline
from Bio.Phylo.Applications import RaxmlCommandline
import os




def length_filter(sequences, minlength=None, maxlength=None):
    lengths = []
    for seq in sequences:
        lengths.append(len(seq.seq))

    if minlength is None:
        minlength = min(lengths)
    if maxlength is None:
        maxlength = max(lengths)
    goodseqs_min = []
    goodseqs_max = []
    for seq in sequences:
        if len(seq.seq) <= maxlength:
            goodseqs_max.append(seq)
        else:
            print("sequence {} is longer than the maximum length and has been removed".format(seq.id))
    for seq in goodseqs_max:
        if len(seq.seq) >= minlength:
            goodseqs_min.append(seq)
        else:
            print("sequence {} is shorter than the minimum length and has been removed".format(seq.id))

    return goodseqs_min


def alignment_gap_standardize(alignment):
    """ This replaces '.' positions in the alignment with '-' """
    changed_alignment = MultipleSeqAlignment([])
    for entry in alignment:
        new_entry = []
        for position in entry.seq:

            if position == ".":
                new_position = "-"
                new_entry.append(new_position)
            else:
                new_position = position
                new_entry.append(new_position)

        new_entry = "".join(new_entry)

        new_seq = Seq(new_entry)
        entry.seq = new_seq
        changed_alignment.append(entry)
    return changed_alignment


def codon_align(protein_alignment, dna_sequences):

    codon_positions = []

    for seq in protein_alignment:
        codons = []
        codon_count = 0
        for position in seq.seq:
            if position != "-":
                codons.append(codon_count)
                codon_count += 3
            else:
                codon_count += 3
        codon_positions.append(codons)

    DNA_sequence_index = 0
    Codon_list_index = 0

    DNA_align_length = len(protein_alignment[1].seq) * 3
    Codon_DNA_alignment = []

    for x in codon_positions:
        blank_alignment = "-" * DNA_align_length
        blank_alignment = MutableSeq(blank_alignment)
        DNA_align = blank_alignment
        DNA_position_index = 0
        for y in x:
            DNA_align[y:y+3] = dna_sequences[DNA_sequence_index].seq[DNA_position_index:(DNA_position_index+3)]
            DNA_position_index += 3

        DNA_align = DNA_align.toseq()
        DNA_align = SeqRecord(DNA_align)
        DNA_align.id = dna_sequences[DNA_sequence_index].id
        DNA_align.name = dna_sequences[DNA_sequence_index].name
        DNA_align.description = dna_sequences[DNA_sequence_index].description
        DNA_align.dbxrefs = dna_sequences[DNA_sequence_index].dbxrefs
        DNA_sequence_index += 1
        Codon_DNA_alignment.append(DNA_align)
        Codon_list_index += 1
    return Codon_DNA_alignment


def unique_seqs(sequences):
    """returns a list of SeqRecord objects with redundant sequences removed"""
    unique_records = []
    checksum_container = []
    for seq in sequences:
        checksum = seguid(seq.seq)
        if checksum not in checksum_container:
            checksum_container.append(checksum)
            unique_records.append(seq)
    return unique_records


def unique_ids(sequences):
    """returns a list of SeqRecord objects with redundant ids renamed"""
    unique_records = []
    checksum_container = []
    redundant_id_count = 0
    for seq in sequences:
        checksum = seguid(seq.id)
        if checksum not in checksum_container:
            checksum_container.append(checksum)
            unique_records.append(seq)
        else:
            print("repeated id detected, adding '.{}' suffix".format(redundant_id_count))
            seq.id = "{}.{}".format(seq.id, redundant_id_count)
            unique_records.append(seq)
            redundant_id_count += 1
    return unique_records


def make_protein_record(nuc_record):
    """Returns a new SeqRecord with the translated sequence (default table)."""

    return SeqRecord(seq=nuc_record.seq.translate(cds=False), id=nuc_record.id, description="translation using def table")


def alignment_screen(protein_alignment, protein_seqs, dna_seqs, cutoff):
    aln = protein_alignment
    alignment_length = (len(aln[0]))
    alignment_iterator = range(0, alignment_length)  # need the +1 because python doesnt include the last index

    iterator_count = 0
    list_of_percent_gaps = []
    for x in alignment_iterator:
        column = aln[:, iterator_count]
        number_of_gaps = column.count('-')
        percent_gaps = number_of_gaps/len(column)
        list_of_percent_gaps.append(percent_gaps)
        iterator_count += 1
    print(list_of_percent_gaps)

    counter = 0
    good_columns = []
    bad_columns = []
    for entry in list_of_percent_gaps:
        if entry < cutoff:
            good_columns.append(counter)
            counter += 1
        else:
            bad_columns.append(counter)
            counter += 1

    prot_DNA_dict = {}
    for sequence in protein_alignment:
        position_counter = 0
        residue_counter = 0
        prot_DNA_dict[sequence.id] = {}
        for position in sequence.seq:
            if position != "-":
                prot_DNA_dict[sequence.id][position_counter] = residue_counter
                residue_counter += 1
            position_counter += 1

    residues_removed = {}
    for entry in prot_DNA_dict:
        # print(entry)
        residues_removed[entry] = []
        for position in prot_DNA_dict[entry]:
            if position in bad_columns:
                residues_removed[entry].append(prot_DNA_dict[entry][position])

    print(residues_removed)
    nucleotides_removed = {}
    for entry in residues_removed:
        nucleotides_removed[entry] = []
        for y in residues_removed[entry]:
            y = [y*3, (y*3) + 1, (y*3) + 2]
            nucleotides_removed[entry].append(y)
        nucleotides_removed[entry] = sum(nucleotides_removed[entry], [])  # flattens the lists of lists into just lists
    for entry in nucleotides_removed:
        print("sequence ID: {}\nNucleotides removed: {}".format(entry, nucleotides_removed[entry]))

    prot_dict = {}
    dna_dict = {}

    for sequence in protein_seqs:
        prot_dict[sequence.id] = sequence.seq

    for sequence in dna_seqs:
        dna_dict[sequence.id] = sequence.seq

    edited_prots = []

    for seq in protein_seqs:
        blank_seq = []
        position_counter = 0
        # print(len(seq.seq))
        print("id = {}".format(seq.id))
        for position in seq.seq:
            if position_counter in residues_removed[seq.id]:

                print("RESIDUE #{} REMOVED".format(position_counter))
                position_counter += 1
            else:
                position_counter += 1
                blank_seq.append(position)
        new_seq = "".join(blank_seq)
        seq.seq = Seq(new_seq)
        edited_prots.append(seq)

    edited_DNAs = []

    for seq in dna_seqs:
        blank_seq = []
        position_counter = 0
        for position in seq.seq:
            if position_counter in nucleotides_removed[seq.id]:
                position_counter += 1
            else:
                position_counter += 1
                blank_seq.append(position)
        new_seq = "".join(blank_seq)
        seq.seq = Seq.Seq(new_seq)
        edited_DNAs.append(seq)

    return edited_prots, edited_DNAs


def primer_coverage(FWDprimer, REVprimer, FWDregion, REVregion):
    """ Returns a pandas dataframe with the number of mismatches for each primer in their corresponding region.
    It also attempts to pull metadata from the sequences (organism, and gene definition) this works if the sequences
    are from fungene.  Row names are sequence IDs. """
    number_mismatches_FWD = {}

    for seq in FWDregion:
        number_mismatches_FWD[seq.id] = 0
        for count, position in enumerate(seq.seq.upper()):
            #print(count, position)
            #print(ambiguous_dna_values[fwdprimer[count]])
            if position not in ambiguous_dna_values[FWDprimer[count]]:  # because you changed name fwdprimer obj change this
                number_mismatches_FWD[seq.id] += 1
    #print(number_mismatches_FWD)

    #print(REVprimerregion)
    REVprimer = REVprimer.reverse_complement()   # because you changed the name of the primer object you need to change this

    number_mismatches_REV = {}

    for seq in REVregion:
        number_mismatches_REV[seq.id] = 0
        for count, position in enumerate(seq.seq.upper()):
            #print(count, position)
            #print(ambiguous_dna_values[revprimer[count]])
            if position not in ambiguous_dna_values[REVprimer[count]]:
                number_mismatches_REV[seq.id] += 1
    print(number_mismatches_REV)

    no_FWD_mismatches = []
    one_FWD_mismatch = []
    more_than_one_FWD_mismatch = []
    for entry in number_mismatches_FWD:
        if number_mismatches_FWD[entry] == 0:
            no_FWD_mismatches.append(entry)
        if number_mismatches_FWD[entry] == 1:
            one_FWD_mismatch.append(entry)
        else:
            more_than_one_FWD_mismatch.append(entry)

    no_REV_mismatches = []
    one_REV_mismatch = []
    more_than_one_REV_mismatch = []
    for entry in number_mismatches_REV:
        if number_mismatches_REV[entry] == 0:
            no_REV_mismatches.append(entry)
        if number_mismatches_REV[entry] == 1:
            one_REV_mismatch.append(entry)
        else:
            more_than_one_REV_mismatch.append(entry)

    #print(no_FWD_mismatches)
    #print(no_REV_mismatches)
    perfect_FWDandREV = []
    perfect_FWD_oneREV = []
    one_FWD_perfect_REV = []
    one_FWD_one_REV = []

    for entry in no_FWD_mismatches:
        if entry in no_REV_mismatches:
            perfect_FWDandREV.append(entry)
        if entry in one_REV_mismatch:
            perfect_FWD_oneREV.append(entry)
    for entry in no_REV_mismatches:
        if entry in one_FWD_mismatch:
            one_FWD_perfect_REV.append(entry)

    for entry in one_FWD_mismatch:
        if entry in one_REV_mismatch:
            one_FWD_one_REV.append(entry)

    number_perfect = len(perfect_FWDandREV)
    number_one_fwd = len(one_FWD_perfect_REV)
    number_one_rev = len(perfect_FWD_oneREV)
    number_one_fwdandrev = len(one_FWD_one_REV)
    total_input = len(number_mismatches_FWD)

    # coverage_of_input =

    print("number of seqs with perfect primer coverage: {}".format(len(perfect_FWDandREV)))
    print("number of seqs with one FWD mismatch only: {}".format(len(one_FWD_perfect_REV)))
    print("number of seqs with one REV mismatch only: {}".format(len(perfect_FWD_oneREV)))
    print("number of seqs with one FWD and one REV mismatch: {}".format(len(one_FWD_one_REV)))
    print("total number of seqs that are likely to amplify: {}".format(len(perfect_FWDandREV)+len(one_FWD_perfect_REV)+len(perfect_FWD_oneREV)+len(one_FWD_one_REV)))

    # print(number_mismatches_FWD)
    mismatches_fwd_rev = {}
    for key in number_mismatches_FWD:
        #print("ID: {}".format(key))
        #print("number of FWD mismatches: {}".format(number_mismatches_FWD[key]))
        #print("number of REV mismatches: {}".format(number_mismatches_REV[key]))
        mismatches_fwd_rev[key] = {}
        mismatches_fwd_rev[key]["FWD_mismatch"] = number_mismatches_FWD[key]
        mismatches_fwd_rev[key]["REV_mismatch"] = number_mismatches_REV[key]
    #print(mismatches_fwd_rev)
    #print(mismatches_fwd_rev.values())

    mismatches_df = pd.DataFrame.from_dict(mismatches_fwd_rev,'index')
    print(mismatches_df)

    # making the metadata whatnot

    metadata = {}
    org_regex = re.compile(r"organism=.*,")
    def_regex = re.compile(r"definition=.*")
    for seq in FWDregion:
        #print(seq.description)
        metadata[seq.id] = {}
        organism = re.search(org_regex, seq.description).group(0)
        organism = (organism[9:(len(organism)-1)])  #gets rid of the "organism=" and "," that my regex pulled out
        definition = re.search(def_regex, seq.description).group(0)
        definition = (definition[11:(len(definition))])  #gets rid of the "definition=" that my regex pulled out
        #print(organism)
        #print(definition)
        metadata[seq.id]["organism"] = organism
        metadata[seq.id]["definition"] = definition

    metadata_df = pd.DataFrame.from_dict(metadata, 'index')
    #print(metadata_df)

    final_df = pd.concat([mismatches_df, metadata_df], axis=1)
    print(final_df)
    return final_df


def DNA_align_screen(DNAalignment, cutoff=0.75):
    aln = DNAalignment
    alignment_length = (len(aln[0]))
    alignment_iterator = range(0, alignment_length)

    iterator_count = 0
    list_of_percent_gaps = []
    for x in alignment_iterator:
        column = aln[:, iterator_count]
        number_of_gaps = column.count('-')
        percent_gaps = number_of_gaps/len(column)
        list_of_percent_gaps.append(percent_gaps)
        iterator_count += 1
    print(list_of_percent_gaps)

    counter = 0
    good_columns = []
    bad_columns = []
    for entry in list_of_percent_gaps:
        if entry < cutoff:
            good_columns.append(counter)
            counter += 1
        else:
            bad_columns.append(counter)
            counter += 1
    edited_DNAs = []

    for seq in aln:
        blank_seq = []
        position_counter = 0
        for position in seq.seq:
            if position_counter in bad_columns:
                position_counter += 1
            else:
                position_counter += 1
                blank_seq.append(position)
        new_seq = "".join(blank_seq)
        seq.seq = Seq(new_seq)
        edited_DNAs.append(seq)
        edited_DNAs = MultipleSeqAlignment(edited_DNAs)
    return edited_DNAs



############################
































###########################


metadata = {}
OTU = {}
reference = {}
functCONF = {}
functOTHER = {}
functUNKOWN = {}



#####  REFERENCE SECTION #####
os.chdir('/home/julian/NexAmp/framecorr/figglin/coproadd/')
allseqs = []

goodseqs = []

buts = SeqIO.parse("finalbutsforT.fasta", "fasta")


ref_seqs = []
for seq in buts:
    new_id = seq.id.replace(".", "_")
    seq.id = new_id
    print(seq.description)
    if "function=confirmed" in seq.description:
        ref_seqs.append(seq)
        functCONF[seq.id] = 'TRUE'
        OTU[seq.id] = 'FALSE'
        functOTHER[seq.id] = 'FALSE'
        functUNKOWN[seq.id] = 'FALSE'
        reference[seq.id] = 'TRUE'
        allseqs.append(seq)
    if "function=other" in seq.description:
        ref_seqs.append(seq)
        functOTHER[seq.id] = 'TRUE'
        OTU[seq.id] = 'FALSE'
        reference[seq.id] = 'TRUE'
        functUNKOWN[seq.id] = 'FALSE'
        functCONF[seq.id] = 'FALSE'
        allseqs.append(seq)


SeqIO.write(ref_seqs, "ref_seqs.fasta", "fasta")

align = AlignIO.read("ref_seqs.align", "fasta")

for seq in allseqs:
    if "Yersinia" in seq.id:
        allseqs.remove(seq)
    if "aminobutyricum" in seq.id:
        allseqs.remove(seq)
    if seq.id == "Clostridium_kluyveri_DSM_555":
        allseqs.remove(seq)
    if "tetani" in seq.id:
        allseqs.remove(seq)

for seq in allseqs:
    if seq.id == "Clostridium_kluyveri_DSM_555":
        allseqs.remove(seq)


for seq in allseqs:
    print(seq.id)

AlignIO.write(align, "ref_seqs.phy", 'phylip-relaxed')

raxml_cline = RaxmlCommandline(sequences="ref_seqs.phy", name='reftreeO2',
                               num_replicates='autoMRE', parsimony_seed='12345', algorithm='a', model='GTRGAMMA',
                               rapid_bootstrap_seed='12345', outgroup='Clostridium_kluyveri_DSM_555')
out, err = raxml_cline()

print(out)
print(err)






# this section should make the OTUseqs alignment and output it to fasta format
'''
This section should do the following:
    1) pull in representative sequences
    2) make the names be OTUXXX
    3) write the new corrected repseqs to a fasta file
    4) align the fasta file with clustalo

I NEED REFERENCE SEQUENCES IN THIS ALIGNMENT TOO
'''




repseqs = SeqIO.parse("/home/julian/NexAmp/framecorr/DvR.QCD.0.03reps.fasta", 'fasta')
repseqs = list(repseqs)
OTUnames = ["{0:03}".format(i) for i in range(93)]

OTUnames = OTUnames[1:len(OTUnames)]

OTUnames2 = []

for entry in OTUnames:
    entry = "Otu{}".format(entry)
    print(entry)
    OTUnames2.append(entry)


count = 0
for seq in repseqs:
    seq.id = OTUnames2[count]
    count += 1
    seq.description = ''
    allseqs.append(seq)
    OTU[seq.id] = 'TRUE'
    reference[seq.id] = 'FALSE'
    functCONF[seq.id] = 'FALSE'
    functUNKOWN[seq.id] = 'TRUE'
    functOTHER[seq.id] = 'FALSE'

### REMOVE GAPS FUNCTION ###
nogaps = []
for seq in repseqs:
    new_seq = []
    for position in seq.seq:
        if position == '-':
            continue
        elif position != '-':
            new_seq.append(position)
    new_seq = ''.join(new_seq)
    new_seq = Seq(new_seq)
    seq.seq = new_seq


SeqIO.write(allseqs, "OTUreps.fasta", 'fasta')

repcline = ClustalOmegaCommandline(infile="OTUreps.fasta", outfile="OTUreps.align.fasta", auto=True)
repcline()
    #result of this section is a fasta alignment of the representative sequences for each OTU

########### END OTU section  ################




reference_df = pd.DataFrame.from_dict(reference, 'index')
OTU_df = pd.DataFrame.from_dict(OTU, 'index')
functCONF_df = pd.DataFrame.from_dict(functCONF, 'index')
functOTHER_df = pd.DataFrame.from_dict(functOTHER, 'index')
functUNKOWN_df = pd.DataFrame.from_dict(functUNKOWN, 'index')



final_df = pd.merge(reference_df, OTU_df, left_index=True, right_index=True, how='outer')
final_df = pd.merge(final_df, functCONF_df, left_index=True, right_index=True, how='outer')
final_df = pd.merge(final_df, functOTHER_df, left_index=True, right_index=True, how='outer')
final_df = pd.merge(final_df, functUNKOWN_df, left_index=True, right_index=True, how='outer')


# this dataframe/table has issues... the headers are not correct...
final_df.to_csv("reftreeMETA.csv", sep="\t")






repref = AlignIO.read("OTUreps.align.fasta", 'fasta')
screen_repref = DNA_align_screen(repref)
AlignIO.write(screen_repref, 'screen_repref.fasta', 'fasta')
repref_cline = ClustalOmegaCommandline(infile="screen_repref.fasta", outfile="screen_repref.align.fasta", auto=True)
print("aligning screened sequences with ClustalO with the command: \n {}".format(repref_cline))
repref_cline()
shorts = AlignIO.read("screen_repref.align.fasta", 'fasta')
AlignIO.write(shorts, "shorts.phy", "phylip-relaxed")


shorts_cline = RaxmlCommandline(sequences="shorts.phy", name='shorts',
                               num_replicates='autoMRE', parsimony_seed='12345', algorithm='a', model='GTRGAMMA',
                               rapid_bootstrap_seed='12345', outgroup='Anaerostipes_caccae_strain_L1-92_')
print("Building tree with RAxML using the following settings: \n {}".format(shorts_cline))
out, err = shorts_cline()
print(out)
print(err)
