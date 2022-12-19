import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
from Bio import SeqIO
import streamlit as st
from io import StringIO
from Bio import pairwise2
from Bio.Data import CodonTable


st.title('DNA Analysis')
uploaded_fasta_file = st.file_uploader('Upload file', type=['fasta', 'fas'])


def gc_content(seq):
    result = float(str(seq).count('G') + str(seq).count('C')) / len(seq) * 100
    return result


def at_content(seq):
    result = float(str(seq).count('A') + str(seq).count('T')) / len(seq) * 100
    return result


def count_dict(seq):
    dict = {}
    for base in seq:
        if base in dict:
            dict[base] += 1
        else:
            dict[base] = 0

    return dict


if uploaded_fasta_file is not None:
    stringio = StringIO(uploaded_fasta_file.getvalue().decode("utf-8"))
    fasta_file = SeqIO.parse(stringio, 'fasta')

    sequence_list = []

    for sequence in fasta_file:
        st.write('Sequence ID:', sequence.id)

        with st.expander('Show More'):
            st.write('Sequence:', sequence.seq)

        sequence_list.append(sequence)

    for sequence in sequence_list:
        with st.expander('Show Amino Acid Details'):
            Amino_Acid = sequence.translate()
            st.write('id:', Amino_Acid.id)
            st.write('annotations:', Amino_Acid.annotations)
            st.write('Sequence:', Amino_Acid.seq)

    if st.checkbox('GC Content'):
        for sequence in sequence_list:
            result = gc_content(sequence.seq)
            st.write('Results: ', result, '%')

    if st.checkbox('AT Content'):
        for sequence in sequence_list:
            result = at_content(sequence.seq)
            st.write('Results: ', result, '%')

    if st.checkbox('Visualization'):
        count = 0
        for sequence in sequence_list:
            count = count_dict(sequence.seq)

        st.write('Bases Count: ', count)
        st.set_option('deprecation.showPyplotGlobalUse', False)
        sns.barplot(x=["A", "G", "T", "C"], y=list(count.values()))
        plt.xlabel('Bases')
        plt.ylabel('Count')
        st.pyplot()

        plt.boxplot(list(count.values()))
        st.pyplot()

st.write('-------------------------------------------------')
st.title('Globel Alignment')
gap_penalty = -2
match_award = 1
mismatch_penalty = -1


def match_score(alpha, beta):
    if alpha == beta:
        return match_award
    else:
        return mismatch_penalty


def needleman_wunsch(seq1, seq2):
    n = len(seq1)
    m = len(seq2)

    score = np.zeros((m + 1, n + 1))
    for i in range(0, m + 1):
        score[i][0] = gap_penalty * i
    # Fill out first row
    for j in range(0, n + 1):
        score[0][j] = gap_penalty * j

    for i in range(1, m + 1):
        for j in range(1, n + 1):
            match = score[i - 1][j - 1] + match_score(seq1[j - 1], seq2[i - 1])
            delete = score[i - 1][j] + gap_penalty
            insert = score[i][j - 1] + gap_penalty
            score[i][j] = max(match, delete, insert)

    align1 = ""
    align2 = ""

    i = m
    j = n

    while i > 0 and j > 0:
        score_current = score[i][j]
        score_diagonal = score[i - 1][j - 1]
        score_up = score[i][j - 1]
        score_left = score[i - 1][j]

        if score_current == score_diagonal + match_score(seq1[j - 1], seq2[i - 1]):
            align1 += seq1[j - 1]
            align2 += seq2[i - 1]
            i -= 1
            j -= 1
        elif score_current == score_up + gap_penalty:
            align1 += seq1[j - 1]
            align2 += '-'
            j -= 1
        elif score_current == score_left + gap_penalty:
            align1 += '-'
            align2 += seq2[i - 1]
            i -= 1

    align1 = align1[::-1]
    align2 = align2[::-1]

    return align1, align2, score


seq1 = st.text_input('Enter the first sequence')
seq2 = st.text_input('Enter the second sequence')
if st.button('Align'):
    align1, align2, score = needleman_wunsch(seq1, seq2)
    st.write(align1)
    st.write(align2)
    st.write('score: ', score[len(seq2)][len(seq1)])

st.write('-----------------------------------------------------')


DNA_file = st.file_uploader('Upload DNA file', type=['fasta', 'fas', 'fna'])

if DNA_file is not None:
    DNA_stringio = StringIO(DNA_file.getvalue().decode("utf-8"))
    DNAsequence = SeqIO.parse(DNA_stringio, "fasta")
    for record in DNAsequence:
        Record = record.seq #Convert DNA into mRNA Sequence
        with st.expander('Show sequence'):
            st.write(Record)
        st.write('Size : ', len(Record))

    DNA_stringio = StringIO(DNA_file.getvalue().decode("utf-8"))
    DNAsequence = SeqIO.read(DNA_stringio, "fasta")
    st.write(DNAsequence)

    #Convert DNA into mRNA Sequence
    DNA = DNAsequence.seq
    mRNA = DNA.transcribe()   #Transcribe a DNA sequence into RNA.

    with st.expander('Show mRNA'):
        st.write(mRNA)

    st.write('mRNA Length:', len(mRNA))



    st.subheader('Table of Codons')
    with st.expander('Show Table'):
        st.text(CodonTable.unambiguous_rna_by_name['Standard'])

    # Obtain Amino Acid Sequence from mRNA
    Amino_Acid = mRNA.translate(table=1, cds=False)
    with st.expander('Show Amino Acid'):
        st.write('Amino Acid : ', Amino_Acid)

    st.write("Length of Protein : ",len(Amino_Acid))
    st.write("Length of Original mRNA : ",len(mRNA))


    Proteins = Amino_Acid.split('*') # * is translated stop codon
    for i in Proteins[:]:
        if len(i) < 20:
            Proteins.remove(i)

    st.subheader('Proteins that have length more that 20 AA')
    with st.expander('Show Proteins'):
        st.write(Proteins)


    poi_list = []
    MW_list = []
    from Bio.SeqUtils import ProtParam
    with st.expander('Show Proteins'):
        for record in Proteins[:]:
            # print("\n")
            X = ProtParam.ProteinAnalysis(str(record))
            POI = X.count_amino_acids()
            poi_list.append(POI)
            MW = X.molecular_weight()
            MW_list.append(MW)
            st.write("Protein of Interest = ", POI)
            st.write("Amino acids percent = ", str(X.get_amino_acids_percent()))
            st.write("Molecular weight = ", MW)
            st.write("Aromaticity = ", X.aromaticity())
            st.write("Flexibility = ", X.flexibility())
            st.write("Isoelectric point = ", X.isoelectric_point())
            st.write("Secondary structure fraction = ", X.secondary_structure_fraction())

    poi_list = poi_list[48]
    st.set_option('deprecation.showPyplotGlobalUse', False)
    plt.bar(poi_list.keys(), list(poi_list.values()), align='center')
    plt.xlabel('Amino Acids')
    plt.ylabel('Counts')
    st.pyplot()

    SARS_file_uploader = st.file_uploader('Upload SARS file', type=['fasta', 'fas'])
    MERS_file_uploader = st.file_uploader('Upload MERS file', type=['fasta', 'fas'])
    COV2_file_uploader = st.file_uploader('Upload COV2 file', type=['fasta', 'fas'])

    if SARS_file_uploader is not None and MERS_file_uploader is not None and COV2_file_uploader is not None:

        SARS_stringio = StringIO(SARS_file_uploader.getvalue().decode("utf-8"))
        MERS_stringio = StringIO(MERS_file_uploader.getvalue().decode("utf-8"))
        COV2_stringio = StringIO(COV2_file_uploader.getvalue().decode("utf-8"))

        SARS_fasta_file = SeqIO.read(SARS_stringio, 'fasta')
        MERS_fasta_file = SeqIO.read(MERS_stringio, 'fasta')
        COV2_fasta_file = SeqIO.read(COV2_stringio, 'fasta')

        # SARS = SeqIO.read("C:\\Users\\Steven20367691\\Desktop\\sars.fasta", "fasta")
        # MERS = SeqIO.read("C:\\Users\\Steven20367691\\Desktop\\mers.fasta", "fasta")
        # COV2 = SeqIO.read("C:\\Users\\Steven20367691\\Desktop\\cov2.fasta", "fasta")

        st.write('Sequence Lengths:')
        st.write('SARS:', len(SARS_fasta_file.seq))
        st.write('COV2:', len(COV2_fasta_file.seq))
        st.write('MERS:', len(MERS_fasta_file.seq))

        st.write(SARS_fasta_file)

        # Alignments using pairwise2 alghoritm
        SARS_Sample = pairwise2.align.globalxx(SARS_fasta_file.seq, DNAsequence.seq, one_alignment_only=True, score_only=True)
        st.write('SARS/DNAsequence Similarity (%):', SARS_Sample / len(SARS_fasta_file.seq) * 100)

        COV_Sample = pairwise2.align.globalxx(COV2_fasta_file.seq, DNAsequence.seq, one_alignment_only=True, score_only=True)
        st.write('MERS/DNAsequence Similarity (%):', COV_Sample / len(MERS_fasta_file.seq) * 100)

        MERS_SAMPLE = pairwise2.align.globalxx(MERS_fasta_file.seq, DNAsequence.seq, one_alignment_only=True, score_only=True)
        st.write('MERS/DNAsequence Similarity (%):', MERS_SAMPLE / len(SARS_fasta_file.seq) * 100)

        # Plot the data
        X = ['SARS/COV2', 'MERS/COV2', 'MERS/SARS']
        Y = [SARS_Sample/ len(SARS_fasta_file.seq) * 100, COV_Sample/ len(MERS_fasta_file.seq)*100, MERS_SAMPLE/len(SARS_fasta_file.seq)*100]

        plt.title('Sequence identity (%)')
        st.set_option('deprecation.showPyplotGlobalUse', False)
        plt.bar(X, Y)
        st.pyplot()
