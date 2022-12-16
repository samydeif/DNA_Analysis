import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
from Bio import SeqIO
import streamlit as st
from io import StringIO

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

