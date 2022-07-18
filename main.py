import streamlit as st
import pandas as pd
import numpy as np
from IPython.core.display import HTML

st.title("Enhanced DPSA")


def pretty_table(data_array, row_labels, col_labels):
    df = pd.DataFrame(data_array, columns=col_labels, index=row_labels)
    table_html = df.to_html()
    return HTML(table_html)


def calculate(seqtitle, seq1, seq2):
    up_arrow = "\u2191"
    left_arrow = "\u2190"
    up_left_arrow = "\u2196"

    # seqtitle = input("Enter sequence pair: ")
    # seq1 = input("Enter sequence 1: ")
    # seq2 = input("Enter sequence 2: ")

    n_rows = len(seq1) + 1
    n_cols = len(seq2) + 1
    row_labels = [label for label in "-" + seq1]
    col_labels = [label for label in "-" + seq2]

    scoring_array = np.full((n_rows, n_cols), 0)
    traceback_array = np.full((n_rows, n_cols), "-")

    arrow = "-"
    gap_penalty = -3
    match_bonus = 2
    mismatch_penalty_transition = -1  # A/G -> G/A || T/C -> C/T
    mismatch_penalty_transversion = -2  # A/G -> C/T || C/T -> A/G

    for row in range(n_rows):
        for col in range(n_cols):
            if row == 0 and col == 0:
                score = 0
                arrow = "-"
            elif row == 0:
                previous_score = scoring_array[row, col - 1]
                score = previous_score + gap_penalty
                arrow = left_arrow
            elif col == 0:
                previous_score = scoring_array[row - 1, col]
                score = previous_score + gap_penalty
                arrow = up_arrow

            else:
                cell_to_the_left = scoring_array[row, col - 1]
                from_left_score = cell_to_the_left + gap_penalty
                above_cell = scoring_array[row - 1, col]
                from_above_score = above_cell + gap_penalty
                diagonal_left_cell = scoring_array[row - 1, col - 1]

                if seq1[row - 1] == seq2[col - 1]:
                    diagonal_left_score = diagonal_left_cell + match_bonus
                elif (
                    (seq1[row - 1] == "A" and seq2[col - 1] == "G")
                    or (seq1[row - 1] == "G" and seq2[col - 1] == "A")
                    or (seq1[row - 1] == "T" and seq2[col - 1] == "C")
                    or (seq1[row - 1] == "C" and seq2[col - 1] == "T")
                ):
                    diagonal_left_score = (
                        diagonal_left_cell + mismatch_penalty_transition
                    )
                else:
                    diagonal_left_score = (
                        diagonal_left_cell + mismatch_penalty_transversion
                    )

                score = max([from_left_score, from_above_score, diagonal_left_score])

                if score == from_left_score:
                    arrow = left_arrow
                elif score == from_above_score:
                    arrow = up_arrow
                elif score == diagonal_left_score:
                    arrow = up_left_arrow

            traceback_array[row, col] = arrow
            scoring_array[row, col] = score

    st.write(f"Pair: {seqtitle}")
    st.write(pretty_table(scoring_array, row_labels, col_labels))
    st.write("\n")
    st.write(pretty_table(traceback_array, row_labels, col_labels))


with st.sidebar.form("Input"):
    seqtitle = st.text_input("Enter sequence pair: ")
    seq1 = st.text_input("Enter sequence 1: ")
    seq2 = st.text_input("Enter sequence 2: ")
    btnProceed = st.form_submit_button("Calculate")

if btnProceed:
    calculate(seqtitle, seq1, seq2)
