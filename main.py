import streamlit as st
import pandas as pd
import numpy as np
from IPython.core.display import HTML
from itertools import product
from collections import deque

st.title("Enhanced DPSA")


def find_match(a, b, match, mismatch) -> list:
    """
    @param a: str
    @param b: str
    @param match: int
    @param mismatch: list(int)
    @return: tuple
    """

    if a == b:
        return match

    if ((a == "A" or a == "G") and (b == "G" or b == "A")) or (
        (a == "T" or a == "C") and (b == "C" or b == "T")
    ):
        return mismatch[0]

    return mismatch[1]


def needleman_wunsch(x: str, y: str, match: int, mismatch: list, gap: int) -> tuple:

    N, M = len(x), len(y)
    # s = lambda a, b: match if a == b else mismatch
    s = lambda a, b: find_match(a, b, match, mismatch)

    DIAG = -1, -1
    LEFT = -1, 0
    UP = 0, -1

    # Create tables F and Ptr
    F = {}
    Ptr = {}

    F[-1, -1] = 0
    for i in range(N):
        F[i, -1] = -i
    for j in range(M):
        F[-1, j] = -j

    option_Ptr = DIAG, LEFT, UP
    for i, j in product(range(N), range(M)):
        option_F = (
            F[i - 1, j - 1] + s(x[i], y[j]),
            F[i - 1, j] + gap,
            F[i, j - 1] + gap,
        )
        F[i, j], Ptr[i, j] = max(zip(option_F, option_Ptr))

    # Work backwards from (N - 1, M - 1) to (0, 0)
    # to find the best alignment.
    alignment = deque()
    i, j = N - 1, M - 1
    while i >= 0 and j >= 0:
        direction = Ptr[i, j]
        if direction == DIAG:
            element = i, j
        elif direction == LEFT:
            element = i, None
        elif direction == UP:
            element = None, j
        alignment.appendleft(element)
        di, dj = direction
        i, j = i + di, j + dj
    while i >= 0:
        alignment.appendleft((i, None))
        i -= 1
    while j >= 0:
        alignment.appendleft((None, j))
        j -= 1

    return list(alignment)


def print_alignment(x, y, alignment):
    st.write("".join("-" if i is None else x[i] for i, _ in alignment))
    st.write("".join("-" if j is None else y[j] for _, j in alignment))


def pretty_table(data_array, row_labels, col_labels):
    df = pd.DataFrame(data_array, columns=col_labels, index=row_labels)
    table_html = df.to_html()
    return HTML(table_html)


def calculate(seqtitle, seq1, seq2):
    up_arrow = "\u2191"
    left_arrow = "\u2190"
    up_left_arrow = "\u2196"

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
    st.write("\n")
    print_alignment(
        seq1,
        seq2,
        needleman_wunsch(
            seq1,
            seq2,
            match_bonus,
            [mismatch_penalty_transition, mismatch_penalty_transversion],
            gap_penalty,
        ),
    )


with st.sidebar.form("Input"):
    seqtitle = st.text_input("Enter sequence pair: ")
    seq1 = st.text_input("Enter sequence 1: ")
    seq2 = st.text_input("Enter sequence 2: ")
    btnProceed = st.form_submit_button("Calculate")

if btnProceed:
    calculate(seqtitle, seq1, seq2)
