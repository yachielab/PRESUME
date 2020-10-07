__version__ = "1.0.0"

__author__ = "\
            Keito Watano <watano.k10.yachielab@gmail.com>,\
            Naoki Konno <naoki@bs.s.u-tokyo.ac.jp>, \
            Nozomu Yachie <nzmyachie@gmail.com>"

__date__ = "2020/9/15"

import pandas as pd
import numpy as np

def sub_mat_parser(infile):
    df = pd.read_csv(infile, index_col=0)

    # substitution matrix
    # = [
    # A: [AtoA, AtoT, AtoG, AtoC]
    # T: [TtoA, TtoT, TtoG, TtoC]
    # G: [GtoA, GtoT, GtoG, GtoC]
    # C: [CtoA, CtoT, CtoG, CtoC]
    # ] * L

    substitution_rate_matrix = [[
        np.array([0, 0, 0, 0]),
        np.array([0, 0, 0, 0]),
        np.array([0, 0, 0, 0]),
        np.array([0, 0, 0, 0])
    ] for i in range(len(df.columns))]

    for idx, row in df.iterrows():
        if idx[0] == "A":
            if idx[1] == "T":
                for row_i, item in enumerate(row):
                    substitution_rate_matrix[row_i][0] = substitution_rate_matrix[row_i][0] + np.array([0, item, 0 , 0])
            elif idx[1] == "G":
                for row_i, item in enumerate(row):
                    substitution_rate_matrix[row_i][0] = substitution_rate_matrix[row_i][0] + np.array([0, 0, item , 0])
            elif idx[1] == "C":
                for row_i, item in enumerate(row):
                    substitution_rate_matrix[row_i][0] = substitution_rate_matrix[row_i][0] + np.array([0, 0, 0 , item])

        elif idx[0] == "T":
            if idx[1] == "A":
                for row_i, item in enumerate(row):
                    substitution_rate_matrix[row_i][1] = substitution_rate_matrix[row_i][1] + np.array([item, 0, 0, 0])
            elif idx[1] == "G":
                for row_i, item in enumerate(row):
                    substitution_rate_matrix[row_i][1] = substitution_rate_matrix[row_i][1] + np.array([0, 0, item, 0])
            elif idx[1] == "C":
                for row_i, item in enumerate(row):
                    substitution_rate_matrix[row_i][1] = substitution_rate_matrix[row_i][1] + np.array([0, 0, 0 , item])

        elif idx[0] == "G":
            if idx[1] == "A":
                for row_i, item in enumerate(row):
                    substitution_rate_matrix[row_i][2] = substitution_rate_matrix[row_i][2] + np.array([item, 0, 0, 0])
            elif idx[1] == "T":
                for row_i, item in enumerate(row):
                    substitution_rate_matrix[row_i][2] = substitution_rate_matrix[row_i][2] + np.array([0, item, 0, 0])
            elif idx[1] == "C":
                for row_i, item in enumerate(row):
                    substitution_rate_matrix[row_i][2] = substitution_rate_matrix[row_i][2] + np.array([0, 0, 0 , item])

        elif idx[0] == "C":
            if idx[1] == "A":
                for row_i, item in enumerate(row):
                    substitution_rate_matrix[row_i][3] = substitution_rate_matrix[row_i][3] + np.array([item, 0, 0, 0])
            elif idx[1] == "T":
                for row_i, item in enumerate(row):
                    try:
                        substitution_rate_matrix[row_i][3] = substitution_rate_matrix[row_i][3] + np.array([0, item, 0, 0])
                    except:
                        print(substitution_rate_matrix[row_i][3])
            elif idx[1] == "G":
                for row_i, item in enumerate(row):
                    substitution_rate_matrix[row_i][3] = substitution_rate_matrix[row_i][3] + np.array([0, 0, 0 , item])
        
    for idx, matrix in enumerate(substitution_rate_matrix):
        substitution_rate_matrix[idx][0][0] =  1 - (matrix[0][1] + matrix[0][2] + matrix[0][3])
        substitution_rate_matrix[idx][1][1] =  1 - (matrix[1][0] + matrix[1][2] + matrix[1][3])
        substitution_rate_matrix[idx][2][2] =  1 - (matrix[2][0] + matrix[2][1] + matrix[2][3])
        substitution_rate_matrix[idx][3][3] =  1 - (matrix[3][0] + matrix[3][1] + matrix[3][2])
    return substitution_rate_matrix