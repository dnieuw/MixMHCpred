#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import numpy as np
import pandas as pd
import os
import importlib.resources as pkg_resources
from typing import List, Optional, Tuple, Dict, Any, Sequence

AA = ['A','C','D','E','F','G','H','I','K','L','M','N','P','Q','R','S','T','V','W','Y']

def creation_PWM_dict_spec(PWMs: List[List[pd.DataFrame]], AAs: List[str] = AA, PL: int = 9) -> List[List[Dict[str, float]]]:
    """Convert a list of PWM dataframes into a list of dictionaries keyed by AA+position.

    Args:
        PWMs: nested list (per-allele) of pandas DataFrames containing PWM values.
        AAs: list of amino-acid single-letter codes.
        PL: peptide length (number of positions).

    Returns:
        Nested list where each inner list contains dicts mapping keys like 'A0' -> score.
    """
    pos = [f'{i}' for i in range(PL)]
    PWMs_dict = []
    
    for PWM  in PWMs:
        PWM_D = []
        for PW in PWM:
            PWM_d = dict({})
            for i in range (len(AAs)):
                for j in range (len(pos)):
                    PWM_d[AAs[i] + pos[j]] = PW.iloc[i,j]
            PWM_D.append(PWM_d)
        PWMs_dict.append(PWM_D)

    return PWMs_dict


def interpolate_v2(x_values: np.ndarray, y_values: np.ndarray, x_new_values: List) -> np.ndarray:
    """Linear interpolation over reversed x/y arrays and clip out-of-range values.

    The function expects x_values and y_values to be sorted; it reverses them and
    performs interpolation for x_new_values.
    """
    x_values_reversed = x_values[::-1]
    y_values_reversed = y_values[::-1]

    idx = np.searchsorted(x_values_reversed, x_new_values, side='left') - 1
    idx = np.clip(idx, 0, len(x_values_reversed) - 2)

    x0 = x_values_reversed[idx]
    x1 = x_values_reversed[idx + 1]
    y0 = y_values_reversed[idx]
    y1 = y_values_reversed[idx + 1]

    y_new = y0 + (x_new_values - x0) * (y1 - y0) / (x1 - x0)

    # Set out-of-range values to the minimum or maximum of y_values
    y_new[x_new_values < x_values_reversed[0]] = y_values_reversed[0]
    y_new[x_new_values > x_values_reversed[-1]] = y_values_reversed[-1]

    return y_new


def ligands_PWMs_scores_spec(
    PWMs_dict: List[List[Dict[str, float]]],
    peptides: List[str],
    alphas: List[List[float]],
    bias: List[List[float]],
    std_dev: List[List[float]],
    k: int,
    perRank: List[Tuple[np.ndarray, np.ndarray]],
) -> Tuple[List[List[float]], List[List[float]]]:
    """Compute corrected scores and interpolated ranks for a set of peptides.

    Args:
        PWMs_dict: nested dictionary representation of PWMs for this peptide length.
        peptides: list of peptide strings.
        alphas: list of alpha weights per allele/component.
        bias: bias shifts used for score correction.
        std_dev: standard deviation used for score normalization.
        k: index offset for length-dependent arrays.
        perRank: list of tuples (scores_array, ranks_array) for interpolation.

    Returns:
        A tuple (allele_scores, allele_ranks) where each is a list (per-allele) of lists of floats.
    """
    allele_scores = []
    allele_ranks = []
    LL = len(peptides[0]) if peptides else 0

    for i in range(len(PWMs_dict)):
        scores_corr=[]

        for x in peptides:
            score = 0
            for z in range(len(alphas[i])):
                s = 1.0
                for j, a in enumerate(x):
                    s *= PWMs_dict[i][z][f'{a}{j}']
                score = score + (alphas[i][z] * s)

            score = np.log(score) / LL
            aa = (score - bias[i][k]) / (std_dev[i][k])
            scores_corr.append(float(aa))

        allele_scores.append(scores_corr)

        allele_ranks.append(interpolate_v2(perRank[i][0], perRank[i][1], scores_corr))

    return allele_scores, allele_ranks
    

def Ligands_scores_spec(
    ligands: pd.DataFrame,
    ligands_L: np.ndarray,
    Alleles: List[str],
    PWMs_dict: List[Any],
    alphas: List[Any],
    bias: List[Any],
    std_dev: List[Any],
    perRank: List[Any],
) -> pd.DataFrame:
    """Compute scores and ranks for all ligands across requested alleles.

    This wraps the detailed per-length scoring function and aggregates results into a DataFrame.
    """
    ligands = ligands.sort_values('length', ascending=True).reset_index()

    SCORES = []
    RANKS = []
    for i in range(len(ligands_L)):
        peptides = list(ligands[ligands['length'] == ligands_L[i]]['Peptide'])
        L_i = ligands_L[i] - 8
        allele_scores, allele_ranks = ligands_PWMs_scores_spec(PWMs_dict[L_i], peptides, alphas[L_i], bias, std_dev, L_i, perRank)

        SCORES.append(allele_scores)
        RANKS.append(allele_ranks)

    for x in range(len(Alleles)):
        
        scores = []
        ranks = []
        for i in range(len(ligands_L)):
            scores = scores + list(SCORES[i][x])
            ranks = ranks + list(RANKS[i][x])

        temp_df = pd.DataFrame({f'Score_{Alleles[x]}': scores, f'%Rank_{Alleles[x]}': ranks})
        ligands = pd.concat([ligands, temp_df], axis=1)

    ligands.insert(loc=3, column='BestAllele', value=ligands[[f'%Rank_{a}' for a in Alleles]].idxmin(axis=1))
    ligands['BestAllele'] = [a[6:] for a in ligands['BestAllele']]
    ligands.insert(loc=2, column='Score_bestAllele', value=ligands[[f'Score_{a}' for a in Alleles]].max(axis=1))
    ligands.insert(loc=4, column='%Rank_bestAllele', value=ligands[[f'%Rank_{a}' for a in Alleles]].min(axis=1))

    ligands = ligands.sort_values('index', ascending=True).reset_index(drop=True)
    ligands = ligands.drop(columns=['index', 'length'])

    return ligands


def get_allele_index() -> List[int]:
    """Read binding site indices from a file and convert to zero-based positions."""
    with pkg_resources.files('MixMHCpred.data.Allele_pos').joinpath('binding_sites.txt').open('r') as f:
        return [int(int(val) - 2) for val in f.read().strip().split(' ')]


def pairwise_score(seq1: str, seq2: str, matrix: Dict[str, float]) -> float:
    """Score two sequences pairwise using a pair-score matrix (dict keyed by pairs)."""
    score = 0.0
    for i in range(len(seq2)):
        pair = seq1[i] + seq2[i]
        if pair in matrix:
            score += matrix[pair]
        else:
            raise ValueError(str(pair) + ' is not in the matrix')
    return score


def best_score(x: str, K_seq: List[str], matrix: Dict[str, float]) -> Tuple[float, int]:
    """Return the best pairwise score and index among a list of sequences K_seq against x."""
    score, index = 0.0, 0

    for pos, y in enumerate(K_seq):
        scores = pairwise_score(x, y, matrix)
        if scores > score:
            score, index = float(scores), pos
    return float(score), int(index)


def sim_sequence(Alleles: List[str], K_seq: List[str], mutant_seq: str, matrix: Dict[str, float]) -> Tuple[str, float]:
    """Find the most similar allele sequence from K_seq to mutant_seq using the provided matrix.

    Returns the allele identifier and a similarity score in [0,1].
    """
    score, index = best_score(mutant_seq, K_seq, matrix)
    Alleles_sim = Alleles[index]
    denom = np.sqrt(pairwise_score(mutant_seq, mutant_seq, matrix) * pairwise_score(K_seq[index], K_seq[index], matrix))
    score = 1 - (score / denom)

    return Alleles_sim, float(score)


def euclidean_distance(pwm_a: np.ndarray, pwm_b: np.ndarray) -> np.ndarray:
    """Compute Euclidean distance between two PWM arrays along axis 0."""
    return np.sqrt(np.sum((pwm_a - pwm_b) ** 2, axis=0))


def distance_to_training(training_Alleles: List[str], predicted_Alleles: List[str]) -> Tuple[List[str], List[float]]:
    """Compute similarity of predicted alleles to training alleles using a BLOSUM-derived matrix.

    Returns a tuple (closest_alleles, similarity_scores).
    """
    Allele_Pos_idx = get_allele_index()
    with pkg_resources.files('MixMHCpred.data').joinpath('MHC_I_sequences.txt').open('r') as f:
        seq_data = pd.read_csv(f, sep=r'\s+')
    for x in training_Alleles:
        if x not in seq_data['Allele'].tolist():
            print(f'{x} not found in sequence database')
    Alleles_seq_training = [seq_data[seq_data['Allele'] == allele]['Sequence'].iloc[0] for allele in training_Alleles]
    Alleles_seq_predicted = [seq_data[seq_data['Allele'] == allele]['Sequence'].iloc[0] for allele in predicted_Alleles]

    with pkg_resources.files('MixMHCpred.data').joinpath('blosum62_update.npy').open('rb') as f:
        dict_load = np.load(f, allow_pickle=True)
    blosum62 = dict_load.item()

    training_sequences = ["".join([allele[z] for z in Allele_Pos_idx]) for allele in Alleles_seq_training]
    predicted_sequences = ["".join([allele[z] for z in Allele_Pos_idx]) for allele in Alleles_seq_predicted]

    close_Alleles = []
    sim_scores = []


    for i in range(len(predicted_sequences)):
        Alleles_sim, score = sim_sequence(training_Alleles, training_sequences, predicted_sequences[i], blosum62)
        close_Alleles.append(Alleles_sim)
        sim_scores.append(float(score))

    return close_Alleles, sim_scores


def normalize_alleles(alleles_csv: str) -> List[str]:
    """Normalize allele strings from user input to the project's format.

    Examples of normalization handled:
    - HLA-A*01:01  -> A0101
    - H-2-Db -> H2Db
    """
    alleles_input = [s.strip() for s in alleles_csv.split(',') if s.strip()]
    alleles_out: List[str] = []
    for h in alleles_input:
        if h.startswith('HLA-'):
            h = h[4:]
        if len(h) > 1 and h[1] == '*':
            h = h[0] + h[2:]
        if len(h) > 3 and h[3] == ':':
            h = h[:3] + h[4:]
        if h.startswith('H-2'):
            h = 'H2' + h[3:]
        h = h.replace(':', '').replace('*', '')
        alleles_out.append(h)
    return alleles_out


def read_peptides(file_input: str) -> List[str]:
    """Read peptide sequences from a simple one-per-line file (FASTA headers ignored)."""
    if not os.path.isfile(file_input):
        raise FileNotFoundError(f"No such file: {file_input}")
    with open(file_input, 'r') as f:
        peptides = [line.strip() for line in f if not line.startswith('>') and line.strip()]
    if not peptides:
        raise ValueError("Empty peptide file")
    return peptides


def validate_peptides(peptides: List[str], Lmin: int = 8, Lmax: int = 14) -> None:
    """Validate peptide characters and lengths. Raises ValueError on failure."""
    amino_acids_map = {
        "A": 0, "C": 1, "D": 2, "E": 3, "F": 4, "G": 5, "H": 6, "I": 7,
        "K": 8, "L": 9, "M": 10, "N": 11, "P": 12, "Q": 13, "R": 14, "S": 15,
        "T": 16, "V": 17, "W": 18, "Y": 19
    }

    for p in peptides:
        try:
            _ = [amino_acids_map[char] for char in p]
        except KeyError as e:
            raise ValueError(f"Error in input sequence - unknown amino acids {e.args[0]}: {p}")
        if not (Lmin <= len(p) <= Lmax):
            raise ValueError(f"Incompatible peptide length: {p}\t{len(p)}. Only peptides of length {Lmin}-{Lmax} are supported")


def load_pwm_data(alleles_in: List[str], L: List[int]):
    """Load PWMs, alphas, perRank, bias and standard_dev for alleles present in the datarary.

    Returns a tuple: (PWMs_pred_dict_list, alphas_list, Alleles_perRank_f, bias, standard_dev)
    """
    PWMs_pred_dict = []
    alphas = []
    Alleles_perRank_f = []
    bias = []
    standard_dev = []

    with pkg_resources.files('MixMHCpred.data').joinpath('alleles_list.txt').open('r') as f:
        alleles_list = pd.read_csv(f, sep='\t')

    if len(alleles_in) > 0:
        PWMs_pred = []
        for l in L:
            with pkg_resources.files('MixMHCpred.data.pwm').joinpath(f'class1_{l}/alphas.txt').open('r') as f:
                ratios = pd.read_csv(f, sep='\t')
            PWMs_pred.append([])
            alphas.append([])

            for xxx in alleles_in:
                PWM = []
                alpha = []
                n = alleles_list[alleles_list['Allele'] == xxx][f'{l}'].iloc[0]
                for i in range(n):
                    with pkg_resources.files('MixMHCpred.data.pwm').joinpath(f'class1_{l}/PWM_{xxx}_{i+1}.csv').open('r') as f:
                        PWM.append(pd.read_csv(f, index_col=0))
                    alpha.append(ratios[ratios['Allele'] == f'{xxx}_{i+1}']['ratio'].iloc[0])
                PWMs_pred[-1].append(PWM)
                alphas[-1].append(alpha)

            PWMs_pred_dict.append(creation_PWM_dict_spec(PWMs_pred[-1], AAs=AA, PL=l))

        Alleles_perRank = []
        for zz in alleles_in:
            with pkg_resources.files('MixMHCpred.data.PerRank').joinpath(f'{zz}.txt').open('r') as f:
                Alleles_perRank.append(pd.read_csv(f, sep='\t'))
    Alleles_perRank_f = [[Alleles_perRank[i]['score'].to_numpy(), Alleles_perRank[i]['rank'].to_numpy()] for i in range(len(alleles_in))]
    with pkg_resources.files('MixMHCpred.data.shifts').joinpath('bias.txt').open('r') as f:
        bias = pd.read_csv(f, sep='\t', index_col=0).loc[alleles_in].to_numpy().tolist()
    with pkg_resources.files('MixMHCpred.data.shifts').joinpath('standard_dev.txt').open('r') as f:
        standard_dev = pd.read_csv(f, sep='\t', index_col=0).loc[alleles_in].to_numpy().tolist()

    return PWMs_pred_dict, alphas, Alleles_perRank_f, bias, standard_dev


def create_header(alleles_to_test: List[str], closest_alleles: List[str], distance_scores: List[float], file_input: str) -> List[str]:
    Alleles_quality_info = [f'{closest_alleles[i]} ({np.round(distance_scores[i],4)})' for i in range(len(closest_alleles))]
    header_comments = [
        "####################",
        "# Output from MixMHCpred (v3.0)",
        f"# Alleles: {', '.join(alleles_to_test)}",
        f"# Closest Allele (distance score): {' -- '.join(Alleles_quality_info)}",
        f"# Input file: {file_input}",
        "# MixMHCpred is freely available for academic users.",
        "# Private companies should contact Nadette Bulgin (nbulgin@lcr.org) at the Ludwig Institute for Cancer Research Ltd for commercial licenses.",
        "# To cite MixMHCpred3.0, please refer to:",
        "# Tadros et al., Predicting MHC-I ligands across alleles and species: How far can we go?, BioRxiv (2024).",
        "####################"
    ]

    return header_comments


def run_MixMHCpred(file_input: str, alleles_raw: str) -> Tuple[List[str], pd.DataFrame]:
    """Main entrypoint for the refactored script."""

    # Normalize alleles
    alleles_to_test = normalize_alleles(alleles_raw)

    # Read peptides and validate
    peptides = read_peptides(file_input)
    validate_peptides(peptides)

    # Prepare ligand DataFrame
    Ligands = pd.DataFrame({'Peptide': peptides})
    Ligands['length'] = Ligands.Peptide.str.len()
    Ligands_L = np.unique(Ligands['length'])

    # Load alleles present in library
    with pkg_resources.files('MixMHCpred.data').joinpath('alleles_list.txt').open('r') as f:
        alleles_list = pd.read_csv(f, sep='\t')
    Alleles = alleles_list['Allele'].tolist()
    Alleles_in = np.sort(list(set(alleles_to_test) & set(Alleles))).tolist()
    Alleles_out = np.sort(list(set(alleles_to_test) - set(Alleles))).tolist()

    L = [8,9,10,11,12,13,14]

    PWMs_pred_dict, alphas, Alleles_perRank_f, bias, standard_dev = load_pwm_data(Alleles_in, L)

    Alleles_to_testt = Alleles_in + Alleles_out

    ligands = Ligands_scores_spec(Ligands, Ligands_L, Alleles_to_testt, PWMs_pred_dict, alphas, bias, standard_dev, Alleles_perRank_f)

    columns_df = ['Peptide','Score_bestAllele','BestAllele','%Rank_bestAllele']
    for x in alleles_to_test:
        columns_df += [f'Score_{x}', f'%Rank_{x}']
    ligands = ligands[columns_df]
    ligands = ligands.round(6)

    closest_alleles, distance_scores = distance_to_training(Alleles, alleles_to_test)

    header_comments = create_header(alleles_to_test, closest_alleles, distance_scores, file_input)
    
    return header_comments, ligands
