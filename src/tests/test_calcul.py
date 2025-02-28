import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import pytest
import pandas as pd
import numpy as np
from scipy.stats import chi2_contingency, ks_2samp

from calcul import (
    compute_freq_table_df, compute_r2_df, chi2_independance_test, 
    compute_dar_specificity_ratio_inter_molecules, compute_dar_specificity_ratio_intra_molecules,
    compute_two_sample_KS_from_df
)

@pytest.fixture
def sample_dataframe():
    return pd.DataFrame({
        'A': [1, 2, 3, 4, 5],
        'B': [2, 3, 4, 5, 6],
        'C': [0, 1, 0, 1, 0]
    })

@pytest.fixture
def sample_dataframe_bin():
    return pd.DataFrame({
        'A': [1, 1, 1, 1, 1],
        'B': [1, 1, 1, 1, 1],
        'C': [0, 1, 0, 1, 0]
    })


@pytest.fixture
def sample_r2_dataframe():
    return pd.DataFrame({'0': [0.1, 0.2, 0.3]}, index=['r1', 'r2', 'r3'])

@pytest.fixture
def sample_df_trmt():
    return pd.DataFrame({'0': [0.15, 0.25, 0.35]}, index=['r1', 'r2', 'r3'])

@pytest.fixture
def empty_r2_dataframe():
    return pd.DataFrame(columns=['val'])

# Test compute_freq_table_df
def test_compute_freq_table_df(sample_dataframe,sample_dataframe_bin):
    freq_table_bin = compute_freq_table_df(sample_dataframe_bin, 'test','mana')
    expected_bin = pd.Series([1.0, 1.0, 0.4], name='test')
    pd.testing.assert_series_equal(freq_table_bin, expected_bin)

    freq_table = compute_freq_table_df(sample_dataframe, 'test','riptide')
    expected = pd.Series([1.0, 1.0, 0.4], name='test')
    pd.testing.assert_series_equal(freq_table, expected)

    pd.testing.assert_series_equal(freq_table_bin, freq_table)

# Test compute_r2_df
def test_compute_r2_df(sample_r2_dataframe, sample_df_trmt, empty_r2_dataframe):
    r2_result = compute_r2_df(empty_r2_dataframe, sample_r2_dataframe, sample_df_trmt)
    assert isinstance(r2_result, pd.DataFrame)
    assert 'val' in r2_result.columns

# Test chi2_independance_test
def test_chi2_independance_test(tmp_path, sample_dataframe):
    stat_dar_file = tmp_path / "test_output.tsv"
    chi2_independance_test(sample_dataframe, sample_dataframe, str(stat_dar_file))
    assert stat_dar_file.exists()

# Test compute_dar_specificity_ratio_inter_molecules (mocked paths & data)
def test_compute_dar_specificity_ratio_inter_molecules(mocker):
    mocker.patch("os.path.exists", return_value=True)
    mocker.patch("pandas.read_csv", return_value=pd.DataFrame(index=['rxn1', 'rxn2']))
    mocker.patch("sampling_coverage.generate_annotation_table", return_value=pd.DataFrame({'Pathway in model': ['p1', 'p2']}))
    mocker.patch("matplotlib.pyplot.savefig")
    df = pd.DataFrame()
    with pytest.raises(Exception):
        result_df = compute_dar_specificity_ratio_inter_molecules("path", "mana", "reps", "dose", df, ["mol1", "mol2"], None, "p1", "p2", "p3", "p4", "p5")
    result_df = compute_dar_specificity_ratio_inter_molecules("path", "mana", "reps", "dose", df, ["mol1", "mol2"], None, "p1", "p2", "p3", "p4", "p5",tag="r2")
    
    assert isinstance(result_df, pd.DataFrame)

# Test compute_dar_specificity_ratio_intra_molecules (mocked paths & data)
def test_compute_dar_specificity_ratio_intra_molecules(mocker):
    mocker.patch("os.path.exists", return_value=True)
    mocker.patch("pandas.read_csv", return_value=pd.DataFrame(index=['rxn1', 'rxn2']))
    mocker.patch("matplotlib.pyplot.savefig")
    df = pd.DataFrame()
    result_df = compute_dar_specificity_ratio_intra_molecules("path", "mana", "mol1", "reps", "dose", df, "images", "tag", "files")
    assert isinstance(result_df, pd.DataFrame)

# Test compute_two_sample_KS_from_df
def test_compute_two_sample_KS_from_df(tmp_path, sample_dataframe):
    stat_dar_file = tmp_path / "ks_output.tsv"
    compute_two_sample_KS_from_df(sample_dataframe, sample_dataframe, str(stat_dar_file))
    assert stat_dar_file.exists()
