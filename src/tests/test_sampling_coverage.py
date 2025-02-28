import os, sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))

import pytest
from unittest import mock
import time
import pandas as pd
import matplotlib.pyplot as plt
import utils
from unittest.mock import patch, MagicMock

from sampling_coverage import (
    r2_iteration,
    run_tools,
    replicat_intra_molecule,
    inter_molecules_and_dar
)


def test_r2_iteration():
    # Define test inputs
    path_samples = "test_samples"
    tool = "test_tool"
    tag = "test_tag"
    model = "test_model"
    dose_level = "high"
    replicates = 3
    molecules = "test_molecule"
    path_dar_files = "test_dar_files"
    
    # Patch external function calls
    with patch("sampling_coverage.run_tools") as mock_run_tools, \
         patch("sampling_coverage.replicat_intra_molecule") as mock_replicat_intra, \
         patch("sampling_coverage.inter_molecules_and_dar") as mock_inter_molecules, \
         patch("sampling_coverage.intra_molecules_and_inter_dar") as mock_intra_molecules, \
         patch("builtins.print") as mock_print:
        
        # Call the function
        r2_iteration(path_samples, tool, tag, model, dose_level, replicates, molecules, path_dar_files)
        
        # Assert that the mocked functions were called with expected arguments
        mock_run_tools.assert_called_once_with(tool, dose_level, replicates, molecules, path_samples, path_dar_files, tag)
        mock_replicat_intra.assert_called_once_with(path_samples, model, dose_level, replicates, molecules, path_dar_files, tag, tool)
        mock_inter_molecules.assert_called_once_with(path_samples, model, dose_level, replicates, molecules, tool, tag, path_dar_files)
        mock_intra_molecules.assert_called_once_with(path_samples, model, dose_level, replicates, molecules, tool, tag, path_dar_files)
        
        # Check if print statements were executed
        assert mock_print.call_count > 0

def test_run_tools_with_mocked_files():
    """
    Test run_tools function with mocked file system interactions.
    """
    tool = "riptide"
    dose_level = ["low", "medium", "high"]
    replicates = ["rep1", "rep2"]
    molecules = ["mol1"]
    path_samples_init = "/fake/path_samples"
    path_dar_files = "/dar_results"
    tag = "r2"
    
    # Mocking file existence check
    with patch("os.path.exists", return_value=True), \
         patch("pandas.read_csv", return_value=pd.DataFrame({"col1": [1, 2], "col2": [3, 4]})), \
         patch("pandas.DataFrame.to_csv") as mock_to_csv, \
         patch("utils.create_directory_if_not_exists") as mock_create_dir, \
         patch("utils.read_sample_file", return_value=(None, pd.DataFrame({"col1": [1,0], "col2": [2,0]}))), \
         patch("calcul.compute_freq_table_df", return_value=pd.DataFrame({"freq_col": [1,0]})), \
         patch("calcul.compute_r2_df", return_value=pd.DataFrame({"r2": [0.9]})), \
         patch("calcul.chi2_independance_test") as mock_chi2_test, \
         patch("calcul.compute_two_sample_KS_from_df") as mock_ks_test:
        
        run_tools(tool, dose_level, replicates, molecules, path_samples_init, path_dar_files, tag)
        
        # Check if directories are created
        mock_create_dir.assert_called()
        
        # Check if files are being saved
        assert mock_to_csv.call_count > 0  # At least one file should be written
        
        # Check if statistical tests are called
        mock_ks_test.assert_called()

    tool = "mana"
    dose_level = ["low", "medium", "high"]
    replicates = ["rep1", "rep2"]
    molecules = ["mol1"]
    path_samples_init = "/fake/path_samples"
    path_dar_files = "/dar_results"
    tag = "r2"
    
    # Mocking file existence check
    with patch("os.path.exists", return_value=True), \
         patch("pandas.read_csv", return_value=pd.DataFrame({"col1": [1, 1], "col2": [1, 0]})), \
         patch("pandas.DataFrame.to_csv") as mock_to_csv, \
         patch("utils.create_directory_if_not_exists") as mock_create_dir, \
         patch("utils.read_sample_file", return_value=(None, pd.DataFrame({"col1": [1,0], "col2": [1,0]}))), \
         patch("calcul.compute_freq_table_df", return_value=pd.DataFrame({"freq_col": [1,0]})), \
         patch("calcul.compute_r2_df", return_value=pd.DataFrame({"r2": [0.9]})), \
         patch("calcul.chi2_independance_test") as mock_chi2_test, \
         patch("calcul.compute_two_sample_KS_from_df") as mock_ks_test:
        
        run_tools(tool, dose_level, replicates, molecules, path_samples_init, path_dar_files, tag)
        
        # Check if directories are created
        mock_create_dir.assert_called()
        
        # Check if files are being saved
        assert mock_to_csv.call_count > 0  # At least one file should be written
        
        # Check if statistical tests are called
        mock_chi2_test.assert_called() 

        print("Test passed: run_tools correctly interacts with file system and computations.")

@pytest.fixture
def mock_dependencies():
    with patch("utils.create_directory_if_not_exists") as mock_mkdir, \
         patch("calcul.compute_dar_specificity_ratio_inter_molecules", return_value=pd.DataFrame({"A": [1, 2, 3]})) as mock_compute_dar, \
         patch("utils.remove_nan_from_list", side_effect=lambda x: x) as mock_remove_nan, \
         patch("sampling_coverage.generate_annotation_table", return_value=pd.DataFrame({"Pathway in model": ["path1", "path2"]})) as mock_annot, \
         patch("matplotlib.pyplot.savefig") as mock_savefig:
        yield {
            "mock_mkdir": mock_mkdir,
            "mock_compute_dar": mock_compute_dar,
            "mock_remove_nan": mock_remove_nan,
            "mock_annot": mock_annot,
            "mock_savefig": mock_savefig,
        }

def test_inter_molecules_and_dar(mock_dependencies):
    path_dar = "test/path"
    model = MagicMock()
    dose_level = ["Low", "High"]
    replicates = ["rep1", "rep2"]
    molecules = ["mol1", "mol2"]
    tool = "mana"
    tag = "chi2"
    path_dar_files = "dar_files/"
    
    with patch("pandas.DataFrame.to_csv") as mock_to_csv:
        inter_molecules_and_dar(path_dar, model, dose_level, replicates, molecules, tool, tag, path_dar_files)
    
    # Check if directories were created
    mock_dependencies["mock_mkdir"].assert_any_call("test/path/inter_molecules/files/")
    mock_dependencies["mock_mkdir"].assert_any_call("test/path/inter_molecules/images/")
    
    # Ensure compute_dar_specificity_ratio_inter_molecules was called for each replicate
    assert mock_dependencies["mock_compute_dar"].call_count == len(replicates) * 2
    
    # Verify DataFrame is saved correctly
    assert mock_to_csv.call_count >= 3  # At least df_r2, df_ks, and df_comparison should be saved
    
    # Ensure annotation generation was called for Venn diagram
    assert mock_dependencies["mock_annot"].call_count >= 2
    
    # Check if Venn diagrams were saved
    assert mock_dependencies["mock_savefig"].call_count >= 1

