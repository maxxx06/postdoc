import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))


import pytest
from unittest import mock
import yaml
import cobra
from riptide_sampling import load_riptide, load_analysis, get_flux_samples 

@pytest.fixture
def mock_config():
    """Fixture pour simuler le fichier de configuration YAML."""
    return {
        "DATA": {
            "path": "mock_data_path",
            "attributes": "mock_attributes",
            "mapped_gene": "mock_mapped_gene",
            "output_transcriptomic_folder": "mock_output_folder/",
            "sacrific_period": "6 hr",
            "dose_level": "High",
            "compound": "amiodarone",
            "replicates": "3"
        },
        "RIPTiDe": {
            "samples": 101,
            "maxfit": True
        },
        "ANALYSIS": {
            "tool": "mana/riptide",
            "pathway_model_description": "mock_model.xml",
            "samples": "sample1,sample2",
            "tag": "tag1/tag2",
            "dose_level": "High/Low",
            "replicates": "3/4",
            "molecules": "amiodarone/valproic_acid",
            "output_dar_files": "mock_output_dar"
        },
        "COMPARISON": {
            "samples_riptide": "mock_samples_riptide",
            "samples_mana": "mock_samples_mana",
            "output": "mock_output"
        }
    }

@pytest.fixture
def mock_dependencies():
    """Fixture pour simuler les dépendances externes."""
    with mock.patch("riptide_sampling.utils.load_model") as load_model, \
         mock.patch("riptide_sampling.utils.create_directory_if_not_exists") as create_directory, \
         mock.patch("riptide_sampling.data_management.data_filter") as data_filter, \
         mock.patch("riptide_sampling.data_management.map_genes") as map_genes, \
         mock.patch("riptide_sampling.riptide.maxfit") as maxfit, \
         mock.patch("riptide_sampling.riptide.contextualize") as contextualize, \
         mock.patch("riptide_sampling.riptide.save_output") as save_output, \
         mock.patch("riptide_sampling.sampling_coverage.r2_iteration") as r2_iteration, \
         mock.patch("riptide_sampling.sampling_coverage.intra_molecules_inter_dar_and_context_methods") as intra_dar:

        load_model.return_value = mock.Mock()  # Simule un modèle COBRA
        create_directory.return_value = True
        data_filter.return_value = "filtered_data"
        map_genes.return_value = {'mapped_data': ['data1', 'data2', 'data3']}
        maxfit.return_value = "riptide_maxfit_object"
        contextualize.return_value = "riptide_contextualize_object"
        save_output.return_value = None
        r2_iteration.return_value = None
        intra_dar.return_value = None

        yield {
            'load_model': load_model,
            'create_directory': create_directory,
            'data_filter': data_filter,
            'map_genes': map_genes,
            'maxfit': maxfit,
            'contextualize': contextualize,
            'save_output': save_output,
            'r2_iteration': r2_iteration,
            'intra_dar': intra_dar
        }

def test_load_riptide(mock_config, mock_dependencies):
    """Test la fonction load_riptide en simulant un fichier de configuration."""
    with mock.patch("builtins.open", mock.mock_open(read_data=yaml.dump(mock_config))):
        load_riptide("mock_config.yaml")

    mock_dependencies['create_directory'].assert_called()
    mock_dependencies['data_filter'].assert_called()
    mock_dependencies['map_genes'].assert_called()
    mock_dependencies['maxfit'].assert_called()

def test_load_analysis(mock_config, mock_dependencies):
    """Test la fonction load_analysis en simulant le chargement du YAML et des appels à sampling_coverage."""
    with mock.patch("builtins.open", mock.mock_open(read_data=yaml.dump(mock_config))), \
         mock.patch("riptide_sampling.cobra.io.read_sbml_model", return_value=mock.Mock()) as read_model:

        load_analysis("mock_config.yaml")

    read_model.assert_called_once_with("mock_model.xml")
    mock_dependencies['r2_iteration'].assert_called()
    mock_dependencies['intra_dar'].assert_called()

def test_get_flux_samples(mock_dependencies):
    """Test la fonction get_flux_samples avec des valeurs fictives."""
    get_flux_samples(
        data_path="mock_data_path",
        attribute_data_path="mock_attributes",
        sacrific_period="6 hr",
        dose_level="High",
        compound_name="amiodarone",
        replicates=3,
        samples=101,
        maxfit=True,
        output_transcriptomic="mock_output_folder/",
        mapped_gene="mock_mapped_gene"
    )

    mock_dependencies['create_directory'].assert_called()
    mock_dependencies['data_filter'].assert_called()
    mock_dependencies['map_genes'].assert_called()
    mock_dependencies['maxfit'].assert_called()
    mock_dependencies['save_output'].assert_called()
