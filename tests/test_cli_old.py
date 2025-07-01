"""
Unit tests for CLI module.
"""

import unittest
import tempfile
import json
import logging
from pathlib import Path
from unittest.mock import patch, Mock, MagicMock
from click.testing import CliRunner
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

from nova.cli import cli, setup_logging


class TestCLI(unittest.TestCase):
    """Test cases for CLI functionality."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.runner = CliRunner()
    
    def test_setup_logging_info_level(self):
        """Test logging setup at INFO level."""
        with patch('logging.basicConfig') as mock_config:
            setup_logging(verbose=False)
            mock_config.assert_called_once()
            args, kwargs = mock_config.call_args
            self.assertEqual(kwargs['level'], logging.INFO)
    
    def test_setup_logging_debug_level(self):
        """Test logging setup at DEBUG level."""
        with patch('logging.basicConfig') as mock_config:
            setup_logging(verbose=True)
            mock_config.assert_called_once()
            args, kwargs = mock_config.call_args
            self.assertEqual(kwargs['level'], logging.DEBUG)
    
    def test_cli_help(self):
        """Test CLI help output."""
        result = self.runner.invoke(cli, ['--help'])
        self.assertEqual(result.exit_code, 0)
        self.assertIn('nova: de novo variant insertion simulator', result.output)
    
    def test_validate_config_valid(self):
        """Test validate-config command with valid configuration."""
        config_file = 'tests/test_data/test_config_small.json'
        result = self.runner.invoke(cli, ['validate-config', config_file])
        self.assertEqual(result.exit_code, 0)
    
    def test_validate_config_invalid(self):
        """Test validate-config command with invalid configuration (wrong data type)."""
        config_file = 'tests/test_data/test_config_broken.json'
        result = self.runner.invoke(cli, ['validate-config', config_file])
        self.assertEqual(result.exit_code, 1)
    
    def test_validate_config_malformed_json(self):
        """Test validate-config command with malformed JSON."""
        config_file = 'tests/test_data/test_config_malformed.json'
        result = self.runner.invoke(cli, ['validate-config', config_file])
        self.assertEqual(result.exit_code, 1)
    
    @patch('nova.cli.ReadSelector')
    @patch('nova.cli.ReadInserter')
    @patch('nova.cli.VariantGenerator')
    @patch('nova.cli.VariantRegistry')
    def test_simulate_command_success(self, mock_registry_class, mock_generator_class, 
                                    mock_inserter_class, mock_selector_class):
        """Test successful simulate command execution."""
        # Setup mocks
        mock_registry = Mock()
        mock_registry_class.return_value = mock_registry
        
        mock_generator = Mock()
        mock_generator_class.return_value = mock_generator
        mock_generator.validate_config.return_value = []  # No validation errors
        mock_generator.generate_from_config.return_value = ['id1', 'id2']
        
        mock_selector = Mock()
        mock_selector_class.return_value = mock_selector
        mock_selector.select_reads.return_value = [
            Mock(read_name='read1'),  # LazyReadReference objects
            Mock(read_name='read2')
        ]
        
        mock_inserter = Mock()
        mock_inserter_class.return_value = mock_inserter
        mock_inserter.insert_streaming.return_value = {'total_insertions': 2}
        mock_inserter.get_insertion_statistics.return_value = {'total_insertions': 2}
        mock_registry.get_statistics.return_value = {'total_sequences': 2}
        
        # Create test files
        config = {'random': {'n': 2, 'length': 100}}
        
        with tempfile.TemporaryDirectory() as temp_dir:
            # Create config file
            config_file = Path(temp_dir) / 'config.json'
            with open(config_file, 'w') as f:
                json.dump(config, f)
            
            # Create dummy BAM file (just needs to exist for CLI validation)
            bam_file = Path(temp_dir) / 'test.bam'
            bam_file.touch()
            
            # Create output directory
            output_dir = Path(temp_dir) / 'output'
            
            # Run simulate command
            result = self.runner.invoke(cli, [
                'simulate', 
                str(bam_file), 
                str(config_file),
                '--output-dir', str(output_dir),
                '--output-prefix', 'test_sim',
                '--random-seed', '42'
            ])
            
            # Check success
            self.assertEqual(result.exit_code, 0)
            # Click captures logging output, not stdout, so we check exit code
    
    @patch('nova.cli.VariantGenerator')
    @patch('nova.cli.VariantRegistry')
    def test_simulate_command_validation_failure(self, mock_registry_class, mock_generator_class):
        """Test simulate command with validation failure."""
        # Setup mocks
        mock_registry = Mock()
        mock_registry_class.return_value = mock_registry
        
        mock_generator = Mock()
        mock_generator_class.return_value = mock_generator
        mock_generator.validate_config.return_value = ['Invalid config error']
        
        # Create test files
        config = {'random': {'n': -1, 'length': 0}}
        
        with tempfile.TemporaryDirectory() as temp_dir:
            config_file = Path(temp_dir) / 'config.json'
            with open(config_file, 'w') as f:
                json.dump(config, f)
            
            bam_file = Path(temp_dir) / 'test.bam'
            bam_file.touch()
            
            result = self.runner.invoke(cli, [
                'simulate', 
                str(bam_file), 
                str(config_file)
            ])
            
            self.assertEqual(result.exit_code, 1)
            # Check for logging output in the result
    
    @patch('nova.cli.VariantGenerator')
    @patch('nova.cli.VariantRegistry')
    def test_simulate_command_no_sequences_generated(self, mock_registry_class, mock_generator_class):
        """Test simulate command when no sequences are generated."""
        # Setup mocks
        mock_registry = Mock()
        mock_registry_class.return_value = mock_registry
        
        mock_generator = Mock()
        mock_generator_class.return_value = mock_generator
        mock_generator.validate_config.return_value = []
        mock_generator.generate_from_config.return_value = []  # No sequences generated
        
        config = {'random': {'n': 0, 'length': 100}}
        
        with tempfile.TemporaryDirectory() as temp_dir:
            config_file = Path(temp_dir) / 'config.json'
            with open(config_file, 'w') as f:
                json.dump(config, f)
            
            bam_file = Path(temp_dir) / 'test.bam'
            bam_file.touch()
            
            result = self.runner.invoke(cli, [
                'simulate', 
                str(bam_file), 
                str(config_file)
            ])
            
            self.assertEqual(result.exit_code, 1)
            # Check for logging output in the result
    
    def test_simulate_command_file_not_found(self):
        """Test simulate command with non-existent files."""
        result = self.runner.invoke(cli, [
            'simulate', 
            'nonexistent.bam', 
            'nonexistent.json'
        ])
        
        self.assertEqual(result.exit_code, 2)  # Click file not found error
    
    def test_cli_verbose_flag(self):
        """Test CLI with verbose flag."""
        result = self.runner.invoke(cli, ['--verbose', '--help'])
        self.assertEqual(result.exit_code, 0)
        # The verbose flag affects logging level, which we can't easily test in Click
    
    def test_cli_without_verbose_flag(self):
        """Test CLI without verbose flag."""
        result = self.runner.invoke(cli, ['--help'])
        self.assertEqual(result.exit_code, 0)
        # The verbose flag affects logging level, which we can't easily test in Click


if __name__ == '__main__':
    unittest.main()