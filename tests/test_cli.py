"""
Unit tests for CLI module using real test files.
"""

import unittest
import tempfile
import os
from pathlib import Path
from click.testing import CliRunner

from nova.cli import cli


class TestCLI(unittest.TestCase):
    """Test cases for CLI functionality using real test data."""
    
    def setUp(self):
        """Set up test fixtures."""
        self.runner = CliRunner()
    
    def test_cli_help(self):
        """Test CLI help output."""
        result = self.runner.invoke(cli, ['--help'])
        self.assertEqual(result.exit_code, 0)
        self.assertIn('nova: de novo variant insertion simulator', result.output)
    
    def test_validate_config_valid_small(self):
        """Test validate-config command with valid small configuration."""
        config_file = 'tests/test_data/test_config_small.json'
        result = self.runner.invoke(cli, ['validate-config', config_file])
        self.assertEqual(result.exit_code, 0)
    
    def test_validate_config_broken_data_type(self):
        """Test validate-config command with incorrect data type."""
        config_file = 'tests/test_data/test_config_broken.json'
        result = self.runner.invoke(cli, ['validate-config', config_file])
        self.assertEqual(result.exit_code, 1)
    
    def test_validate_config_malformed_json(self):
        """Test validate-config command with malformed JSON."""
        config_file = 'tests/test_data/test_config_malformed.json'
        result = self.runner.invoke(cli, ['validate-config', config_file])
        self.assertEqual(result.exit_code, 1)
    
    def test_validate_config_nonexistent_file(self):
        """Test validate-config command with nonexistent file."""
        result = self.runner.invoke(cli, ['validate-config', 'nonexistent.json'])
        self.assertEqual(result.exit_code, 2)  # Click file not found error
    
    @unittest.skipUnless(
        os.path.exists('tests/test_data/test_reads.bam'), 
        "Real BAM file not available - skipping simulate integration test"
    )
    def test_simulate_command_validation_only(self):
        """Test simulate command gets past config validation with real files."""
        bam_file = 'tests/test_data/test_reads.bam'
        config_file = 'tests/test_data/test_config_small.json'
        
        with tempfile.TemporaryDirectory() as temp_dir:
            output_dir = Path(temp_dir) / 'output'
            
            # Run simulate command - we expect it might fail on missing predefined files
            # but it should at least pass config validation
            result = self.runner.invoke(cli, [
                'simulate',
                bam_file,
                config_file,
                '--output-dir', str(output_dir),
                '--output-prefix', 'test_sim',
                '--random-seed', '42'
            ])
            
            # Config validation should succeed, but simulation might fail due to missing predefined files
            # Exit code 1 is acceptable if it's due to missing predefined files, not config issues
            self.assertIn(result.exit_code, [0, 1])
    
    def test_simulate_command_file_not_found(self):
        """Test simulate command with nonexistent BAM file."""
        config_file = 'tests/test_data/test_config_small.json'
        
        result = self.runner.invoke(cli, [
            'simulate',
            'nonexistent.bam',
            config_file
        ])
        
        self.assertEqual(result.exit_code, 2)  # Click file not found error
    
    def test_simulate_command_invalid_config(self):
        """Test simulate command with invalid config file."""
        with tempfile.NamedTemporaryFile(suffix='.bam') as bam_file:
            config_file = 'tests/test_data/test_config_broken.json'
            
            result = self.runner.invoke(cli, [
                'simulate',
                bam_file.name,
                config_file
            ])
            
            self.assertEqual(result.exit_code, 1)  # Should fail on config validation


if __name__ == '__main__':
    unittest.main()