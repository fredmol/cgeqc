import os
import logging
import subprocess
import sys

from cgeqc.qc_report import create_qc_report
from cgeqc.qc_config import KMA_DEFAULTS

class TrimRunner():
    def __init__(self, input_file, output_dir, name, pipeline_type="bacterial", parameters=None):
        self.logger = logging.getLogger(__name__)
        self.logger.setLevel(logging.INFO)
        self.input_file = input_file
        self.output_dir = output_dir
        self.trimmed_name = f"{name}"
        self.trimmed_output_path = os.path.join(self.output_dir, self.trimmed_name)
        self.pipeline_type = pipeline_type
        self.parameters = parameters if parameters is not None else KMA_DEFAULTS

    def run(self):
        """Runs KMA trim on input file"""
        # Make sure output directory exists
        os.makedirs(self.output_dir, exist_ok=True)
        
        trim_cmd = (
            f"kma trim "
            f"-ml {self.parameters['min_length']} "
            f"-xl {self.parameters['max_length']} "
            f"-mp {self.parameters['min_phred']} "
            f"-mi {self.parameters['min_internal_phred']} "
            f"-eq {self.parameters['min_average_quality']} "
            f"-5p {self.parameters['trim_5_prime']} "
            f"-3p {self.parameters['trim_3_prime']} "
            f"-qc "
            f"-i {self.input_file} "
            f"-o {self.trimmed_output_path}"
        )
        
        self.logger.info(f"Running KMA trim with command: {trim_cmd}")
        
        # Execute command
        ret = os.system(trim_cmd)
        if ret != 0:
            error_msg = f"Error: KMA trim failed with return code {ret}"
            print(error_msg)
            sys.exit(error_msg)
            
        expected_output = f"{self.trimmed_output_path}.fq"
        if not os.path.exists(expected_output):
            error_msg = f"Error: KMA trim failed to create output file {expected_output}"
            print(error_msg)
            sys.exit(error_msg)
        
        # Generate QC report if json exists
        json_output = f"{self.trimmed_output_path}.json"
        if os.path.exists(json_output):
            try:
                print(f"Generating QC report...")
                qc_report_path = create_qc_report(
                    json_output, 
                    self.output_dir, 
                    self.trimmed_name,
                    self.pipeline_type,
                    self.parameters
                )
                self.logger.info(f"Generated QC report: {qc_report_path}")
                print(f"QC report generated: {qc_report_path}")
            except Exception as e:
                error_msg = f"Failed to generate QC report: {str(e)}"
                self.logger.error(error_msg)
                print(f"WARNING: {error_msg}")
        else:
            print(f"WARNING: KMA did not generate a JSON file. QC report cannot be created.")
        
        return expected_output