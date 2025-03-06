"""
Configuration file for QC thresholds and parameters
"""

# Quality assessment thresholds for bacterial data
QC_THRESHOLDS = {
    "GOOD": {
        "min_quality": 15,
        "min_coverage": 50,
        "min_read_length": 3000,
    },
    "FAIR": {
        "min_quality": 12,
        "min_coverage": 20,
        "min_read_length": 2000,
    },
    "POOR": {
        "min_quality": 12,  # Less than this is poor
        "min_coverage": 10,  # Less than this is poor
        "min_read_length": 1000,  # Less than this is unusually short for ONT
    },
    "GC_CONTENT": {
        "min": 25, # bacterial genome range: 25 - 75%
        "max": 75
    }
}

# Quality assessment thresholds for viral data
QC_THRESHOLDS_VIRUS = {
    "GOOD": {
        "min_quality": 15,
        "min_bp_count": 1000000,  # 1 Mbp for good coverage of most viral genomes
    },
    "FAIR": {
        "min_quality": 12,
        "min_bp_count": 500000,  # 500 Kbp for adequate coverage
    },
    "POOR": {
        "min_quality": 12,  # Less than this is poor
        "min_bp_count": 100000,  # 100 Kbp - less than this is poor for viral analysis
    }
}

# Quality assessment thresholds for metagenomic data
QC_THRESHOLDS_META = {
    "GOOD": {
        "min_quality": 15,
        "min_bp_count": 10000000,  # 10 Mbp for good community representation
    },
    "FAIR": {
        "min_quality": 12,
        "min_bp_count": 5000000,  # 5 Mbp for some diversity coverage
    },
    "POOR": {
        "min_quality": 12,  # Less than this is poor
        "min_bp_count": 1000000,  # 1 Mbp - less than this is poor for metagenomic analysis
    }
}

# Default KMA parameters - these will be used if no custom parameters are provided
KMA_DEFAULTS = {
    "min_length": 16,
    "max_length": 2147483647,
    "min_phred": 20,
    "min_internal_phred": 0,
    "min_average_quality": 10,
    "trim_5_prime": 0,
    "trim_3_prime": 0
}

# Pipeline-specific trim parameters - kept for backward compatibility
# These should be moved to CGELabs in future versions
BACTERIAL_TRIM_DEFAULTS = {
    "min_length": 500,
    "max_length": 2147483647,
    "min_phred": 20,
    "min_internal_phred": 0,
    "min_average_quality": 10,
    "trim_5_prime": 0,
    "trim_3_prime": 0
}

VIRAL_TRIM_DEFAULTS = {
    "min_length": 100,
    "max_length": 500000, # Mpox is around 280k
    "min_phred": 20,
    "min_internal_phred": 0,
    "min_average_quality": 10,
    "trim_5_prime": 0,
    "trim_3_prime": 0
}

METAGENOMIC_TRIM_DEFAULTS = {
    "min_length": 500,
    "max_length": 2147483647,
    "min_phred": 20,
    "min_internal_phred": 0,
    "min_average_quality": 10,
    "trim_5_prime": 0,
    "trim_3_prime": 0
}

# For backward compatibility
TRIM_DEFAULTS = KMA_DEFAULTS

def get_thresholds(pipeline_type="bacterial"):
    """Get the appropriate QC thresholds based on pipeline type."""
    if pipeline_type == "viral":
        return QC_THRESHOLDS_VIRUS
    elif pipeline_type == "metagenomic":
        return QC_THRESHOLDS_META
    else:
        return QC_THRESHOLDS