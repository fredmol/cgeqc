import yaml
import os
import sys
import re

# Get the absolute path to the directory containing this script
script_dir = os.path.dirname(os.path.realpath(__file__))
# Add the script directory to sys.path
sys.path.insert(0, script_dir)

# Read version from version.py
version_path = os.path.join(script_dir, 'cgeqc', 'version.py')
with open(version_path, 'r') as f:
    version_file = f.read()
    version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]", version_file, re.M)
    if version_match:
        __version__ = version_match.group(1)
    else:
        raise RuntimeError("Unable to find version string in cgeqc/version.py")

# Ensure conda directory exists
os.makedirs('conda', exist_ok=True)

data = {
    "package": {
        "name": "cgeqc",
        "version": __version__
    },
    "source": {
        "url": f"https://bitbucket.org/genomicepidemiology/cgeqc/get/{__version__}.tar.gz",
    },
    "build": {
        "number": 0,
        "noarch": "python",
        "script": "{{ PYTHON }} -m pip install . --no-deps --ignore-installed -vvv"
    },
    "requirements": {
        "host": [
            "python >=3.6",
            "pip"
        ],
        "run": [
            "python >=3.6",
            "kma >=1.4.9",
            "weasyprint >=63.0",
            "jinja2 >=3.1.4",
            "matplotlib >=3.9.2",
            "scipy >=1.10.0",
            "numpy >=1.24.0"
        ]
    },
    "about": {
        "home": "https://bitbucket.org/genomicepidemiology/cgeqc",
        "summary": "CGE Quality Control Tool for sequence data",
        "license": "Apache-2.0"
    }
}

# Convert the data to YAML and print it
yaml_str = yaml.dump(data, sort_keys=False)

with open('conda/meta.yaml', 'w') as f:
    f.write(yaml_str)

print(f"Created conda/meta.yaml for version {__version__}")