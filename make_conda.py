import yaml
import os
import sys

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)), '')] + sys.path

import cgeqc.version as version

data = {
    "package": {
        "name": "cgeqc",
        "version": version.__version__
    },
    "source": {
        "url": "https://bitbucket.org/genomicepidemiology/cgeqc/archive/refs/tags/{}.tar.gz".format(version.__version__),
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
            "python >=3.6"
        ]
    },
    "about": {
        "home": "https://bitbucket.org/genomicepidemiology/cgeqc",
        "summary": "CGE quality control tool",
        "license": "Apache-2.0"
    }
}

os.system('mkdir -p conda')
yaml_str = yaml.dump(data, sort_keys=False)

with open('conda/meta.yaml', 'w') as f:
    f.write(yaml_str)
