from setuptools import find_packages, setup

requirements = """
  numpy
  pandas
  scipy
  statsmodels
  scikit-learn
  progressbar2
  scipy
  requests
""".split()

setup(
    name='func_e',
    packages=find_packages(),
    url='https://systemsgenetics.github.io/FUNC-E/',
    version='2.0.1',
    description='FUNC-E is a python library and script for functional enrichment of gene lists. It follows a similar approach to that of DAVID (https://david.ncifcrf.gov/) in that it performs enrichment analysis using a Fisher\'s test but then clusters enriched annotations using Kappa Statistics. FUNC-E allows the user to provide their own annotation lists. It is fully executable on a UNIX command-line.',
    author='Ficklin and Feltus computational Labs (Washington State University & Clemson University)',
    license='GNU General Public License v3.0',
    python_requires='>=3.6',
    install_requires=requirements,
    tests_require=['pytest'],
    entry_points={'console_scripts': [
        'FUNC-E = func_e.cmd:func_e',
        'FUNC-E-terms = func_e.cmd:getTerms'
    ]},
)
