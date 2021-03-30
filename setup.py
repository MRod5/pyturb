# https://github.com/pypa/sampleproject

from setuptools import setup, find_packages


setup(name='pyTurb', version='0.4.0', author='Marcos Rodriguez', description='Gas Turbine solver', url='https://github.com/MRod5/pyturb', python_requires='>=3.5', install_requires=['numpy', 'scipy'],
      packages=find_packages(where='src'), package_dir={'': 'src'}, include_package_data=True, zip_safe=False)

