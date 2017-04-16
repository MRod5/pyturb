from setuptools import setup, find_packages

setup(name='pyTurb', version='0.1a0', python_requires='>=3.5', install_requires=['numpy', 'scipy'],
      packages=find_packages(''), package_dir={'': 'src'}, include_package_data=True, zip_safe=False)

# http://www.siafoo.net/article/77