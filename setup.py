from setuptools import setup, find_packages

setup(name='d3flux',
      version='0.2.5',
      description='A d3.js-based metabolic visualization tool for cobrapy',
      url='https://github.com/pstjohn/d3flux',
      download_url='https://github.com/pstjohn/d3flux',
      author='Peter St. John',
      author_email='peter.stjohn@nrel.gov',
      license='MIT',
      packages=find_packages(),
      install_requires=['pandas', 'cobra', 'jinja2', 'ipython', 'csscompressor'],
      package_data={'d3flux': ['templates/*']},
      )
