from setuptools import setup

setup(name='d3flux',
      version='0.1',
      description='A d3.js-based metabolic visualization tool for cobrapy',
      url='http://github.com/pstjohn/d3flux',
      author='Peter St. John',
      author_email='peter.stjohn@nrel.gov',
      license='MIT',
      packages=['d3flux'],
      install_requires=['pandas', 'cobra', 'jinja2', 'ipython'],
      )
