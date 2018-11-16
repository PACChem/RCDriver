""" Install RCDriver
"""
from distutils.core import setup


setup(name="rc_driver",
      scripts=["rc_driver.py"],
      package_dir={'me_parser': 'external/me_parser/me_parser'},
      packages=['rcd', 'me_parser'])
