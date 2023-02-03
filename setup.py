# -*- coding: utf-8 -*-
"""
Created on Fri June 6 3 16:40:44 2022

@author: danielgodinez
"""
from setuptools import setup, find_packages, Extension


setup(
    name="AstroDynamics",
    version="0.1",
    author="Daniel Godines",
    author_email="danielgodinez123@gmail.com",
    description="Numerical Methods for Hydrodynamical Systems",
    license='GPL-3.0',
    url = "https://github.com/Professor-G/AstroDynamics",
    classifiers=[
		'Development Status :: 5 - Production/Stable',
		'Intended Audience :: Developers',
		'Topic :: Software Development :: Build Tools',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3)',
		'Programming Language :: Python :: 3',	   
],
    packages=find_packages('.'),
    install_requires = ['numpy', 'matplotlib', 'progress'],
    python_requires='>=3.7,<4',
    include_package_data=True,
    test_suite="nose.collector",
)
