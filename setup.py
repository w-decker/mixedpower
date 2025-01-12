from setuptools import setup, find_packages

setup(
    name="mixedpower",
    version="1.0.0",
    author="Will Decker",
    author_email="will.decker@gatech.edu",
    description="Power for linear mixed effects models in Python",
    url="https://github.com/w-decker/mixedpower",
    packages=find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering"
    ]

)