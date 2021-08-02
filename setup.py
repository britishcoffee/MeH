import setup, find_packages

VERSION = '0.0.1' 
DESCRIPTION = 'A bioinformatic tool for profiling methylation heterogeneity genomewide'
LONG_DESCRIPTION = ''

with open("README.md", "r") as fh:
    long_description = fh.read()

# Setting up
setup(
        name="MeHscr", 
        version=VERSION,
        author="Ya-Ting Sabrina Chang",
        author_email="<ytchang.sabrina@gmail.com>",
        url="https://github.com/britishcoffee/Methylationhet",
	license='MIT',
        packages=find_packages(),
        description=DESCRIPTION,
        long_description=LONG_DESCRIPTION,
        install_requires=['pysam==0.16.0.1','numpy==1.16.6','pandas==0.24.2','joblib'], # add any additional packages that 
        # needs to be installed along with your package. Eg: 'caer'
        
        keywords=['python', 'first package'],
        classifiers= [
            "Intended Audience :: Science/Research",
            "Programming Language :: Python :: 2",
            "Programming Language :: Python :: 3",
	    "License :: OSI Approved :: MIT License",
            "Operating System :: MacOS :: MacOS X",
            "Operating System :: POSIX :: Linux",
	
        ]
)
