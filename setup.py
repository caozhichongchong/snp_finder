from setuptools import setup

setup(
    name="snp_finder",
    packages=['snp_finder'],
    version="0.0",
    description="identify and summarize SNPs in clonal populations",
    author='Anni Zhang',
    author_email='anniz44@mit.edu',
    url='https://github.com/caozhichongchong/snp_finder',
    keywords=['population genetics', 'WGS', 'SNPs', 'selection', 'parallel evolution'],
    license='MIT',
    install_requires=[
    'biopython',
    'argparse',
    'glob2',
    'statistics',
    'copy',
    'itertools',
    'random',
    'datetime'
    ],
    include_package_data=True,
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    package_dir={'snp_finder': 'snp_finder'},
    package_data={'snp_finder': ['bin/*','*.py']},
    entry_points={'console_scripts': ['snp_finder = snp_finder.__main__:main']},
    classifiers=[
        'Programming Language :: Python :: 3'
    ]
)
