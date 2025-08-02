from setuptools import setup, find_packages

setup(
    name="metagenomics_analysis",
    version="0.1.0",
    packages=find_packages(),
    install_requires=[
        "argparse",
        "staramr",
        "bioservices",
        "matplotlib",
        "pandas",
        "seaborn",
        "xlsxwriter",
        "matplotlib-venn",
    ],
    scripts=[
        "metagenomics_pipeline.py",
        "genomic_pipeline.py",
        "visualize.py",
        "main.py",
    ],
    author="Your Name",
    author_email="your.email@example.com",
    description="A pipeline for metagenomics analysis, variant calling, and visualization",
    long_description=open("README.md").read() if os.path.exists("README.md") else "",
    long_description_content_type="text/markdown",
    url="https://github.com/yourusername/metagenomics_analysis",
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.10",
)