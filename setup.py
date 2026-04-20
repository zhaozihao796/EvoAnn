from setuptools import setup, find_packages

setup(
    name="EvoAnn",
    version="1.0.0",
    description="A high-performance, GFF-based variant annotation and sequence processing suite tailored for evolutionary genomics.",
    packages=find_packages(),
    python_requires=">=3.10",
    install_requires=[
        "numpy",
        "biopython",
        "cyvcf2",
    ],
    entry_points={
        "console_scripts": [
            "evoann=src.main:main",
        ],
    },
)