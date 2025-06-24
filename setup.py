from setuptools import setup, find_packages

setup(
    name="nova",
    version="0.1.0",
    description="De novo variant insertion simulator for testing structural variant detection",
    author="voidshapes",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    python_requires=">=3.8",
    install_requires=[
        "pysam>=0.19.0",
        "numpy>=1.20.0",
        "biopython>=1.79",
        "click>=8.0.0",
    ],
    extras_require={
        "test": [
            "pytest>=6.0.0",
            "pytest-cov>=2.12.0",
        ],
    },
    entry_points={
        "console_scripts": [
            "nova=nova.cli:main",
        ],
    }
)