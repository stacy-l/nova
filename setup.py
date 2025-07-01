from setuptools import setup, find_packages

setup(
    name="nova",
    version="0.1.0",
    description="de novo variant simulator",
    long_description="de novo variant simulator",
    author="voidshapes",
    packages=find_packages(where="src"),
    package_dir={"": "src"},
    python_requires=">=3.10",
    install_requires=[
        "setuptools==80.9.0", # note: flagged for deprecation, dependency of pyranges
        "biopython==1.85",
        "pysam==0.23.3",
        "pandas==2.3.0",
        "numpy==2.3.1",
        "pyranges==0.1.4",
        "click==8.2.1",
        "matplotlib==3.10.3",
        "seaborn==0.13.2",
    ],
    extras_require={
        "test": [
            "pytest==8.4.1",
            "pytest-cov==6.2.1",
        ],
    },
    entry_points={
        "console_scripts": [
            "nova=nova.cli:main",
        ],
    }
)