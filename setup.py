import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="gptenrichment",
    version="0.0.1",
    author="Alexander Lachmann",
    author_email="alexander.lachmann@mssm.edu",
    description="Package for over-representation analysis using LLM.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/maayanlab/blitzgsea",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Apache Software License",
        "Operating System :: OS Independent",
    ],
    package_data={
        "gptenrichment": ["data/*"]
    },
    include_package_data=True,
    install_requires=[
        'pandas>=1.1.5',
        'numpy',
        'scikit-learn',
        'progress',
        'tqdm',
        'statsmodels'
    ],
    python_requires='>=3.6',
)
