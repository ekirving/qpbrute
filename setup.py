import setuptools

with open("README.md", "r") as f:
    long_description = f.read()

with open("requirements.txt") as f:
    required = [x for x in f.read().splitlines() if not x.startswith("#")]

setuptools.setup(
    name="qpbrute",
    version="0.2",
    author="Evan K. Irving-Pease",
    author_email="evan.irvingpease@gmail.com",
    description="Heuristic algorithm for automatically fitting admixture graphs",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ekirving/qpbrute",
    package_dir={"qpbrute": "qpbrute", "qpbrute.rscript": "qpbrute/rscript"},
    packages=["qpbrute", "qpbrute.rscript"],
    include_package_data=True,
    entry_points={
        "console_scripts": [
            "qpBrute = qpbrute.qpbrute:qpbrute",
            "qpBayes = qpbrute.qpbayes:qpbayes",
            # "qpCluster = qpbrute.qpcluster:qpcluster",
        ]
    },
    classifiers=["Programming Language :: Python :: 3", "License :: MIT License"],
    license="MIT",
    install_requires=required,
)
