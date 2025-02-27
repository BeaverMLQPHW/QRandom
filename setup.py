from setuptools import setup, find_packages

setup(
    name="qrandom",
    version="0.1.0",
    author="Maksim",
    author_email="qragacontact@gmail.com",
    description="A quantum random number generator based on Qiskit and other services that provide their own API. Contains almost all the functionality of the random library in Python.",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/BeaverMLQPHW/QRandom",
    packages=find_packages(),
    install_requires=[
        "qiskit",
        "qiskit-aer",
        "qiskit-ibm-runtime",
        "requests",
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: Apache 2.0",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.7",
)