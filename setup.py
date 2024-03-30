from setuptools import find_packages, setup

print(find_packages(include=["memdof", "memdof.*"]))
setup(
    name="memdof",
    version="0.1.1",
    description="A python package to determine the degrees of freedom for a single lipid membrane.",
    url="https://github.com/NoxelS/memdof",
    author="Noel Schwabenland",
    author_email="noel@lusi.uni-sb.de",
    license="MIT",
    packages=find_packages(include=["memdof", "memdof.*"]),
    install_requires=["numpy>=1.21.5", "biopython>=1.81"],
    classifiers=[
        "Development Status :: 1 - Planning",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Physics",
        "License :: OSI Approved :: MIT License",
        "Operating System :: POSIX :: Linux",
        "Programming Language :: Python :: 3.9",
        "Typing :: Typed",
    ],
)
