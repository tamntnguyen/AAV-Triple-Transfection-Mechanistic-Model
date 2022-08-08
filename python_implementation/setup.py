# pip install --editable .
# python setup.py bdist_wheel
from setuptools import find_packages, setup

setup(
    name="aav_ode",
    version="0.1.1",
    description="ODE model of AAV transfection plasmid production",
    long_description_content_type="text/markdown",
    package_dir={"": "src"},  # Packages are subfolders of "src" directory
    packages=find_packages(where="src"),
    package_data={"": ["data/*"]},
    python_requires="> 3.6",
    install_requires=["pandas", "numpy", "plotly", "click", "scipy"],
    classifiers=[
        "Topic :: Scientific/Engineering :: Gene Therapy",
        "Topic :: Scientific/Engineering :: Adeno-associated Virus",
        "Programming Language :: Python :: 3",
    ],
)
