import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="deeppbs",
    version="0.0.1",
    author="Raktim Mitra, Jared Sagendorf, Jinsen Li, Yibei Jiang and Remo Rohs",
    author_email="jsagendorf@gmail.com",
    description="Geometry-invariant deep learning on protein-DNA structures for interpretable prediction of binding specificity across protein families",
    long_description=long_description,
    long_description_content_type="text/markdown",
    include_package_data=True,
    url="https://github.com/timkartar/DeepPBS",
    packages=['deeppbs', 'deeppbs.nn', 'deeppbs._data'],
    classifiers=[
	"Development Status :: 3 - Alpha",
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)
