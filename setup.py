import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()


setuptools.setup(
    name='cgr_view',
    version='0.0.1dev2',
    url='https://github.com/TeamMacLean/cgr_view',
    packages=setuptools.find_packages(),
    long_description_content_type="text/markdown",
    license='LICENSE.txt',
    author='Dan MacLean',
    author_email='dan.maclean@tsl.ac.uk',
    description='A package for rendering Chaos Game Representations of DNA sequences',
    python_requires='>=3.7',
    install_requires = [
        "numpy >= 1.17",
        "matplotlib >= 3.1.0",
        "scipy >= 1.3.1",
        "pillow >= 6.2.1",
        "scikit-image >= 0.15.0",
        "biopython >= 1.7"
    ]
)
