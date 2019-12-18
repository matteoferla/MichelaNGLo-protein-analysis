import setuptools
##from setuptools.command.install import install

#with open("README.md", "r") as fh:
#    long_description = fh.read()
long_description = 'See github readme.md'

setuptools.setup(
    name="michelanglo_protein",
    version="0.4.0",
    author="Matteo Ferla",
    author_email="matteo@well.ox.ac.uk",
    description="protein module for Michelanglo",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/matteoferla/MichelaNGLo-protein-module",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ], install_requires=['Bio', 'requests_ftp'
    ]
)

'''
class PostInstallCommand(install):
    """Post-installation for installation mode."""
    def run(self):
        ###
        install.run(self)

'''
