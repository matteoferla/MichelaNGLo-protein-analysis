import setuptools
#from setuptools.command.install import install

#with open("README.md", "r") as fh:
#    long_description = fh.read()
long_description = 'See github readme.md'

setuptools.setup(
    name="michelanglo_protein",
    version="0.4.1",
    author="Matteo Ferla",
    author_email="matteo@well.ox.ac.uk",
    description="protein module for Michelanglo",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/matteoferla/MichelaNGLo-protein-module",
    packages=setuptools.find_packages(),
    include_package_data=True,
    package_data={'michelanglo_protein.analyse.pyrosetta_modifier.params': ['*.params'],},
    classifiers=[
        'Development Status :: 5 - Production/Stable', #
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Chemistry',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
        "Operating System :: OS Independent",
    ],
    install_requires=['Bio', 'requests_ftp', 'typing-extensions', 'pyrosetta-help', 'ConsurfDB-client-API']
)

'''
class PostInstallCommand(install):
    """Post-installation for installation mode."""
    def run(self):
        ##
        install.run(self)

'''
