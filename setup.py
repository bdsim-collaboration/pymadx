from setuptools import setup, find_packages

try:
    import pypandoc
    long_description = pypandoc.convert_file("README.md", "rst")
except:
    print ("Warning: pypandoc module not found, could not convert Markdown to reStructuredText." )
    long_description = ""


setup(
    name='pymadx',
    version='1.8.2',
    packages=find_packages(exclude=["docs", "tests", "obsolete"]),
    install_requires=["matplotlib>=1.7.1",
                      "numpy >= 1.4",
                      "pytransport"],
    setup_requires=["pytest-runner"],
    tests_require=["pytest"],
    python_requires=">=3.7.*",

    author='JAI@RHUL',
    author_email='laurie.nevay@rhul.ac.uk',
    description="Write MADX models and load MADX output.",
    long_description=long_description,
    url='https://bitbucket.org/jairhul/pymadx/',
    license='GPL3',
    keywords=['madx', 'accelerator', 'twiss', 'ptc'],
)
