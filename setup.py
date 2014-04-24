import ez_setup
ez_setup.use_setuptools()
from setuptools import setup

# from mpld3
def get_version(path):
    """Get the version info from the mpld3 package without importing it"""
    import ast

    with open(path) as init_file:
        module = ast.parse(init_file.read())

    version = (ast.literal_eval(node.value) for node in ast.walk(module)
               if isinstance(node, ast.Assign)
               and node.targets[0].id == "__version__")
    try:
        return next(version)
    except StopIteration:
        raise ValueError("version could not be located")


setup(name='bwameth',
      version=get_version("bwameth.py"),
      description="align BS-Seq reads with bwa mem",
      py_modules=['bwameth'],
      author="Brent Pedersen",
      author_email="bpederse@gmail.com",
      license="MIT",
      install_requires=["toolshed"],
      long_description=open('README.md').read(),
      classifiers=["Topic :: Scientific/Engineering :: Bio-Informatics"],
      scripts=['bwameth.py']
  )
