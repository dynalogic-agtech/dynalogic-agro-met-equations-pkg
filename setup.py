from setuptools import setup, find_packages

setup(
    name='agro-met-equations-dynalogic',
    version='0.0.3',
    description=(
        'Library for agrometeorological calculation like evapotranspiration, water balance, degree days, etc'
    ),
    long_description=open('README.md').read(),
    long_description_content_type="text/markdown",
    author='Bruno Ducraux',
    author_email='bducraux@dynalogic.net',
    license='BSD 3-Clause',
    url='https://github.com/dynalogic-agtech/dynalogic-agro-met-equations-pkg/',
    packages=find_packages(),
    include_package_data=True,
    test_suite='tests',
    classifiers=[
        "Programming Language :: Python :: 3",
        'License :: OSI Approved :: BSD License',
        'Operating System :: OS Independent',
        'Intended Audience :: Science/Research',
        'Development Status :: 3 - Alpha',
        'Topic :: Scientific/Engineering',
        'Topic :: Scientific/Engineering :: Atmospheric Science',
    ],
    python_requires='>=3.6',
)
