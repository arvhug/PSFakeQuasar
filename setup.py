from setuptools import setup, find_packages

setup(
    name="PSFakeQuasar",
    version="0.1.0",
    description="Tools for PSF modelling and photometry",
    author="Arvind Hughes",
    packages=find_packages(),
    install_requires=[
        'numpy>=1.20',
        'astropy>=5.0',
        'photutils>=1.5',
        'scipy>=1.7'
    ],
    python_requires='>=3.7',
    keywords="astronomy psf quasars euclid simulation",
    url="https://github.com/arvhug/PSFakeQuasar",
    license="GPLv3", 
    classifiers=[
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
    ],
)