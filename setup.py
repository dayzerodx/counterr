from setuptools import setup

setup(
    name="counterr",
    version="0.1",
    packages=["counterr"],
    entry_points={
        'console_scripts': ['counterr=counterr.counterr:main'],
    },

    # Metadata
    author="DZD",
    author_email="jae@dayzerodiagnostics.com",
    description="A sequencing error profiler"
)