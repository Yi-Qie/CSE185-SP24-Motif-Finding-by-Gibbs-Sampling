from setuptools import setup

setup(
    name='gsmf',
    version='0.1',
    py_modules=['gsmf'],
    entry_points={
        'console_scripts': [
            'gsmf = gsmf:main',
        ],
    },
)
