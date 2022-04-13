from setuptools import setup, find_packages


with open('README.md') as readme_file:
    readme = readme_file.read()


requirements = [
    'setuptools',
    'tqdm',
    'click',
    'pandas>=1.3',
    'pyranges>=0.0.71',
    'matplotlib',
]

setup_requirements = ['pytest-runner', ]

test_requirements = ['pytest']

setup(
    name='cerberus',
    version='0.0.idk',

    author="M. Hasan Çelik, Fairlie Reese",
    author_email='muhammedhasancelik@gmail.com, fairlie.reese@gmail.com',
    url='https://github.com/mortazavilab/cerberus',

    keywords=['genomics', 'long read RNA-seq'],
    description="TODO",

    classifiers=[
        'Natural Language :: English',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9'
    ],
    license="MIT license",
    long_description=readme + '\n',
    long_description_content_type='text/markdown',

    install_requires=requirements,
    setup_requires=setup_requirements,

    entry_points='''
        [console_scripts]
        cerberus=cerberus.main:cli
    ''',
    packages=find_packages(include=['cerberus*']),
    include_package_data=True,

    test_suite='tests',
    tests_require=test_requirements,
)
