# -*- coding:utf-8 -*-
from setuptools import setup
# from setuptools import find_packages
# from OHDLF_pro.__main__ import __VERSION__
VERSION = "1.0.0"

def readme():
    with open('README.md', encoding='utf-8') as f:
        content = f.read()
    return content


setup(
    name='ohdlf-pro',
    version=VERSION,
    description='OHDLF-pro: a pipeline designed to filter and address orthologous gene heterogeneity, duplication, and loss',
    long_description=readme(),  # 长文描述
    long_description_content_type='text/markdown',  # 长文描述的文本格式
    author=['caijunhao'],
    author_email='cjunhao416@gmail.com',
    url='https://github.com/qunjie-zhang/OHDLF-pro',
    keywords=["Phylogenetic",'bioinformatics','pipeline'],
    packages=["OHDLF_pro", "OHDLF_pro.utils"],
    # packages=find_packages(),
    python_requires='>=3.5',
    classifiers=[
            'Development Status :: 5 - Production/Stable',
            'License :: OSI Approved :: Apache Software License',
            'Operating System :: OS Independent',
            'Programming Language :: Python :: 3.5',
            'Programming Language :: Python :: 3.6',
            'Programming Language :: Python :: 3.7',
            'Programming Language :: Python :: 3.8',
            'Programming Language :: Python :: 3.9',
            'Programming Language :: Python :: 3.10',
            'Programming Language :: Python :: 3.11',
        ],
    install_requires=[
        "biopython",
        "tqdm",
        "treeswift",
        "importlib-metadata; python_version < '3.8'",
    ],
    entry_points={
        'console_scripts': [
            'OHDLF-pro = OHDLF_pro.main:run',# 格式: 命令名 = 包名.文件名:函数名
        ]
    },
    include_package_data = True
)