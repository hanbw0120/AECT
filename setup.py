from setuptools import setup

setup(name='AECT',
      version='0.9',
      description='Autoencoder on Cell-free DNA TSS coverage profile',
      packages=['ae'],
      scripts=['AECT.py'],

      author='Han Bowei',
      author_email='hanbw0120@foxmail.com',
      url='https://github.com/hanbw0120/AECT', 

      install_requires=['numpy>=1.7',
                        'keras>=2.0.8',
                        'pandas',
                        'tensorflow==1.15.4'
                        ],
      python_requires='>3.6.0',
     )
