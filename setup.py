from setuptools import setup, find_packages

setup(
    name="q2-epitope",
    version='1',
    packages=find_packages(),
    package_data={},
    author="",
    author_email="",
    description="",
    license="",
    url="",
    entry_points={
        "qiime2.plugins": ["q2-epitope=q2_epitope.plugin_setup:plugin"]
    },
    zip_safe=False,
)
