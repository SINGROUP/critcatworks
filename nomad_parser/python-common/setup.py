"""
This is a setup script for those that are familiar with setuptooll. You can
also just setup the paths and necessary dependencies in any other way you wish.
"""
from setuptools import setup, find_packages


#===============================================================================
def main():
    # Start package setup
    setup(
        name="nomadcore",
        version="0.1",
        description="Tools for NOMAD parser development.",
        package_dir={'': 'common/python'},
        packages=find_packages("common/python"),
        package_data={
            'nomadcore.unit_conversion': ['*.txt'],
            'nomadcore.md_data_access': ['test/*'],
            'nomadcore.metainfo_storage': ['*.txt'],
        }
    )

# Run main function by default
if __name__ == "__main__":
    main()
