from setuptools import setup, find_packages

# Function to read requirements from requirements.txt
def parse_requirements(filename):
    with open(filename, 'r') as f:
        return f.read().splitlines()

setup(
    name="CompactObject",
    version="1.9.8",
    packages=find_packages(),
    install_requires=parse_requirements('requirements.txt'),  # Automatically adds requirements
    # other setup parameters
)