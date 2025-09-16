#!/bin/bash
# Activation script for ssiamb development environment
# Usage: source activate_env.sh

echo "Activating ssiamb conda environment..."
conda activate ./env

echo "Environment activated!"
echo "Python version: $(python --version)"
echo "Environment location: $(which python)"
echo ""
echo "Core packages installed:"
echo "  - typer $(python -c 'import typer; print(typer.__version__)')"
echo "  - pandas $(python -c 'import pandas; print(pandas.__version__)')"
echo "  - numpy $(python -c 'import numpy; print(numpy.__version__)')"
echo "  - pysam $(python -c 'import pysam; print(pysam.__version__)')"
echo ""
echo "External bioinformatics tools:"
echo "  - minimap2 $(minimap2 --version 2>/dev/null || echo 'not found')"
echo "  - bwa-mem2 $(bwa-mem2 version 2>/dev/null | tail -1 || echo 'not found')"
echo "  - bcftools $(bcftools --version 2>/dev/null | head -1 | cut -d' ' -f2 || echo 'not found')"
echo "  - samtools $(samtools --version 2>/dev/null | head -1 | cut -d' ' -f2 || echo 'not found')"
echo "  - mosdepth $(mosdepth --version 2>/dev/null || echo 'not found')"
echo "  - bbtools/callvariants.sh $(test -f $(which callvariants.sh 2>/dev/null) && echo 'available' || echo 'not found')"
echo ""
echo "Development environment ready for ssiamb development!"