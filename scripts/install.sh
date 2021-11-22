#!/usr/bin/env bash
current_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
. "$current_dir/info.sh" && \
python setup.py bdist_wheel && \
pip install --force dist/pydockrmsd-${version}-*.whl && \
printf "pydockrmsd " && \
python -c "import pydockrmsd; print(pydockrmsd.__version__, end='');" && \
echo " correctly installed";
cd examples && python3 crystal_example.py