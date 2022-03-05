#!/usr/bin/env bash
current_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
. "$current_dir/info.sh" && \
# python -m build . && \
python setup.py bdist_wheel && \
# python -m cibuildwheel --platform linux \
# pip install --force dist/pydockrmsd-${version}-py3-none-any.whl && \
pip install --force dist/pydockrmsd-${version}-cp38-cp38-linux_x86_64.whl && \
printf "pydockrmsd " && \
python -c "import pydockrmsd; print(pydockrmsd.__version__, end='');" && \
echo " correctly installed";
# rm -r build/
