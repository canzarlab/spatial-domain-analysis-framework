if [ ! -d workflow/scripts/CCST ] ; then
    git clone https://github.com/xiaoyeye/CCST workflow/scripts/CCST
fi
if [ ! -d workflow/scripts/SEDR ] ; then
    git clone https://github.com/JinmiaoChenLab/SEDR workflow/scripts/SEDR
fi
pip install spaceflow --no-deps
