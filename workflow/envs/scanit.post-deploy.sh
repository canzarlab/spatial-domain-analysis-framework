if [ ! -d workflow/scripts/SCAN-IT ] ; then
    git clone https://github.com/zcang/SCAN-IT workflow/scripts/SCAN-IT
fi
pip install -e workflow/scripts/SCAN-IT --force-reinstall
