if [ ! -d workflow/scripts/PROST ] ; then
    git clone https://github.com/Tang-Lab-super/PROST.git workflow/scripts/PROST ;
fi
pip install --no-deps -e workflow/scripts/PROST
