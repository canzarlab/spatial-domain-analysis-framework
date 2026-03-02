if [ ! -d workflow/scripts/STAGATE ] ; then
    git clone https://github.com/QIFEIDKN/STAGATE workflow/scripts/STAGATE ;
fi
pip install -e workflow/scripts/STAGATE ;
pip install --no-deps scikit-misc
