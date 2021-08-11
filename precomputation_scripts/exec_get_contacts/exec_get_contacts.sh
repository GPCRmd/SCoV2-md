#!/bin/bash

python -u fix_missing_chainindex.py
python -u  covid_get_contacts_dynfreq_ori.py   
python -u covid_flareplot_from_getcontacts_ori.py
