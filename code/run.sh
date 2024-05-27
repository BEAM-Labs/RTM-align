filex=4v9k_AV.pdb
filey=4v9m_AV.pdb

./RTMalign ./${filex} ./${filey} -f

python parse_pdb.py --pdb_file_x ${filex} --pdb_file_y ${filey}