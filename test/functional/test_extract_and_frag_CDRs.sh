test -e ssshtest || wget -q https://raw.githubusercontent.com/ryanlayer/ssshtest/master/ssshtest
. ssshtest

CDR_file_to_create='../Datasets/CDR_pdb_files/1f58_H1.pdb'
run basic_extract python ../../src/extract_CDRs_from_PDB.py
                    
assert_in_stdout 0
assert_exit_code 0
assert_equal $CDR_file_to_create $( ls $CDR_file_to_create )


fragment_file_to_create='../Datasets/CDR_fragments/1f58_H1_frag_26_29.pdb'
run basic_extract python ../../src/fragment_CDRs.py
                    
assert_in_stdout 0
assert_exit_code 0
assert_equal $fragment_file_to_create $( ls $fragment_file_to_create )