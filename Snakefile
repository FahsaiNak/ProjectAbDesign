rule all:
    input: 
        'Datasets/AbAg_processed.csv'

rule gen_pdb90:
    input:
        'Datasets/somePDB.csv'
    output:
        'gen_pdb90.temp'
    run:
        shell('cd run && bash PDB90.sh')
        shell('touch gen_pdb90.temp')

rule gen_cdrfrag:
    input:
        'gen_pdb90.temp' # It should be something like 'Datasets/chothia_pdb_files/{x}.pdb', but I don't know how to deal with wildcards.
    output:
        'gen_cdrfrag.temp'
    run:
        shell('cd run && bash CDR_fragment_database.sh')
        shell('touch gen_cdrfrag.temp')

rule target_prep:
    input: 
        'gen_pdb90.temp'
    output:
        'target_prep.temp'
    run:
        shell('cd run && bash Prep_targets_AbAg.sh')
        shell('touch target_prep.temp')

rule query_prep:
    input:
        'gen_cdrfrag.temp'
    output:
        'query_prep.temp'
    run:
        shell('cd run && bash Prep_queries_AbAg.sh')
        shell('touch query_prep.temp')

rule find_Ablike:
    input:
        'query_prep.temp', 'target_prep.temp'
    output: 
        'find_Ablike.temp'
    run:
        shell('cd run && bash Ablike.sh')
        shell('cd run && bash get_matchInfo.sh')
        shell('touch find_Ablike.temp')

rule find_Aglike:
    input: 
        'find_Ablike.temp'
    output: 
        'Datasets/AbAg.pkl'
    run: 
        shell('cd run && bash Aglike.sh')
        shell('rm find_Ablike.temp target_prep.temp query_prep.temp gen_pdb90.temp')

rule get_AbAg:
    input:
        'Datasets/AbAg.pkl'
    output:
        'Datasets/AbAg_processed.csv'
    run:
        shell('cd run && bash get_AbAg.sh')
