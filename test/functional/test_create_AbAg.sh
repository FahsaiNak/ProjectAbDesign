test -e ssshtest || wget -q https://raw.githubusercontent.com/ryanlayer/ssshtest/master/ssshtest
. ssshtest

run collect_Ablike_7js3 python ../../src/collect_Ablike.py --pdb 7js3 --data_dir ../Datasets/AbAg/
assert_equal ../Datasets/AbAg/7js3.pkl $( ls ../Datasets/AbAg/7js3.pkl )
assert_exit_code 0

run collect_Ablike_4d47 python ../../src/collect_Ablike.py --pdb 4d47 --data_dir ../Datasets/AbAg/
assert_equal ../Datasets/AbAg/7js3.pkl $( ls ../Datasets/AbAg/7js3.pkl )
assert_exit_code 0

run collect_Ablike_err python ../../src/collect_Ablike.py --pdb err1 --data_dir ../Datasets/AbAg/
assert_exit_code 1

run find_Aglike_7js3 python ../../src/find_Aglike.py --file ../Datasets/AbAg/7js3.pkl
assert_equal ../Datasets/AbAg/7js3.AbAg.pkl $( ls ../Datasets/AbAg/7js3.AbAg.pkl )

run find_Aglike_4d47 python ../../src/find_Aglike.py --file ../Datasets/AbAg/4d47.pkl
assert_equal ../Datasets/AbAg/4d47.AbAg.pkl $( ls ../Datasets/AbAg/4d47.AbAg.pkl )

run find_Aglike_err python ../../src/find_Aglike.py --file ../Datasets/AbAg/err1.pkl
assert_exit_code 1

run create_AbAg python ../../src/create_AbAg.py --dir ../Datasets/AbAg
assert_equal ../Datasets/AbAg.pkl $( ls ../Datasets/AbAg.pkl )

run create_AbAg_err python ../../src/create_AbAg.py --dir ../Datasets/AbAg_err
assert_exit_code 1
