from os.path import join, dirname, abspath

ngsdb_dir = dirname(dirname(abspath(__file__)))

ngsdb_code = join(ngsdb_dir, 'ngsdb_code')
ngsdb_data = join(ngsdb_dir, 'ngsdb_data')

