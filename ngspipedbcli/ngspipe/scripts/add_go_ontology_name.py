import sys
import re

obofile = sys.argv[1]
annofile = sys.argv[2]

go2name = dict()

with open(obofile, 'r') as f:
    for line in f:
        line = line.rstrip()
        matObj = re.match(r'^\[(\S+)\]$', line)
        if matObj:
            if matObj.group(1) == "Term":
                id_line = f.readline().rstrip()
                name_line = f.readline().rstrip()
                namespace_line = f.readline().rstrip()
                matObj_id = re.match(r'^id: (GO:\d+)$', id_line)
                matObj_name = re.match(r'^name: (.*)$', name_line)
                matObj_namespace = re.match(r'^namespace: (.*)$', namespace_line)
                #continue
                if matObj_id:
                    id = matObj_id.group(1)
                else:
                    sys.stderr.write("id error\n")
                    sys.exit(-1)
                if matObj_name:
                    name = matObj_name.group(1)
                if matObj_namespace:
                    namespace = matObj_namespace.group(1)
                else:
                    sys.stderr.write("name error\n")
                    sys.exit(-1)
                go2name[id] = [name, namespace]

with open(annofile, 'r') as f:
    for line in f:
        line = line.rstrip()
        if '\tGO:' in line:
            gene, go = line.split('\t')
            if go in go2name.keys():
                print(gene, go, go2name[go][0], go2name[go][1], sep='\t')
            
