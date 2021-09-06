# make network object and load file
from clustergrammer import Network
net = Network()

b = "cluster.txt"
d= "cluster.json"

net.load_file(b)

# calculate clustering using default parameters
net.cluster()

# save visualization JSON to file for use by front end
net.write_json_to_file('viz','cluster.json')