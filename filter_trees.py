import sys
from dendropy import TreeList

tc_file = sys.argv[1]
outgroup_file = sys.argv[2]
do_filter = sys.argv[3] == 'filter'
output_file = sys.argv[4]


tl = TreeList.get(path=tc_file, schema='newick')

with open(outgroup_file) as f:
    outgroup_taxa = [x.strip() for x in f.readlines()]

print outgroup_taxa
cnt = 0
bad_trees = [False for _ in tl]

for i, tree in enumerate(tl):
    present_outgroup_taxa = filter(lambda t: tree.find_node_with_taxon_label(t) is not None, outgroup_taxa)
    if len(present_outgroup_taxa):
        node = tree.mrca(taxon_labels=present_outgroup_taxa)
        if node:
            if not (len(node.leaf_nodes()) == len(present_outgroup_taxa) or len(tree.leaf_nodes()) - len(node.leaf_nodes()) == len(present_outgroup_taxa)):
                print 'Tree', i, ' is fxxked up.  There are ', len(node.leaf_nodes()), 'or', len(tree.leaf_nodes()) - len(node.leaf_nodes()), \
                    ' in outgroup, but it supposed to be ', len(present_outgroup_taxa)
                cnt += 1
                bad_trees[i] = True


print cnt, 'out of', len(tl)

if do_filter:
    filtered_tl = TreeList([t for i, t in enumerate(tl) if not bad_trees[i]])
    filtered_tl.write(path=output_file, schema='newick')

