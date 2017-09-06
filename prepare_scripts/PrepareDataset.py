
# coding: utf-8

# In[8]:

import ete3
import os, sys
from pymongo import MongoClient
from bson.objectid import ObjectId
from pymongo.errors import DuplicateKeyError
from pprint import pprint
import json, re, time


# In[9]:

# Import configuration
config_filename = sys.argv[1]
print 'Reading configuration', config_filename, '...'
config = json.load(open(config_filename))


# In[10]:

# Open MongoDB connection
print 'Connecting MongoDB...'
mongo_client = MongoClient('localhost', 27017)
mongo_db = mongo_client.visphy_dev
tree_col = mongo_db.tree
branch_col = mongo_db.branch
input_group_col = mongo_db.inputGroup
entity_col = mongo_db.entity


# In[11]:

process_start_time = time.clock()

# Add input group to MongoDB
print 'Adding dataset...'
all_input_groups = input_group_col.find(projection={'inputGroupId': True})
max_id = 0
for d in all_input_groups:
    if d['inputGroupId'] > max_id:
        max_id = d['inputGroupId']
input_group_id = max_id + 1
print input_group_id

has_missing_taxa = config.get('hasMissingTaxa', True)
input_group_data = {
        'inputGroupId': input_group_id, 
        'title': config['datasetTitle'], 
        'trees': [], 
        'description': config['datasetDescription'],
        'supportRange': config.get("supportRange", [0, 100]),
        'hasMissingTaxa': has_missing_taxa,
        'withParalogs': config.get("paralogs", False)
        }
input_group_col.insert_one(input_group_data)


# In[12]:

# Open tree files
trees = {}
tree_cnt = 0
tree_col.delete_many({'inputGroupId': input_group_id})
default_reference_tree = None

def insert_tree(tree_file, tree_name, tree_type=None, data=None):
    global tree_cnt
    tid = 't' + str(tree_cnt)
    new_tree = {
            'tid': tid,
            'inputGroupId': input_group_id,
            'newickString': data or tree_file.readline().strip(),
            'name': tree_name, 
            'type': tree_type,
            'rfDistance': {},
            'rootBranch': 'b0'
            }
    insert_res = tree_col.insert_one(new_tree)
    trees[tid] = new_tree
    tree_cnt += 1
    return tid
    
print 'Reading tree files and inserting trees to DB...'
if config['mode'] == 1:
    default_reference_tree = insert_tree(open(config['speciesTreePath']), re.search('.*/(\S*).tre$', config['speciesTreePath']).group(1), 'species')     
        
    for gene_num in open(config['genelistPath']):
        tree_path = config['geneTreesPath'].format(gene_num.strip())
        tree_name = re.search('.*/(\S*)/(\S*)$', tree_path).group(1)
        insert_tree(open(tree_path), tree_name, 'gene')
elif config['mode'] == 0:    
    name_regex = re.compile(config.get('treeNameRegex', r'.*'))
    if config['referenceTree']:
        default_reference_tree = insert_tree(open(config['referenceTree']), re.search('.*/(\S*).tre$', config['referenceTree']).group(1), 'species')     

    for dirpath, dirnames, filenames in os.walk(config['treesPath']):
        for filename in filenames:
            if not filename.startswith('.'):
                tree_name = name_regex.search(filename)
                if tree_name:
                    tree_id = insert_tree(open(os.path.join(dirpath, filename)), tree_name.group(), config['treeType'])
                    if default_reference_tree is None:
                            default_reference_tree = tree_id
elif config['mode'] == 2:
    default_reference_tree = insert_tree(open(config['referenceTreePath']), config['referenceTreePath'])
    i = 1
    for line in open(config['treeCollectionPath']):
        insert_tree(None, 'tree_' + str(i), data=line.strip())
        i += 1



# In[13]:

# Put all tree ids to input group document in Mongo
input_group_col.find_one_and_update({'inputGroupId': input_group_id}, {'$set': {'trees': [trees[k]['tid'] for k in trees], 
        'defaultReferenceTree': default_reference_tree}}) 


# In[14]:

# Parse the trees
entities = {}
ent_names = {}
eid_cnt = 0
entity_col.delete_many({'inputGroupId': input_group_id})
branch_col.delete_many({'inputGroupId': input_group_id})

entity_name_regex = None
if 'entityNameRegex' in config:
    entity_name_regex = re.compile(config['entityNameRegex'])


print 'Parsing tree data, inserting entities and branches...'
for tid, tval in trees.iteritems():
    tree_format = config.get('eteNewickFormatReferenceTree', 2) if tid == default_reference_tree else config.get('eteNewickFormat', 2)
    root = None
    print 'Parsing tree', tid
    try:
        root = ete3.Tree(tval['newickString'], format=tree_format)
    except Exception as e:
        print 'fall back to format 0'
        root = ete3.Tree(tval['newickString'], format=0)
    print tid
    # Make sure the tree is strictly bifurcated
    root.resolve_polytomy(recursive=True)
    tval['ete_tree'] = root
    # Process the leaves
    for leaf in root:
        leaf_name = entity_name_regex.match(leaf.name).group() if entity_name_regex is not None else leaf.name
        if leaf_name not in ent_names:
            ent_id = 'e' + str(eid_cnt)
            entity_col.insert_one({'name': leaf_name, 'eid': ent_id, 'inputGroupId': input_group_id})
            ent_names[leaf_name] = ent_id
            eid_cnt += 1
        leaf.add_feature('entity_id', ent_names[leaf_name])
    tree_col.find_one_and_update({'inputGroupId': input_group_id, 'tid': tid}, {'$set': {'entities': [x.entity_id for x in root]}})

    # Process the branches
    branches = []
    root_branch_id = 'b0'
    root.add_feature('branch_id', root_branch_id)
    bid_cnt = 0

    root.ladderize()

    # Assign a bid to each node first
    for node in root.traverse('levelorder'):
        bid = 'b' + str(bid_cnt)
        node.add_feature('branch_id', bid)
        bid_cnt += 1

    # Insert to DB
    for node in root.traverse():
        branch_data = {
            'inputGroupId': input_group_id,
            'bid': node.branch_id, 'tid': tid,
            'length': node.dist, 'support': node.support, 
            'isLeaf': node.is_leaf(), 
            }
        if not node.is_root():
            branch_data['parent'] = node.up.branch_id
        if not node.is_leaf(): 
            branch_data['left'] = node.children[0].branch_id
            branch_data['right'] = node.children[1].branch_id
        else:
            branch_data['entity'] = node.entity_id
        branches.append(branch_data)

    branch_col.insert_many(branches)



# In[15]:

# Calculate RF distance between any two trees
print 'Calculating RF distances...'
for tid1, tval1 in trees.iteritems():
    for tid2, tval2 in trees.iteritems():
        if tid1 < tid2:
            res = tval1['ete_tree'].robinson_foulds(tval2['ete_tree'])
            tval1['rfDistance'][tid2] = res[0]
            tval2['rfDistance'][tid1] = res[0]


# In[16]:

# Update rfDistance in MongoDB
for tid, tval in trees.iteritems():
    tree_col.find_one_and_update({'inputGroupId': input_group_id, 'tid': tid}, {'$set': {'rfDistance': tval['rfDistance']}})
    

# In[17]:

for tid, tval in trees.iteritems():
    for n in tval['ete_tree'].traverse():
        n.add_feature('entity_set', set([x.entity_id for x in n]))


# In[ ]:

# Find corresponding node of "node" from t
#   missing_from_node is the set of entities that are absent in node but present in the other tree
#   missing_from_t is similar
# Return the corresponding node, the jaccard index
def find_corr_node(node, t, missing_from_node=set(), missing_from_t=set(), with_paralogs=False):
    comparing_entities= node.entity_set - missing_from_t
    l = len(comparing_entities)
    if l == 0: 
        return None, 0, None

    # Init
    max_jaccard = 0
    num_entities_diff = 0
    inter = 0
    corr_node = None
    is_ingroup = True
    epsilon = 0.001
    intersection = {}
    intersection_out = {}

    has_no_intersection = lambda x: len(comparing_entities.intersection(x.entity_set - missing_from_node)) == 0
    # Checking ingroup
    for x in t.traverse('postorder', is_leaf_fn=has_no_intersection):
        if has_no_intersection(x):
            intersection[x.branch_id] = 0
            continue
        elif x.is_leaf():
            intersection[x.branch_id] = 1 if x.entity_id not in missing_from_node and x.entity_id in comparing_entities else 0
            # intersection[x.branch_id] = len(comparing_entities.intersection(x.entity_set - missing_from_node)
        elif with_paralogs:
            # Cannot use the memoization technique if there are paralogs
            intersection[x.branch_id] = len(comparing_entities.intersection(x.entity_set - missing_from_node))
        else:
            intersection[x.branch_id] = intersection[x.children[0].branch_id] + intersection[x.children[1].branch_id]

        d = float(intersection[x.branch_id]) / float(l + len(x.entity_set - missing_from_node) - intersection[x.branch_id])
        if d > max_jaccard or (abs(d - max_jaccard) < epsilon and intersection[x.branch_id] > inter):
            max_jaccard = d
            corr_node = x
            inter = intersection[x.branch_id]
            if abs(d - 1.0) < epsilon:
                return corr_node, max_jaccard, is_ingroup

    # Checking outgroup
    full = t.entity_set - missing_from_node
    for x in t.traverse('preorder'):
        if x.is_root():
            intersection_out[x.branch_id] = 0
            continue
        elif with_paralogs:
            intersection_out[x.branch_id] = len(comparing_entities.intersection(full - x.entity_set))
        else:
            sister = x.get_sisters()[0]
            intersection_out[x.branch_id] = intersection_out[x.up.branch_id] + intersection.get(sister.branch_id, 0)

        d = float(intersection_out[x.branch_id]) / float(l + len(full - x.entity_set) - intersection_out[x.branch_id])
        if d > max_jaccard or (abs(d - max_jaccard) < epsilon and intersection_out[x.branch_id] > inter):
            max_jaccard = d
            corr_node = x
            inter = intersection_out[x.branch_id]
            is_ingroup = False
            if abs(d - 1.0) < epsilon:
                return corr_node, max_jaccard, is_ingroup

    return corr_node, max_jaccard, is_ingroup



# In[ ]:

# For every pair of tree, iterate over every node in trees
# i = 0
print 'Getting corresponding branches for the reference tree...'
with_paralogs = True if 'paralogs' in config and config['paralogs'] else False
# Process the default_referent_tree first
# def cmp_func(a, b):
#     if a[0] == default_reference_tree:
#         return -1
#     elif b[0] == default_reference_tree:
#         return 1
#     else:
#         return 0

# for tid1, tval1 in sorted(trees.iteritems(), cmp=cmp_func):
#    print i, tid1,
#    sys.stdout.flush()
#
#    i += 1

def get_corresponding_nodes_for_tree(tid):
    start_time = time.clock()
    t1 = trees[tid]['ete_tree']
    for n in t1.traverse('levelorder'):
        corr_without_missing = {}
        corr = {}
        for tid2, tval2 in trees.iteritems():
            if tid != tid2:
                t2 = tval2['ete_tree']

                # get the missing part
                # the entities missing from tree 1
                missing1 = t2.entity_set - t1.entity_set
                missing2 = t1.entity_set - t2.entity_set
            
                c, max_jaccard, is_ingroup = find_corr_node(n, t2, with_paralogs=with_paralogs)
                corr[tid2] = {
                    'bid': None if c is None else c.branch_id,
                    'jac': max_jaccard,
                    'in': is_ingroup
                }

                if has_missing_taxa:
                    c, max_jaccard, is_ingroup = find_corr_node(n, t2, missing1, missing2, with_paralogs=with_paralogs)
                    corr_without_missing[tid2] = {
                        'bid': None if c is None else c.branch_id,
                        'jac': max_jaccard,
                        'in': is_ingroup
                    }

        # n.add_feature('corr_branches', corr)
        branch_col.find_one_and_update({'inputGroupId': input_group_id, 'tid': tid, 'bid': n.branch_id}, 
                {'$set': {'cb': corr, 'cb2': corr_without_missing}})
    print time.clock() - start_time

get_corresponding_nodes_for_tree(default_reference_tree)
             


# In[ ]:

print 'Time elapsed: ', time.clock() - process_start_time
print 'Done.'

