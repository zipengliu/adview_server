import config
import os
from celery import Celery
from dendropy import Tree, TreeList
from dendropy.calculate import treecompare
from database import Connection
import time

PREPROCESS_STEPS = [
    'Loading data',
    'Re-rooting',
    'Calculating Robinson-Foulds distance',
    'Finding matches for each branch',
    'Building indices for branches and taxa and storing to database',
]

connection = Connection(config.ENVIRONMENT)
connection.test_connection()
app = Celery('preprocess', broker=config.CELERY_BROKER_URL, backend=config.CELERY_RESULT_BACKEND)


def getTid(idx):
    return 't' + str(idx)


def get_corresponding_nodes(node, tree, missing_from_node=set(), missing_from_tree=set()):
    comparing_entities = node.entity_set - missing_from_tree
    numComp = len(comparing_entities)
    if numComp == 0:
        return None, 0
    max_jaccard = 0
    inter = 0
    corr_node = None
    epsilon = 0.001
    intersection = {}

    def has_no_intersection(x):
        return len(comparing_entities.intersection(x.entity_set) - missing_from_node) == 0

    for x in tree.postorder_node_iter():
        if has_no_intersection(x):
            intersection[x.bid] = 0
            continue
        elif x.is_leaf():
            intersection[x.bid] = 1 if x.taxon not in missing_from_node and x.taxon in comparing_entities else 0
        else:
            intersection[x.bid] = sum([intersection[c.bid] for c in x.child_node_iter()], 0)

        d = float(intersection[x.bid]) / float(numComp + len(x.entity_set - missing_from_node) - intersection[x.bid])
        if d > max_jaccard or (abs(d - max_jaccard) < epsilon and intersection[x.bid] > inter):
            max_jaccard = d
            corr_node = x
            inter = intersection[x.bid]
            if abs(d - 1.0) < epsilon:
                return corr_node, max_jaccard

    return corr_node, max_jaccard


@app.task(bind=True)
def preprocess_dataset(self, input_group_id, outgroup_string):
    outgroup_taxa = [x.strip() for x in str(outgroup_string).split(',')]

    # Heavy lifting starts here...
    ###### Read information of this dataset from DB and FS
    self.update_state(state='PROGRESS', meta={'steps': PREPROCESS_STEPS, 'current': 0})

    dataset = connection.input_group.find_one({'inputGroupId': input_group_id})
    dataset_dir = os.path.join(config.DATA_FILES_PATH, str(input_group_id))
    reference_tree = Tree.get(path=os.path.join(dataset_dir, dataset['referenceTreeFileName']), schema='newick')
    reference_tree.tid = getTid(0)
    tn = reference_tree.taxon_namespace
    tree_collection = TreeList.get(path=os.path.join(dataset_dir, config.TREE_COLLECTION_FILENAME), schema='newick', taxon_namespace=tn)
    for i, tree in enumerate(tree_collection):
        tree.tid = getTid(i + 1)
    tree_collection_names_path = os.path.join(dataset_dir, config.TREE_COLLECTION_NAMES_FILENAME)

    tree_collection_names = [l.strip() for l in open(tree_collection_names_path)] \
        if os.path.exists(tree_collection_names_path) else \
        ['tree_' + str(i) for i in xrange(1, len(tree_collection) + 1)]


    ####### Re-root if needed
    self.update_state(state='PROGRESS', meta={'steps': PREPROCESS_STEPS, 'current': 1})

    def reroot_by_outgroup(tree, outgroup_taxa):
        node = tree.mrca(taxon_labels=outgroup_taxa)
        if node:
            # tree.to_outgroup_posiiton(node, update_bipartitions=True)
            tree.reroot_at_node(node, update_bipartitions=True)

    if len(outgroup_taxa):
        if not dataset['isReferenceRooted']:
            reroot_by_outgroup(reference_tree, outgroup_taxa)
        if not dataset['isTCRooted']:
            for tree in tree_collection:
                reroot_by_outgroup(tree, outgroup_taxa)


    ######### Tree distances
    total_pairs = (len(tree_collection) + 1) * len(tree_collection) / 2
    self.update_state(state='PROGRESS', meta={'steps': PREPROCESS_STEPS, 'current': 2, 'progress': {'done': 0, 'total': total_pairs}})
    # Init
    done_pairs = 0
    reference_tree.distances = {}
    for t in tree_collection:
        t.distances = {}

    for i, t1 in enumerate(tree_collection):
        d = treecompare.symmetric_difference(reference_tree, t1)
        reference_tree.distances[t1.tid] = d
        t1.distances[reference_tree.tid] = d
        for j in xrange(i + 1, len(tree_collection)):
            t2 = tree_collection[j]
            d = treecompare.symmetric_difference(t1, t2)
            t1.distances[t2.tid] = d
            t2.distances[t1.tid] = d
        done_pairs += len(tree_collection) - i
        self.update_state(state='PROGRESS', meta={'steps': PREPROCESS_STEPS, 'current': 2,
                                                  'progress': {'done': done_pairs, 'total': total_pairs}})


    ####### Find matches
    total = sum([1 for _ in reference_tree], 0)
    self.update_state(state='PROGRESS', meta={'steps': PREPROCESS_STEPS, 'current': 3,
                                              'progress': {'done': 0, 'total': total}})
    # Init
    done = 0
    all_entities = set(tn)
    reference_tree.entities = set([node.taxon for node in reference_tree.leaf_node_iter()])
    for tree in tree_collection:
        tree.entities = set([node.taxon for node in tree.leaf_node_iter()])

    # Assign branch id to each branch
    def assign_bid(tree):
        i = 0
        for node in tree.levelorder_node_iter():
            bid = 'b' + str(i)
            i += 1
            node.bid = bid
    assign_bid(reference_tree)
    for tree in tree_collection:
        assign_bid(tree)

    # Cache the entity set for each branch in tree collection (for performance reason)
    for node in reference_tree.levelorder_node_iter():
        node.entity_set = set([l.taxon for l in node.leaf_iter()])
    for tree in tree_collection:
        for node in tree.levelorder_node_iter():
            node.entity_set = set([l.taxon for l in node.leaf_iter()])

    for node in reference_tree:
        node.corresponding_nodes = {}
        for tree in tree_collection:
            if node.is_leaf():
                matched_taxon_node = tree.find_node_for_taxon(node.taxon)
                node.corresponding_nodes[tree.tid] = {'jac': 1, 'bid': matched_taxon_node.bid if matched_taxon_node else None}
            elif node.parent_node is None:
                inter = len(reference_tree.entities.intersection(tree.entities))
                node.corresponding_nodes[tree.tid] = {
                    'jac': float(inter) / (len(reference_tree.entities) + len(tree.entities) - inter),
                    'bid': tree.seed_node.bid
                }
            else:
                matched_node, max_jaccard = \
                    get_corresponding_nodes(node, tree, missing_from_node=all_entities - reference_tree.entities,
                                            missing_from_tree=all_entities - tree.entities)
                node.corresponding_nodes[tree.tid] = {
                    'jac': max_jaccard,
                    'bid': matched_node.bid if matched_node is not None else None
                }
        done += 1
        self.update_state(state='PROGRESS', meta={'steps': PREPROCESS_STEPS, 'current': 3,
                                                  'progress': {'done': done, 'total': total}})


    ######## Build index
    self.update_state(state='PROGRESS', meta={'steps': PREPROCESS_STEPS, 'current': 4})

    # all taxa
    connection.entity.insert_many([{
        'name': taxon.label,
        'eid': 'e' + str(i),
        'inputGroupId': input_group_id}
        for i, taxon in enumerate(tn)])

    def insert_tree(tree, name):
        data = {
            'tid': tree.tid,
            'inputGroupId': input_group_id,
            'newickString': tree.as_string('newick'),
            'name': name,
            'rfDistance': tree.distances,
            'rootBranch': 'b0',
            'entities': ['e' + str(tn.accession_index(leaf.taxon)) for leaf in tree.leaf_node_iter()]
        }
        branches = []
        for node in tree.levelorder_node_iter():
            d = {
                'bid': node.bid,
                'tid': tree.tid,
                'inputGroupId': input_group_id,
                'length': node.edge_length or 0
            }
            if node.is_leaf():
                d['entity'] = 'e' + str(tn.accession_index(node.taxon))
            else:
                support_val = 0 if not node.label else float(node.label)
                cn = node.child_nodes()
                d.update({
                    'left': cn[0].bid,
                    'right': cn[1].bid,
                    'support': support_val,
                })
                # Multifurcation is likely to happen at the root node because of the re-rooting
                if len(cn) == 3:
                    d['rightmost'] = cn[2].bid
            if tree == reference_tree:
                d['cb'] = node.corresponding_nodes
            branches.append(d)
        data['branches'] = branches

        # The data for a tree might be large so we do not use bulk insert
        connection.tree.insert_one(data)
        connection.branch.insert_many(branches)

    reference_tree.ladderize()
    insert_tree(reference_tree, dataset['referenceTreeFileName'])

    for tno, tree in enumerate(tree_collection):
        tree.ladderize()
        insert_tree(tree, tree_collection_names[tno])

    return {'input_group_id': input_group_id}





