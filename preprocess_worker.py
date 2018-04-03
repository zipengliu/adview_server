import config
import os
from celery import Celery
from dendropy import Tree, TreeList, Node
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

connection = Connection(config.ENVIRONMENT, connect=False)
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
def preprocess_dataset(self, input_group_id, outgroup_string, is_updating=False):
    outgroup_taxa = filter(lambda x: len(x) > 0, [x.strip() for x in str(outgroup_string).split(',')])
    print 'Outgroup taxa:', outgroup_taxa

    # Heavy lifting starts here...
    ###### Read information of this dataset from DB and FS
    self.update_state(state='PROGRESS', meta={'steps': PREPROCESS_STEPS, 'current': 0})

    dataset = connection.input_group.find_one({'inputGroupId': input_group_id})
    dataset_dir = os.path.join(config.DATA_FILES_PATH, str(input_group_id))
    reference_tree = Tree.get(path=os.path.join(dataset_dir, dataset['referenceTreeFileName']), schema='newick',
                              rooting='force-rooted' if dataset['isReferenceRooted'] else 'force-unrooted')
    reference_tree.tid = getTid(0)
    tn = reference_tree.taxon_namespace
    tree_collection = TreeList.get(path=os.path.join(dataset_dir, config.TREE_COLLECTION_FILENAME), schema='newick',
                                   taxon_namespace=tn,
                                   rooting='force-rooted' if dataset['isTCRooted'] else 'force-unrooted')
    for i, tree in enumerate(tree_collection):
        tree.tid = getTid(i + 1)
    tree_collection_names_path = os.path.join(dataset_dir, config.TREE_COLLECTION_NAMES_FILENAME)

    tree_collection_names = [l.strip() for l in open(tree_collection_names_path)] \
        if os.path.exists(tree_collection_names_path) else \
        ['tree_' + str(i) for i in xrange(1, len(tree_collection) + 1)]


    ####### Re-root if needed
    self.update_state(state='PROGRESS', meta={'steps': PREPROCESS_STEPS, 'current': 1})

    # A special treatment for the root to resolve a tri-nary root branch
    # The outgroup (first child) must be singled out
    def resolve_polytomy_at_root(tree, outgroup_taxa):
        root = tree.seed_node
        if root.num_child_nodes() == 2:
            return
        if root.num_child_nodes() > 3:
            tree.resolve_polytomies()
            tree.print_plot()
        assert root.num_child_nodes() == 3
        sorted_children = sorted(root.child_nodes(), key=lambda t: len(t.leaf_nodes()))

        middle_node = sorted_children[1]
        found = False
        for x in middle_node.leaf_iter():
            if x in outgroup_taxa:
                found = True
                break
        c1 = None
        c2 = None
        if found:
            # merge the first two
            c1 = root._child_nodes[0]
            c2 = root._child_nodes[1]
        else:
            # merge the 2nd and 3rd
            c1 = root._child_nodes[1]
            c2 = root._child_nodes[2]
        root.remove_child(c1)
        root.remove_child(c2)
        nn1 = Node()
        nn1.edge.length = 0.0
        nn1.add_child(c1)
        nn1.add_child(c2)
        root.add_child(nn1)

        # n = node.num_child_nodes()
        # i = 0
        # while n > 2:
        #     nn1 = Node()
        #     nn1.edge.length = 0.0
        #     c1 = sorted_children[i]
        #     c2 = sorted_children[i + 1]
        #     # c1 = node._child_nodes[n - 2]
        #     # c2 = node._child_nodes[n - 1]
        #     node.remove_child(c1)
        #     node.remove_child(c2)
        #     nn1.add_child(c1)
        #     nn1.add_child(c2)
        #     node.add_child(nn1)
        #     n -= 2
        #     i += 2

    def reroot_by_outgroup(tree, outgroup_taxa):
        present_outgroup_taxa = filter(lambda t: tree.find_node_with_taxon_label(t) is not None, outgroup_taxa)
        if len(present_outgroup_taxa):
            tree.resolve_polytomies()
            node = tree.mrca(taxon_labels=present_outgroup_taxa)
            if node and node.parent_node:
                # tree.reroot_at_edge(node.edge)
                tree.to_outgroup_position(node, update_bipartitions=False)
                # tree.reroot_at_node(node, update_bipartitions=True)

            # This is not what I want: it is going to group the outgroup and a part of the ingroup together first
            # because it joins pairs of children in the order given
            # tree.resolve_polytomies(update_bipartitions=True)

            tree.is_rooted = True
            resolve_polytomy_at_root(tree, present_outgroup_taxa)

            # Turns out I have to call this function because it needs to resolve polytomies deep inside the tree,
            #   not just the root!
            tree.resolve_polytomies()
            tree.update_bipartitions()
        else:
            # All outgroup taxa is missing
            # This should rarely happen!
            print 'Tree with no outgroup taxa'

            tree.reroot_at_midpoint()

            assert tree.seed_node is not None
            tree.is_rooted = True
            tree.resolve_polytomies(update_bipartitions=True)

        assert tree.seed_node.num_child_nodes() < 3

    if len(outgroup_taxa):
        if not dataset['isReferenceRooted']:
            reroot_by_outgroup(reference_tree, outgroup_taxa)
        if not dataset['isTCRooted']:
            for tree in tree_collection:
                reroot_by_outgroup(tree, outgroup_taxa)
    else:
        assert dataset['isReferenceRooted'], 'Cannot deal with unrooted trees without an outgroup'
        assert dataset['isTCRooted'], 'Cannot deal with unrooted trees without an outgroup'

    if dataset['isReferenceRooted']:
        reference_tree.resolve_polytomies(update_bipartitions=True)

    if dataset['isTCRooted']:
        for tree in tree_collection:
            tree.resolve_polytomies(update_bipartitions=True)

    reference_tree.ladderize()
    for tree in tree_collection:
        tree.ladderize()


    ######### Tree distances
    total_pairs = (len(tree_collection) + 1) * len(tree_collection) / 2
    self.update_state(state='PROGRESS', meta={'steps': PREPROCESS_STEPS, 'current': 2, 'progress': {'done': 0, 'total': total_pairs}})
    # Init
    done_pairs = 0
    reference_tree.distances = {}
    for t in tree_collection:
        t.distances = {}

    if len(tree_collection) > config.TREE_DISTANCE_THRESHOLD:
        print '#trees exceed limit.  Do not calcualte RF distances'
    else:
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
    print 'Start finding corresponding branches'

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
    print 'Start stroing restuls to DB'
    self.update_state(state='PROGRESS', meta={'steps': PREPROCESS_STEPS, 'current': 4})

    # If we need to update the outgroup, remove that dataset first
    connection.entity.remove({'inputGroupId': input_group_id})
    connection.branch.remove({'inputGroupId': input_group_id})
    connection.tree.remove({'inputGroupId': input_group_id})

    # all taxa
    connection.entity.insert_many([{
        'name': taxon.label,
        'eid': 'e' + str(i),
        'inputGroupId': input_group_id}
        for i, taxon in enumerate(tn)])

    def find_outgroup_branch_bid(tree, outgroup_taxa):
        # avoid non-monophyletic outgroup
        # assume the the first child is the outgroup (all taxa under ourgroup branch must be outgroup taxa)
        present_outgroup_taxa = filter(lambda t: tree.find_node_with_taxon_label(t) is not None, outgroup_taxa)
        node = tree.seed_node.child_nodes()[0]
        if len(node.leaf_nodes()) > len(present_outgroup_taxa):
            return None
        for x in node.leaf_iter():
            if x not in present_outgroup_taxa:
                return None
        return node.bid

    def insert_tree(tree, name):
        data = {
            'tid': tree.tid,
            'inputGroupId': input_group_id,
            'newickString': tree.as_string('newick'),
            'name': name,
            'rfDistance': tree.distances,
            'rootBranch': 'b0',
            'entities': ['e' + str(tn.accession_index(leaf.taxon)) for leaf in tree.leaf_node_iter()],
        }
        if len(outgroup_taxa):
            # node = tree.mrca(taxon_labels=present_outgroup_taxa) if len(present_outgroup_taxa) else None
            # assert(node.parent_node == tree.seed_node)
            data['outgroupBranch'] = find_outgroup_branch_bid(tree, outgroup_taxa)

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
                if not node.label:
                    support_val = 0
                    label = ''
                else:
                    try:
                        support_val = float(node.label)
                        label = ''
                    except ValueError:
                        support_val = 0
                        label = node.label
                cn = node.child_nodes()
                d.update({
                    'left': cn[0].bid,
                    'right': cn[1].bid,
                    'support': support_val,
                })
                if len(label):
                    d['label'] = label

                # Multifurcation is likely to happen at the root node because of the re-rooting
                if len(cn) == 3:
                    tree.print_plot()
                    print node
                assert len(cn) < 3, 'Multifurcation at {} in tree {} (root= {})'.format(node.bid, tree.tid, tree.seed_node.bid)
                # if len(cn) == 3:
                #     d['rightmost'] = cn[2].bid
            if tree == reference_tree:
                d['cb'] = node.corresponding_nodes
            branches.append(d)

        connection.insert_tree(data)
        connection.insert_branches(branches)

    insert_tree(reference_tree, dataset['referenceTreeFileName'])

    for tno, tree in enumerate(tree_collection):
        insert_tree(tree, tree_collection_names[tno])

    connection.input_group.find_one_and_update({'inputGroupId': input_group_id},
                                               {'$set': {'outgroupTaxa': ['e' + str(tn.accession_index(tn.get_taxon(t)))
                                                                          for t in outgroup_taxa]}})
    print 'DONE'

    return {'input_group_id': input_group_id}





